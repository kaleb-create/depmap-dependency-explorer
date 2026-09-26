import csv
import json
import os
import sqlite3
import tempfile
import threading
import unittest
from pathlib import Path
from contextlib import closing
from unittest.mock import patch

import dependency_stratifiers as ds
from scripts.build_hpv_dependency_data import build_cohort_catalog
from stratifier_jobs import StratifierJobs


def group(**fields):
    return dict(lineages=[], diseases=[], subtypes=[], include_terms=[], exclude_terms=[],
                annotations=[], annotation_logic="all") | fields


def model(index, lineage="Pancreas", disease="Pancreatic Adenocarcinoma"):
    return {"ModelID": f"M{index}", "CellLineName": f"Cell{index}", "OncotreeLineage": lineage,
            "OncotreePrimaryDisease": disease, "OncotreeSubtype": "", "CCLEName": f"CELL{index}",
            "StrippedCellLineName": f"CELL{index}", "RRID": f"RRID{index}"}


def plan(**fields):
    return dict(source_kind="metadata", negative_mode="complement", unsupported_reason="",
                positive=group(lineages=["Pancreas"]), negative=group(), dataset={"format": "metadata_only"}) | fields


class CohortTests(unittest.TestCase):
    def setUp(self):
        self.models = [model(i) for i in range(3)] + [model(i, "Lung", "Lung Adenocarcinoma") for i in range(3, 6)]
        self.models += [model(6, "Skin", "Non-Cancerous"), model(7, "", "")]

    def test_metadata_complement_excludes_non_cancer_and_unknown(self):
        spec = plan()
        positive, negative, *_ = ds.resolve_cohorts(self.models, spec, {})
        self.assertEqual(positive, {"M0", "M1", "M2"})
        self.assertEqual(negative, {"M3", "M4", "M5"})
        self.assertEqual(spec["dataset"]["name"], "DepMap Model.csv")

    def test_catalog_export_requires_source_records_for_negatives(self):
        definition = dict(id="kras", category="Top driver mutations by prevalence",
                          label="KRAS", positive_label="Mutant", negative_label="Non-mutant",
                          source="DepMap driver mutations", positive_model_ids={"M0"})
        copies = {gene: {"M0", "M1", "M2"} for gene in ("PTEN", "RB1", "SMAD4", "NF1", "KEAP1", "MDM2", "BRCA1", "BRCA2")}
        catalog = build_cohort_catalog([definition], {"M0", "M1", "M2"}, {"M0", "M1"}, set(), copies)
        self.assertEqual(catalog[0]["negative_ids"], ["M1"])
        self.assertTrue(catalog[0]["limitations"])

    def test_metadata_fields_are_intersection(self):
        self.assertFalse(ds.model_matches(self.models[0], group(lineages=["Pancreas"], diseases=["Lung Adenocarcinoma"])))

    def test_categorical_metadata_is_exact_and_missing_is_excluded(self):
        for i, row in enumerate(self.models[:6]):
            row["Sex"] = "Female" if i < 3 else "Male"
        spec = plan(positive=group(attributes=[{"field": "Sex", "values": ["Female"]}]),
                    negative=group(attributes=[{"field": "Sex", "values": ["Male"]}]), negative_mode="explicit")
        positive, negative, *_ = ds.resolve_cohorts(self.models + [model(9)], spec, {})
        self.assertEqual(positive, {"M0", "M1", "M2"})
        self.assertEqual(negative, {"M3", "M4", "M5"})
        with self.assertRaisesRegex(ValueError, "attribute"):
            ds.metadata_scope(self.models, group(attributes=[{"field": "Sex", "values": ["guessed"]}]))

    def test_unknown_or_empty_definition_does_not_broaden(self):
        for positive in (group(lineages=["Typo disease"]), group()):
            with self.assertRaises(ValueError):
                ds.resolve_cohorts(self.models, plan(positive=positive), {})

    def test_catalog_complement_uses_known_negatives(self):
        catalog = {"loss": {"positive_ids": ["M0", "M1", "M2"], "negative_ids": ["M3", "M4", "M5"], "source": "Known calls"}}
        comparison = plan(source_kind="catalog", positive=group(annotations=[{"analysis_id": "loss", "side": "positive"}]))
        models = self.models + [model(8)]
        positive, negative, *_ = ds.resolve_cohorts(models, comparison, catalog)
        self.assertNotIn("M8", positive | negative)
        self.assertEqual(negative, {"M3", "M4", "M5"})

    def test_composite_preserves_unknown_status(self):
        catalog = {
            "a": {"positive_ids": ["M0", "M1"], "negative_ids": ["M2", "M3"]},
            "b": {"positive_ids": ["M0", "M2"], "negative_ids": ["M1"]},
        }
        rule = group(annotations=[{"analysis_id": x, "side": "positive"} for x in ("a", "b")])
        selected, known = ds.annotation_scope(rule, catalog, {"M0", "M1", "M2", "M3"})
        self.assertEqual(selected, {"M0"})
        self.assertEqual(known, {"M0", "M1", "M2"})

    def test_overlap_fails(self):
        with self.assertRaisesRegex(ValueError, "overlap"):
            ds.resolve_cohorts(self.models, plan(negative_mode="explicit"), {})

    def test_external_absence_is_not_negative_and_conflicts_excluded(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "source.csv"
            path.write_text("ModelID,Status\nM0,yes\nM1,no\nM2,unknown\nM3,\nM4,yes\nM4,no\n")
            positive, negative, mapping = ds.map_external_dataset(self.models,
                {"identifier_column": "ModelID", "group_column": "Status", "positive_values": ["yes"], "negative_values": ["no"]},
                {"cached_path": str(path), "format": "csv"})
            self.assertEqual(positive, {"M0"})
            self.assertEqual(negative, {"M1"})
            self.assertEqual(mapping["conflicting_models_excluded"], 1)
            positive, negative, _ = ds.map_external_dataset(self.models,
                {"group_column": "Status", "positive_values": ["yes"], "negative_values": []},
                {"cached_path": str(path), "format": "csv"})
            self.assertEqual(negative, {"M1"})
            self.assertEqual(positive, {"M0"})

    def test_numeric_json_zero_is_a_valid_class(self):
        self.assertEqual(ds.normalized(0), "0")
        self.assertEqual(ds.normalized(None), "")

    def test_ambiguous_alias_not_assigned_arbitrarily(self):
        self.models[1]["CellLineName"] = self.models[0]["CellLineName"]
        lookup = ds.model_identifier_lookup(self.models)
        self.assertNotIn("cell0", lookup)
        self.assertEqual(lookup["m0"], "M0")

    def test_nonfinite_measurements_are_not_counted(self):
        self.assertEqual(ds.values_at(["nan", "inf", "-inf", "-1.2"], range(4)), [-1.2])
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "matrix.csv"
            path.write_text("ModelID,GOOD,BAD\nM0,-1,nan\nM1,-1.1,inf\nM2,-.9,-1\nM3,0,.1\nM4,.1,.2\nM5,-.1,0\n")
            rows, _ = ds.row_pair_differential(str(path), {"M0", "M1", "M2"}, {"M3", "M4", "M5"})
            self.assertEqual([r["gene"] for r in rows], ["GOOD"])


class JobTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.dbpath = str(Path(self.temp.name) / "jobs.db")
        self.release = threading.Event()
        self.started = threading.Event()
        with closing(self.connect()) as db, db:
            db.execute("CREATE TABLE results (id INTEGER PRIMARY KEY, prompt TEXT)")

    def tearDown(self):
        self.release.set()
        if hasattr(self, "jobs") and self.jobs.thread:
            self.jobs.thread.join(5)
        self.temp.cleanup()

    def connect(self):
        conn = sqlite3.connect(self.dbpath, timeout=10)
        conn.row_factory = sqlite3.Row
        return conn

    def save(self, db, prompt, *_):
        return db.execute("INSERT INTO results(prompt) VALUES (?)", (prompt,)).lastrowid

    def compute(self, prompt, progress):
        progress("Computing CRISPR")
        self.started.set()
        self.release.wait(5)
        if prompt == "fail":
            raise ValueError("No classified negatives")
        return {}, {}, {}

    def make_jobs(self):
        self.jobs = StratifierJobs(self.connect, lambda db, sql, args=(): db.execute(sql, args),
                                  self.compute, self.save, self.temp.name)

    def test_background_status_and_duplicate_submission(self):
        self.make_jobs()
        job = self.jobs.enqueue("Pancreas", "same-request")
        self.assertTrue(self.started.wait(2))
        self.assertEqual(self.jobs.get(job["id"])["stage"], "Computing CRISPR")
        again = self.jobs.enqueue("Pancreas", "same-request")
        self.assertEqual(job["id"], again["id"])
        self.release.set()
        self.jobs.thread.join(5)
        self.assertEqual(self.jobs.get(job["id"])["status"], "complete")
        with closing(self.connect()) as db, db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM results").fetchone()[0], 1)

    def test_failure_keeps_prompt_and_message(self):
        self.make_jobs()
        self.release.set()
        with self.assertLogs("stratifier_jobs", level="ERROR"):
            job = self.jobs.enqueue("fail", "failure")
            self.jobs.thread.join(5)
        result = self.jobs.get(job["id"])
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["error"], "No classified negatives")
        self.assertEqual(result["prompt"], "fail")

    def test_interrupted_job_resumes(self):
        self.make_jobs()
        with patch.object(self.jobs, "start"):
            job = self.jobs.enqueue("Pancreas", "resume")
        with closing(self.connect()) as db, db:
            db.execute("UPDATE stratifier_jobs SET status='running', attempts=1")
        self.release.set()
        self.jobs.start()
        self.jobs.thread.join(5)
        result = self.jobs.get(job["id"])
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["attempts"], 2)


if __name__ == "__main__":
    unittest.main()
