import csv
import json
import math
import os
import ssl
import urllib.parse
import urllib.request
from datetime import datetime, timezone


BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE_DIR, "data")
DATA_DIR = os.environ.get("DEPMAP_DATA_DIR", DATA_DIR)
STATIC_DATA_DIR = os.environ.get("DEPMAP_STATIC_DATA_DIR", os.path.join(BASE_DIR, "static", "data"))
FILE_INDEX = os.path.join(DATA_DIR, "depmap_files.csv")
SUMMARY_JSON = os.path.join(STATIC_DATA_DIR, "hpv_dependency_summary.json")
ANALYSIS_DATA_DIR = os.path.join(STATIC_DATA_DIR, "dependency_analyses")

DOWNLOAD_INDEX_URL = "https://depmap.org/portal/api/download/files"
CELLOSAURUS_SEARCH_URL = "https://api.cellosaurus.org/search/cell-line"
CRISPR_RELEASE = "DepMap Public 26Q1"
RNAI_RELEASE = "DEMETER2 Data v6"
HPV_TRANSFORMANT_QUERY = "transformant:papillomavirus"
MTAP_COPY_NUMBER_URL = "https://ndownloader.figshare.com/files/34008428"
MTAP_COPY_NUMBER_FILENAME = "CCLE_gene_cn_22Q1.csv"
GLOBAL_SIGNATURE_URL = "https://ndownloader.figshare.com/files/46500361"
GLOBAL_SIGNATURE_FILENAME = "OmicsSignatures_24Q2.csv"
EXPRESSION_URL = "https://ndownloader.figshare.com/files/46490878"
EXPRESSION_FILENAME = "OmicsExpressionProteinCodingGenesTPMLogp1_24Q2.csv"
MTAP_DELETION_THRESHOLD = 0.4
DEEP_DELETION_THRESHOLD = 0.3
AMPLIFICATION_THRESHOLD = 2.0
TOP_DRIVER_COUNT = 20
MIN_GROUP_VALUES = 3
MIN_GROUP_COVERAGE = 0.5

AMPLIFICATION_GENES = [
    "MYC",
    "CCNE1",
    "CCND1",
    "ERBB2",
    "EGFR",
    "MDM2",
    "MET",
    "MYCN",
    "CDK4",
    "FGFR1",
    "FGFR2",
    "FGFR3",
    "SOX2",
    "TP63",
]
DELETION_GENES = ["CDKN2A", "CDKN2B", "MTAP", "PTEN", "RB1", "SMAD4", "NF1", "BRCA1", "BRCA2", "KEAP1"]
EPITHELIAL_MARKERS = ["CDH1", "EPCAM", "KRT8", "KRT18", "KRT19", "OCLN"]
MESENCHYMAL_MARKERS = ["VIM", "CDH2", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1"]
NEUROENDOCRINE_MARKERS = ["ASCL1", "NEUROD1", "INSM1", "CHGA", "SYP", "NCAM1"]

ADDITIONAL_MUTATION_ANALYSES = [
    ("egfr", "EGFR mutant vs non-mutant", ["EGFR"]),
    ("pi3k", "PI3K mutant vs non-mutant", ["PIK3CA", "PIK3R1"]),
    ("hras", "HRAS mutant vs non-mutant", ["HRAS"]),
    ("ctnnb1", "CTNNB1 mutant vs non-mutant", ["CTNNB1"]),
    ("brca", "BRCA1/2 mutant vs non-mutant", ["BRCA1", "BRCA2"]),
]


def ensure_dirs() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(STATIC_DATA_DIR, exist_ok=True)
    os.makedirs(ANALYSIS_DATA_DIR, exist_ok=True)


def download(url: str, path: str) -> None:
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return
    ctx = ssl._create_unverified_context()
    with urllib.request.urlopen(url, context=ctx) as src, open(path, "wb") as dst:
        while True:
            chunk = src.read(1024 * 1024)
            if not chunk:
                break
            dst.write(chunk)


def read_index() -> list[dict[str, str]]:
    download(DOWNLOAD_INDEX_URL, FILE_INDEX)
    with open(FILE_INDEX, newline="") as f:
        return list(csv.DictReader(f))


def find_file(index: list[dict[str, str]], release: str, filename: str) -> dict[str, str]:
    for row in index:
        if row["release"] == release and row["filename"] == filename:
            return row
    raise RuntimeError(f"Missing DepMap file {filename} in {release}")


def parse_gene(label: str) -> tuple[str, str | None]:
    label = label.strip()
    if label.endswith(")") and " (" in label:
        gene, entrez = label.rsplit(" (", 1)
        return gene, entrez[:-1]
    return label, None


def model_lookup(model_path: str) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    with open(model_path, newline="") as f:
        for row in csv.DictReader(f):
            lookup[row["ModelID"]] = {
                "model_id": row["ModelID"],
                "cell_line": row["CellLineName"] or row["StrippedCellLineName"],
                "stripped_name": row["StrippedCellLineName"],
                "ccle_name": row["CCLEName"],
                "cellosaurus": row["RRID"],
                "lineage": row["OncotreeLineage"],
                "disease": row["OncotreePrimaryDisease"],
            }
    return lookup


def fetch_hpv_transformants() -> dict[str, list[str]]:
    ctx = ssl._create_unverified_context()
    accessions: dict[str, list[str]] = {}
    for start in range(0, 10000, 500):
        params = urllib.parse.urlencode(
            {
                "q": HPV_TRANSFORMANT_QUERY,
                "fields": "id,ac,transformant",
                "format": "json",
                "rows": "500",
                "start": str(start),
            }
        )
        with urllib.request.urlopen(f"{CELLOSAURUS_SEARCH_URL}?{params}", context=ctx) as src:
            batch = json.load(src)["Cellosaurus"]["cell-line-list"]
        if not batch:
            break
        for cell_line in batch:
            labels = [
                transformant.get("xref", {}).get("label", "")
                for transformant in cell_line.get("transformant-list", [])
                if "papillomavirus" in transformant.get("xref", {}).get("label", "").lower()
            ]
            if not labels:
                continue
            primary_accession = next(
                (a["value"] for a in cell_line.get("accession-list", []) if a.get("type") == "primary"),
                None,
            )
            if primary_accession:
                accessions[primary_accession] = labels
    return accessions


def is_driver_mutation(row: dict[str, str]) -> bool:
    return (
        row.get("IsDefaultEntryForModel") == "Yes"
        and (
            row.get("Hotspot") == "True"
            or row.get("HessDriver") in {"True", "Y"}
            or row.get("OncogeneHighImpact") == "True"
            or row.get("TumorSuppressorHighImpact") == "True"
        )
    )


def mutation_groups(
    mutation_path: str,
) -> tuple[dict[str, set[str]], list[str], dict[str, set[str]], dict[str, set[str]]]:
    by_gene: dict[str, set[str]] = {}
    likely_biallelic_lof: dict[str, set[str]] = {}
    with open(mutation_path, newline="") as f:
        for row in csv.DictReader(f):
            gene = row.get("HugoSymbol")
            if not gene or row.get("IsDefaultEntryForModel") != "Yes":
                continue
            if is_driver_mutation(row):
                by_gene.setdefault(gene, set()).add(row["ModelID"])
            try:
                allele_fraction = float(row.get("AF") or 0)
            except ValueError:
                allele_fraction = 0
            if row.get("LikelyLoF") == "True" and allele_fraction >= 0.5:
                likely_biallelic_lof.setdefault(gene, set()).add(row["ModelID"])
    top_driver_genes = [
        gene
        for gene, _ in sorted(
            by_gene.items(),
            key=lambda item: (-len(item[1]), item[0]),
        )[:TOP_DRIVER_COUNT]
    ]
    groups = {
        key: set().union(*(by_gene.get(gene, set()) for gene in genes))
        for key, _, genes in ADDITIONAL_MUTATION_ANALYSES
    }
    groups.update({gene.lower(): by_gene[gene] for gene in top_driver_genes})
    groups["alk"] = by_gene.get("ALK", set())
    return groups, top_driver_genes, by_gene, likely_biallelic_lof


def copy_number_groups(
    copy_number_path: str,
) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, dict[str, float]]]:
    selected_genes = sorted(set(AMPLIFICATION_GENES + DELETION_GENES))
    groups = {f"{gene.lower()}_amp": set() for gene in AMPLIFICATION_GENES}
    groups.update({f"{gene.lower()}_del": set() for gene in DELETION_GENES})
    measured = {gene: set() for gene in selected_genes}
    values_by_model: dict[str, dict[str, float]] = {}
    with open(copy_number_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        indexes = {
            parse_gene(label)[0]: i
            for i, label in enumerate(header)
            if parse_gene(label)[0] in selected_genes
        }
        missing = set(selected_genes) - set(indexes)
        if missing:
            raise RuntimeError(
                f"Missing copy-number genes in {os.path.basename(copy_number_path)}: {', '.join(sorted(missing))}"
            )
        for row in reader:
            model_id = row[0]
            model_values: dict[str, float] = {}
            for gene, index in indexes.items():
                if index >= len(row) or row[index] in {"", "NA"}:
                    continue
                value = float(row[index])
                model_values[gene] = value
                measured[gene].add(model_id)
                if gene in AMPLIFICATION_GENES and value >= AMPLIFICATION_THRESHOLD:
                    groups[f"{gene.lower()}_amp"].add(model_id)
                if gene in DELETION_GENES and value < DEEP_DELETION_THRESHOLD:
                    groups[f"{gene.lower()}_del"].add(model_id)
            values_by_model[model_id] = model_values

    groups["mtap_del"] = {
        model_id
        for model_id, values in values_by_model.items()
        if values.get("MTAP", math.inf) < MTAP_DELETION_THRESHOLD
    }
    groups["cdkn2ab_del"] = groups["cdkn2a_del"] | groups["cdkn2b_del"]
    groups["mtap_cdkn2a_codeletion"] = {
        model_id
        for model_id, values in values_by_model.items()
        if values.get("MTAP", math.inf) < MTAP_DELETION_THRESHOLD
        and values.get("CDKN2A", math.inf) < MTAP_DELETION_THRESHOLD
    }
    groups["cdkn2a_only_deletion"] = {
        model_id
        for model_id, values in values_by_model.items()
        if values.get("CDKN2A", math.inf) < MTAP_DELETION_THRESHOLD
        and values.get("MTAP", -math.inf) >= MTAP_DELETION_THRESHOLD
    }
    return groups, measured, values_by_model


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("Cannot calculate a quantile from an empty list")
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def global_signature_groups(
    signature_path: str,
) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, dict[str, object]], dict[str, float]]:
    rows: dict[str, dict[str, float]] = {}
    with open(signature_path, newline="") as f:
        for row in csv.DictReader(f):
            model_id = row.get("") or row.get("ModelID")
            if not model_id:
                continue
            parsed = {}
            for column in ["MSIScore", "Ploidy", "WGD", "LoHFraction", "CIN", "Aneuploidy"]:
                try:
                    parsed[column] = float(row[column])
                except (KeyError, TypeError, ValueError):
                    pass
            rows[model_id] = parsed

    thresholds = {}
    for column, key in [("CIN", "cin"), ("Aneuploidy", "aneuploidy"), ("LoHFraction", "loh")]:
        values = [row[column] for row in rows.values() if column in row]
        thresholds[f"{key}_low"] = quantile(values, 0.25)
        thresholds[f"{key}_high"] = quantile(values, 0.75)

    groups = {
        "msi_high": {m for m, row in rows.items() if row.get("MSIScore", -math.inf) >= 20},
        "wgd": {m for m, row in rows.items() if row.get("WGD") == 1},
        "cin_high": {m for m, row in rows.items() if row.get("CIN", -math.inf) >= thresholds["cin_high"]},
        "aneuploidy_high": {
            m for m, row in rows.items() if row.get("Aneuploidy", -math.inf) >= thresholds["aneuploidy_high"]
        },
        "hyperploid": {m for m, row in rows.items() if row.get("Ploidy", -math.inf) >= 3},
        "loh_high": {m for m, row in rows.items() if row.get("LoHFraction", -math.inf) >= thresholds["loh_high"]},
    }
    eligible = {
        "msi_high": {m for m, row in rows.items() if "MSIScore" in row},
        "wgd": {m for m, row in rows.items() if "WGD" in row},
        "cin_high": {
            m
            for m, row in rows.items()
            if "CIN" in row and (row["CIN"] <= thresholds["cin_low"] or row["CIN"] >= thresholds["cin_high"])
        },
        "aneuploidy_high": {
            m
            for m, row in rows.items()
            if "Aneuploidy" in row
            and (
                row["Aneuploidy"] <= thresholds["aneuploidy_low"]
                or row["Aneuploidy"] >= thresholds["aneuploidy_high"]
            )
        },
        "hyperploid": {
            m
            for m, row in rows.items()
            if "Ploidy" in row and (1.5 <= row["Ploidy"] <= 2.5 or row["Ploidy"] >= 3)
        },
        "loh_high": {
            m
            for m, row in rows.items()
            if "LoHFraction" in row
            and (
                row["LoHFraction"] <= thresholds["loh_low"]
                or row["LoHFraction"] >= thresholds["loh_high"]
            )
        },
    }
    extras = {
        analysis_id: {
            model_id: {"grouping_note": ", ".join(f"{k}={v:.3g}" for k, v in rows[model_id].items())}
            for model_id in positives
        }
        for analysis_id, positives in groups.items()
    }
    return groups, eligible, extras, thresholds


def expression_state_groups(
    expression_path: str,
    model_info: dict[str, dict[str, str]],
) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, dict[str, object]], dict[str, float]]:
    required_genes = sorted(set(EPITHELIAL_MARKERS + MESENCHYMAL_MARKERS + NEUROENDOCRINE_MARKERS))
    scores: dict[str, dict[str, float]] = {}
    with open(expression_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        indexes = {
            parse_gene(label)[0]: i
            for i, label in enumerate(header)
            if parse_gene(label)[0] in required_genes
        }
        missing = set(required_genes) - set(indexes)
        if missing:
            raise RuntimeError(f"Missing expression markers: {', '.join(sorted(missing))}")
        for row in reader:
            if not row:
                continue
            model_id = row[0]
            values = {gene: float(row[indexes[gene]]) for gene in required_genes if row[indexes[gene]] not in {"", "NA"}}
            epithelial = [values[g] for g in EPITHELIAL_MARKERS if g in values]
            mesenchymal = [values[g] for g in MESENCHYMAL_MARKERS if g in values]
            neuroendocrine = [values[g] for g in NEUROENDOCRINE_MARKERS if g in values]
            if epithelial and mesenchymal:
                scores.setdefault(model_id, {})["EMT"] = sum(mesenchymal) / len(mesenchymal) - sum(epithelial) / len(epithelial)
            if neuroendocrine:
                scores.setdefault(model_id, {})["NE"] = sum(neuroendocrine) / len(neuroendocrine)

    carcinoma_terms = ("carcinoma", "adenocarcinoma", "squamous", "epithelial", "mesothelioma")
    emt_models = {
        model_id
        for model_id, score in scores.items()
        if "EMT" in score
        and model_id in model_info
        and any(term in model_info[model_id]["disease"].lower() for term in carcinoma_terms)
    }
    neuroendocrine_lineages = {"Lung", "Prostate", "Bowel", "Pancreas", "Esophagus/Stomach"}
    ne_models = {
        model_id
        for model_id, score in scores.items()
        if "NE" in score
        and model_id in model_info
        and model_info[model_id]["lineage"] in neuroendocrine_lineages
    }
    thresholds = {
        "emt_low": quantile([scores[m]["EMT"] for m in emt_models], 0.25),
        "emt_high": quantile([scores[m]["EMT"] for m in emt_models], 0.75),
        "ne_low": quantile([scores[m]["NE"] for m in ne_models], 0.25),
        "ne_high": quantile([scores[m]["NE"] for m in ne_models], 0.75),
    }
    groups = {
        "emt_high": {m for m in emt_models if scores[m]["EMT"] >= thresholds["emt_high"]},
        "neuroendocrine_high": {m for m in ne_models if scores[m]["NE"] >= thresholds["ne_high"]},
    }
    eligible = {
        "emt_high": {
            m for m in emt_models if scores[m]["EMT"] <= thresholds["emt_low"] or scores[m]["EMT"] >= thresholds["emt_high"]
        },
        "neuroendocrine_high": {
            m for m in ne_models if scores[m]["NE"] <= thresholds["ne_low"] or scores[m]["NE"] >= thresholds["ne_high"]
        },
    }
    extras = {
        "emt_high": {m: {"grouping_note": f"EMT score={scores[m]['EMT']:.3f}"} for m in groups["emt_high"]},
        "neuroendocrine_high": {
            m: {"grouping_note": f"neuroendocrine score={scores[m]['NE']:.3f}"}
            for m in groups["neuroendocrine_high"]
        },
    }
    return groups, eligible, extras, thresholds


def fusion_groups(fusion_path: str) -> dict[str, set[str]]:
    groups = {
        "bcr_abl": set(),
        "alk": set(),
        "ret_fusion": set(),
        "ros1_fusion": set(),
        "ews_fusion": set(),
        "kmt2a_fusion": set(),
        "ntrk_fusion": set(),
        "pax_foxo1_fusion": set(),
        "fgfr_fusion": set(),
    }
    with open(fusion_path, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("IsDefaultEntryForModel") != "Yes":
                continue
            model_id = row["ModelID"]
            genes = {
                parse_gene(row.get("Gene1", ""))[0].upper(),
                parse_gene(row.get("Gene2", ""))[0].upper(),
            }
            if {"BCR", "ABL1"} <= genes:
                groups["bcr_abl"].add(model_id)
            if "ALK" in genes:
                groups["alk"].add(model_id)
            if "RET" in genes:
                groups["ret_fusion"].add(model_id)
            if "ROS1" in genes:
                groups["ros1_fusion"].add(model_id)
            if "EWSR1" in genes:
                groups["ews_fusion"].add(model_id)
            if "KMT2A" in genes:
                groups["kmt2a_fusion"].add(model_id)
            if genes & {"NTRK1", "NTRK2", "NTRK3"}:
                groups["ntrk_fusion"].add(model_id)
            if genes & {"PAX3", "PAX7"} and "FOXO1" in genes:
                groups["pax_foxo1_fusion"].add(model_id)
            if genes & {"FGFR1", "FGFR2", "FGFR3"}:
                groups["fgfr_fusion"].add(model_id)
    return groups


def row_differentials(
    matrix_path: str,
    analyses: dict[str, set[str]],
    eligible_models: dict[str, set[str]] | None = None,
) -> tuple[dict[str, list[dict[str, object]]], dict[str, dict[str, object]]]:
    eligible_models = eligible_models or {}
    sums = {key: None for key in analyses}
    sum_squares = {key: None for key in analyses}
    counts = {key: None for key in analyses}
    positive_rows = {key: [] for key in analyses}
    eligible_sums = {key: None for key in eligible_models}
    eligible_sum_squares = {key: None for key in eligible_models}
    eligible_counts = {key: None for key in eligible_models}
    eligible_rows = {key: 0 for key in eligible_models}

    with open(matrix_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        labels = header[1:]
        total_sums = [0.0] * len(labels)
        total_sum_squares = [0.0] * len(labels)
        total_counts = [0] * len(labels)
        total_rows = 0
        for key in analyses:
            sums[key] = [0.0] * len(labels)
            sum_squares[key] = [0.0] * len(labels)
            counts[key] = [0] * len(labels)
        for key in eligible_models:
            eligible_sums[key] = [0.0] * len(labels)
            eligible_sum_squares[key] = [0.0] * len(labels)
            eligible_counts[key] = [0] * len(labels)

        for row in reader:
            if not row:
                continue
            total_rows += 1
            model_id = row[0]
            parsed = []
            for value in row[1:]:
                if value == "" or value.upper() == "NA":
                    parsed.append(None)
                else:
                    try:
                        parsed.append(float(value))
                    except ValueError:
                        parsed.append(None)

            for i, value in enumerate(parsed):
                if value is None:
                    continue
                total_sums[i] += value
                total_sum_squares[i] += value * value
                total_counts[i] += 1

            for key, positives in analyses.items():
                if key in eligible_models:
                    if model_id not in eligible_models[key]:
                        continue
                    eligible_rows[key] += 1
                    for i, value in enumerate(parsed):
                        if value is None:
                            continue
                        eligible_sums[key][i] += value
                        eligible_sum_squares[key][i] += value * value
                        eligible_counts[key][i] += 1
                if model_id not in positives:
                    continue
                positive_rows[key].append(model_id)
                for i, value in enumerate(parsed):
                    if value is None:
                        continue
                    sums[key][i] += value
                    sum_squares[key][i] += value * value
                    counts[key][i] += 1

    datasets: dict[str, list[dict[str, object]]] = {}
    for key in analyses:
        rows = []
        analysis_row_count = eligible_rows[key] if key in eligible_models else total_rows
        analysis_sums = eligible_sums[key] if key in eligible_models else total_sums
        analysis_sum_squares = (
            eligible_sum_squares[key] if key in eligible_models else total_sum_squares
        )
        analysis_counts = eligible_counts[key] if key in eligible_models else total_counts
        negative_row_count = analysis_row_count - len(positive_rows[key])
        min_positive_n = min_values_for_group(len(positive_rows[key]))
        min_negative_n = min_values_for_group(negative_row_count)
        for i, label in enumerate(labels):
            pos_n = counts[key][i]
            neg_n = analysis_counts[i] - pos_n
            if pos_n < min_positive_n or neg_n < min_negative_n:
                continue
            positive_average = sums[key][i] / pos_n
            negative_sum = analysis_sums[i] - sums[key][i]
            negative_average = negative_sum / neg_n
            raw_difference = positive_average - negative_average
            effect_size = hedges_g_from_moments(
                sums[key][i],
                sum_squares[key][i],
                pos_n,
                negative_sum,
                analysis_sum_squares[i] - sum_squares[key][i],
                neg_n,
            )
            if effect_size is None:
                continue
            gene, entrez = parse_gene(label)
            rows.append(
                {
                    "gene": gene,
                    "entrez": entrez,
                    "label": label,
                    "score": round(effect_size, 6),
                    "raw_difference": round(raw_difference, 6),
                    "positive_average": round(positive_average, 6),
                    "negative_average": round(negative_average, 6),
                    "positive_n": pos_n,
                    "negative_n": neg_n,
                }
            )
        add_ranks(rows)
        datasets[key] = rows

    included = {
        key: {
            "positive": sorted(positive_rows[key]),
            "negative_n": (eligible_rows[key] if key in eligible_models else total_rows) - len(positive_rows[key]),
            "min_positive_n": min_values_for_group(len(positive_rows[key])),
            "min_negative_n": min_values_for_group(
                (eligible_rows[key] if key in eligible_models else total_rows) - len(positive_rows[key])
            ),
        }
        for key in analyses
    }
    return datasets, included


def column_differentials(
    matrix_path: str,
    analyses: dict[str, set[str]],
    eligible_models: dict[str, set[str]] | None = None,
) -> tuple[dict[str, list[dict[str, object]]], dict[str, dict[str, object]]]:
    eligible_models = eligible_models or {}
    datasets: dict[str, list[dict[str, object]]] = {key: [] for key in analyses}
    included: dict[str, dict[str, object]] = {}

    with open(matrix_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        indexes: dict[str, tuple[list[int], list[int], list[str]]] = {}
        minimums: dict[str, tuple[int, int]] = {}
        for key, positives in analyses.items():
            pos_indexes = []
            neg_indexes = []
            pos_columns = []
            for i, column in enumerate(header[1:], start=1):
                if key in eligible_models and column not in eligible_models[key]:
                    continue
                if column in positives:
                    pos_indexes.append(i)
                    pos_columns.append(column)
                else:
                    neg_indexes.append(i)
            indexes[key] = (pos_indexes, neg_indexes, pos_columns)
            minimums[key] = (min_values_for_group(len(pos_indexes)), min_values_for_group(len(neg_indexes)))
            included[key] = {
                "positive": sorted(pos_columns),
                "negative_n": len(neg_indexes),
                "min_positive_n": minimums[key][0],
                "min_negative_n": minimums[key][1],
            }

        for row in reader:
            gene, entrez = parse_gene(row[0])
            all_values = values_at(row, range(1, len(row)))
            all_sum = sum(all_values)
            all_sum_squares = sum(value * value for value in all_values)
            all_count = len(all_values)
            for key, (pos_indexes, neg_indexes, _) in indexes.items():
                positive_values = values_at(row, pos_indexes)
                positive_sum = sum(positive_values)
                positive_sum_squares = sum(value * value for value in positive_values)
                if key in eligible_models:
                    negative_values = values_at(row, neg_indexes)
                    negative_sum = sum(negative_values)
                    negative_sum_squares = sum(value * value for value in negative_values)
                    negative_count = len(negative_values)
                else:
                    negative_sum = all_sum - positive_sum
                    negative_sum_squares = all_sum_squares - positive_sum_squares
                    negative_count = all_count - len(positive_values)
                min_positive_n, min_negative_n = minimums[key]
                if len(positive_values) < min_positive_n or negative_count < min_negative_n:
                    continue
                positive_average = positive_sum / len(positive_values)
                negative_average = negative_sum / negative_count
                raw_difference = positive_average - negative_average
                effect_size = hedges_g_from_moments(
                    positive_sum,
                    positive_sum_squares,
                    len(positive_values),
                    negative_sum,
                    negative_sum_squares,
                    negative_count,
                )
                if effect_size is None:
                    continue
                datasets[key].append(
                    {
                        "gene": gene,
                        "entrez": entrez,
                        "label": row[0],
                        "score": round(effect_size, 6),
                        "raw_difference": round(raw_difference, 6),
                        "positive_average": round(positive_average, 6),
                        "negative_average": round(negative_average, 6),
                        "positive_n": len(positive_values),
                        "negative_n": negative_count,
                    }
                )

    for rows in datasets.values():
        add_ranks(rows)
    return datasets, included


def values_at(row: list[str], indexes) -> list[float]:
    values = []
    for i in indexes:
        if i >= len(row):
            continue
        value = row[i]
        if value == "" or value.upper() == "NA":
            continue
        try:
            values.append(float(value))
        except ValueError:
            continue
    return values


def hedges_g_from_moments(
    positive_sum: float,
    positive_sum_squares: float,
    positive_n: int,
    negative_sum: float,
    negative_sum_squares: float,
    negative_n: int,
) -> float | None:
    if positive_n < 2 or negative_n < 2:
        return None
    positive_average = positive_sum / positive_n
    negative_average = negative_sum / negative_n
    positive_variance = max(
        0.0,
        (positive_sum_squares - positive_sum * positive_sum / positive_n) / (positive_n - 1),
    )
    negative_variance = max(
        0.0,
        (negative_sum_squares - negative_sum * negative_sum / negative_n) / (negative_n - 1),
    )
    degrees_of_freedom = positive_n + negative_n - 2
    pooled_variance = (
        (positive_n - 1) * positive_variance + (negative_n - 1) * negative_variance
    ) / degrees_of_freedom
    if pooled_variance <= 1e-12:
        return 0.0 if abs(positive_average - negative_average) <= 1e-12 else None
    small_sample_correction = 1.0 - 3.0 / (4.0 * degrees_of_freedom - 1.0)
    return small_sample_correction * (positive_average - negative_average) / math.sqrt(pooled_variance)


def add_ranks(rows: list[dict[str, object]]) -> None:
    rows.sort(key=lambda x: (float(x["score"]), str(x["gene"])))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank


def min_values_for_group(group_size: int) -> int:
    return max(MIN_GROUP_VALUES, math.ceil(max(0, group_size) * MIN_GROUP_COVERAGE))


def compact_rows(rows: list[dict[str, object]]) -> list[list[object]]:
    return [
        [
            row["gene"],
            row["score"],
            row["positive_average"],
            row["negative_average"],
            row["positive_n"],
            row["negative_n"],
            row["rank"],
            row["raw_difference"],
        ]
        for row in rows
    ]


def analysis_models(model_info: dict[str, dict[str, str]], model_ids: set[str], extra: dict[str, object] | None = None) -> list[dict[str, object]]:
    rows = []
    for model_id in sorted(model_ids):
        if model_id not in model_info:
            continue
        item: dict[str, object] = dict(model_info[model_id])
        if extra and model_id in extra:
            item.update(extra[model_id])
        rows.append(item)
    rows.sort(key=lambda x: str(x["cell_line"]))
    return rows


def main() -> None:
    ensure_dirs()
    index = read_index()
    files = {
        "model": find_file(index, CRISPR_RELEASE, "Model.csv"),
        "crispr": find_file(index, CRISPR_RELEASE, "CRISPRGeneEffect.csv"),
        "rnai": find_file(index, RNAI_RELEASE, "D2_combined_gene_dep_scores.csv"),
        "mutations": find_file(index, CRISPR_RELEASE, "OmicsSomaticMutations.csv"),
        "fusions": find_file(index, CRISPR_RELEASE, "OmicsFusionFiltered.csv"),
    }

    local_paths = {}
    for key, meta in files.items():
        local_path = os.path.join(DATA_DIR, meta["filename"])
        download(meta["url"], local_path)
        local_paths[key] = local_path
    local_paths["copy_number_22q1"] = os.path.join(DATA_DIR, MTAP_COPY_NUMBER_FILENAME)
    download(MTAP_COPY_NUMBER_URL, local_paths["copy_number_22q1"])
    local_paths["global_signatures_24q2"] = os.path.join(DATA_DIR, GLOBAL_SIGNATURE_FILENAME)
    download(GLOBAL_SIGNATURE_URL, local_paths["global_signatures_24q2"])
    local_paths["expression_24q2"] = os.path.join(DATA_DIR, EXPRESSION_FILENAME)
    download(EXPRESSION_URL, local_paths["expression_24q2"])

    model_info = model_lookup(local_paths["model"])
    hpv_transformants = fetch_hpv_transformants()
    hpv_models = {
        model_id
        for model_id, info in model_info.items()
        if info["cellosaurus"] in hpv_transformants
    }
    hpv_extra = {
        model_id: {"grouping_note": ", ".join(hpv_transformants[info["cellosaurus"]])}
        for model_id, info in model_info.items()
        if info["cellosaurus"] in hpv_transformants
    }

    mutation_sets, top_driver_genes, mutations_by_gene, likely_biallelic_lof = mutation_groups(
        local_paths["mutations"]
    )
    fusion_sets = fusion_groups(local_paths["fusions"])
    copy_sets, copy_measured, copy_values = copy_number_groups(local_paths["copy_number_22q1"])
    signature_sets, signature_eligible, signature_extras, signature_thresholds = global_signature_groups(
        local_paths["global_signatures_24q2"]
    )
    expression_sets, expression_eligible, expression_extras, expression_thresholds = expression_state_groups(
        local_paths["expression_24q2"], model_info
    )

    all_models = set(model_info)
    analysis_defs: list[dict[str, object]] = []

    def add_analysis(
        analysis_id: str,
        label: str,
        positive_label: str,
        negative_label: str,
        source: str,
        category: str,
        positive_models: set[str],
        *,
        eligible_models: set[str] | None = None,
        prevalence_models: set[str] | None = None,
        prevalence_denominator: str = "all DepMap cancer models",
        extra: dict[str, object] | None = None,
        prevalence_rank: int | None = None,
    ) -> None:
        positives = set(positive_models) & all_models
        eligible = set(eligible_models) & all_models if eligible_models is not None else None
        if eligible is not None:
            positives &= eligible
        if len(positives) < MIN_GROUP_VALUES:
            print(
                f"Skipping {label}: {len(positives)} positive models is below "
                f"the minimum cohort size of {MIN_GROUP_VALUES}"
            )
            return
        denominator = set(prevalence_models) & all_models if prevalence_models is not None else (
            eligible if eligible is not None else all_models
        )
        definition: dict[str, object] = {
            "id": analysis_id,
            "label": label,
            "positive_label": positive_label,
            "negative_label": negative_label,
            "source": source,
            "category": category,
            "positive_model_ids": positives,
            "positive_models": analysis_models(model_info, positives, extra),
            "prevalence_total": len(denominator),
            "prevalence_denominator": prevalence_denominator,
        }
        if eligible is not None:
            definition["eligible_model_ids"] = eligible
        if prevalence_rank is not None:
            definition["prevalence_rank"] = prevalence_rank
        analysis_defs.append(definition)

    add_analysis(
        "hpv",
        "HPV+ vs HPV-",
        "HPV+",
        "HPV-",
        "Cellosaurus transformant-list",
        "Viral and expression-defined states",
        hpv_models,
        extra=hpv_extra,
    )
    add_analysis(
        "msi_high",
        "MSI-high vs microsatellite stable",
        "MSI-high",
        "Microsatellite stable",
        "DepMap Public 24Q2 OmicsSignatures.csv; MSIScore >= 20",
        "Viral and expression-defined states",
        signature_sets["msi_high"],
        eligible_models=signature_eligible["msi_high"],
        prevalence_models=signature_eligible["msi_high"],
        prevalence_denominator="models with MSIScore",
        extra=signature_extras["msi_high"],
    )
    add_analysis(
        "emt_high",
        "Mesenchymal-high vs epithelial-high",
        "Mesenchymal-high",
        "Epithelial-high",
        (
            "DepMap Public 24Q2 expression; top vs bottom quartile of mean mesenchymal markers "
            "minus mean epithelial markers within carcinoma models"
        ),
        "Viral and expression-defined states",
        expression_sets["emt_high"],
        eligible_models=expression_eligible["emt_high"],
        prevalence_models={
            model_id
            for model_id in expression_eligible["emt_high"] | expression_sets["emt_high"]
            if model_id in model_info
        },
        prevalence_denominator="high-confidence epithelial or mesenchymal carcinoma models",
        extra=expression_extras["emt_high"],
    )
    add_analysis(
        "neuroendocrine_high",
        "Neuroendocrine-high vs neuroendocrine-low",
        "Neuroendocrine-high",
        "Neuroendocrine-low",
        (
            "DepMap Public 24Q2 expression; top vs bottom quartile of ASCL1, NEUROD1, INSM1, "
            "CHGA, SYP and NCAM1 within lung, prostate and gastrointestinal models"
        ),
        "Viral and expression-defined states",
        expression_sets["neuroendocrine_high"],
        eligible_models=expression_eligible["neuroendocrine_high"],
        prevalence_denominator="high-confidence neuroendocrine-high or -low models",
        extra=expression_extras["neuroendocrine_high"],
    )

    global_state_definitions = [
        ("wgd", "Whole-genome duplication vs no WGD", "WGD", "No WGD", signature_eligible["wgd"]),
        ("cin_high", "Chromosomal instability high vs low", "CIN-high", "CIN-low", signature_eligible["cin_high"]),
        (
            "aneuploidy_high",
            "Aneuploidy high vs low",
            "Aneuploidy-high",
            "Aneuploidy-low",
            signature_eligible["aneuploidy_high"],
        ),
        (
            "hyperploid",
            "Hyperploid vs near-diploid",
            "Ploidy >= 3",
            "Near-diploid",
            signature_eligible["hyperploid"],
        ),
        (
            "loh_high",
            "Genomic LOH-high vs low (HRD proxy)",
            "LOH-high",
            "LOH-low",
            signature_eligible["loh_high"],
        ),
    ]
    for analysis_id, label, positive_label, negative_label, eligible in global_state_definitions:
        add_analysis(
            analysis_id,
            label,
            positive_label,
            negative_label,
            "DepMap Public 24Q2 OmicsSignatures.csv; extreme quartiles used for continuous signatures",
            "Global genomic states",
            signature_sets[analysis_id],
            eligible_models=eligible,
            prevalence_denominator="models in the defined genomic-state comparison",
            extra=signature_extras[analysis_id],
        )

    amplification_labels = {
        "MYC": "MYC",
        "CCNE1": "CCNE1",
        "CCND1": "CCND1",
        "ERBB2": "ERBB2 (HER2)",
        "EGFR": "EGFR",
        "MDM2": "MDM2",
        "MET": "MET",
        "MYCN": "MYCN",
        "CDK4": "CDK4",
        "FGFR1": "FGFR1",
        "FGFR2": "FGFR2",
        "FGFR3": "FGFR3",
        "SOX2": "SOX2",
        "TP63": "TP63",
    }
    for gene in AMPLIFICATION_GENES:
        add_analysis(
            f"{gene.lower()}_amp",
            f"{amplification_labels[gene]} amplified vs non-amplified",
            f"{gene} amplified",
            f"{gene} non-amplified",
            f"DepMap Public 22Q1 CCLE_gene_cn.csv; ploidy-relative linear copy number >= {AMPLIFICATION_THRESHOLD}",
            "Copy-number amplifications",
            copy_sets[f"{gene.lower()}_amp"],
            eligible_models=copy_measured[gene],
            prevalence_models=copy_measured[gene],
            prevalence_denominator=f"models with {gene} copy-number data",
        )

    loss_sets = {
        "pten_loss": copy_sets["pten_del"] | mutations_by_gene.get("PTEN", set()),
        "rb1_loss": copy_sets["rb1_del"] | mutations_by_gene.get("RB1", set()),
        "smad4_loss": copy_sets["smad4_del"] | mutations_by_gene.get("SMAD4", set()),
        "nf1_loss": copy_sets["nf1_del"] | mutations_by_gene.get("NF1", set()),
        "keap1_loss": copy_sets["keap1_del"] | mutations_by_gene.get("KEAP1", set()),
    }
    add_analysis(
        "mtap_del",
        "MTAP deleted vs intact",
        "MTAP deleted",
        "MTAP intact",
        "DepMap Public 22Q1 CCLE_gene_cn.csv; MTAP linear copy number < 0.4",
        "Tumor-suppressor and pathway loss",
        copy_sets["mtap_del"],
        eligible_models=copy_measured["MTAP"],
        prevalence_models=copy_measured["MTAP"],
        prevalence_denominator="models with MTAP copy-number data",
    )
    add_analysis(
        "cdkn2ab_del",
        "CDKN2A/B deep deletion vs intact",
        "CDKN2A/B deep deletion",
        "CDKN2A/B intact",
        f"DepMap Public 22Q1 CCLE_gene_cn.csv; either gene linear copy number < {DEEP_DELETION_THRESHOLD}",
        "Tumor-suppressor and pathway loss",
        copy_sets["cdkn2ab_del"],
        eligible_models=copy_measured["CDKN2A"] & copy_measured["CDKN2B"],
        prevalence_models=copy_measured["CDKN2A"] & copy_measured["CDKN2B"],
        prevalence_denominator="models with CDKN2A and CDKN2B copy-number data",
    )
    for analysis_id, gene in [
        ("pten_loss", "PTEN"),
        ("rb1_loss", "RB1"),
        ("smad4_loss", "SMAD4"),
        ("nf1_loss", "NF1"),
    ]:
        add_analysis(
            analysis_id,
            f"{gene} loss vs retained",
            f"{gene} loss",
            f"{gene} retained",
            (
                f"DepMap driver mutations plus Public 22Q1 deep deletion "
                f"(linear copy number < {DEEP_DELETION_THRESHOLD})"
            ),
            "Tumor-suppressor and pathway loss",
            loss_sets[analysis_id],
            eligible_models=copy_measured[gene],
            prevalence_models=copy_measured[gene],
            prevalence_denominator=f"models with {gene} copy-number data",
        )

    tp53_loss = mutations_by_gene.get("TP53", set())
    rb1_loss = loss_sets["rb1_loss"]
    kras_mutant = mutations_by_gene.get("KRAS", set())
    composite_definitions = [
        (
            "tp53_rb1_dual_loss",
            "TP53/RB1 dual loss vs other models",
            "TP53/RB1 dual loss",
            "Other models",
            tp53_loss & rb1_loss,
            None,
            "DepMap driver mutations plus RB1 deep deletion",
        ),
        (
            "kras_stk11",
            "KRAS + STK11 altered vs KRAS alone",
            "KRAS/STK11 co-altered",
            "KRAS mutant / STK11 wild-type",
            kras_mutant & mutations_by_gene.get("STK11", set()),
            kras_mutant,
            "DepMap explicit driver mutations; restricted to KRAS-mutant models",
        ),
        (
            "kras_keap1",
            "KRAS + KEAP1 altered vs KRAS alone",
            "KRAS/KEAP1 co-altered",
            "KRAS mutant / KEAP1 wild-type",
            kras_mutant & loss_sets["keap1_loss"],
            kras_mutant,
            "DepMap driver mutations plus KEAP1 deep deletion; restricted to KRAS-mutant models",
        ),
        (
            "kras_tp53",
            "KRAS + TP53 altered vs KRAS alone",
            "KRAS/TP53 co-altered",
            "KRAS mutant / TP53 wild-type",
            kras_mutant & tp53_loss,
            kras_mutant,
            "DepMap explicit driver mutations; restricted to KRAS-mutant models",
        ),
        (
            "mtap_cdkn2a_codeletion",
            "MTAP/CDKN2A co-deletion vs CDKN2A-only deletion",
            "MTAP/CDKN2A co-deletion",
            "CDKN2A-only deletion",
            copy_sets["mtap_cdkn2a_codeletion"],
            copy_sets["mtap_cdkn2a_codeletion"] | copy_sets["cdkn2a_only_deletion"],
            "DepMap Public 22Q1 copy number; both genes < 0.4 versus CDKN2A < 0.4 with MTAP retained",
        ),
        (
            "mdm2_amp_tp53_wt",
            "MDM2 amplification in TP53 wild-type models",
            "MDM2 amplified / TP53 WT",
            "MDM2 non-amplified / TP53 WT",
            copy_sets["mdm2_amp"] - tp53_loss,
            copy_measured["MDM2"] - tp53_loss,
            "DepMap Public 22Q1 MDM2 copy number; restricted to TP53 driver-wild-type models",
        ),
    ]
    for analysis_id, label, positive_label, negative_label, positives, eligible, source in composite_definitions:
        add_analysis(
            analysis_id,
            label,
            positive_label,
            negative_label,
            source,
            "Composite genotypes",
            positives,
            eligible_models=eligible,
            prevalence_denominator="models in the defined composite comparison" if eligible else "all DepMap cancer models",
        )

    brca_mutant = mutations_by_gene.get("BRCA1", set()) | mutations_by_gene.get("BRCA2", set())
    brca_biallelic_proxy = (
        likely_biallelic_lof.get("BRCA1", set())
        | likely_biallelic_lof.get("BRCA2", set())
        | copy_sets["brca1_del"]
        | copy_sets["brca2_del"]
    )
    add_analysis(
        "brca_biallelic_proxy",
        "BRCA1/2 likely biallelic loss vs other BRCA alterations",
        "Likely biallelic BRCA1/2 loss",
        "Other BRCA1/2 driver alteration",
        (
            "DepMap LikelyLoF mutation with allele fraction >= 0.5 or deep BRCA1/2 deletion; "
            "a cell-line biallelic-loss proxy"
        ),
        "Composite genotypes",
        brca_biallelic_proxy,
        eligible_models=brca_mutant | brca_biallelic_proxy,
        prevalence_denominator="models with a BRCA1/2 driver alteration or deep deletion",
    )
    add_analysis(
        "nrf2_pathway",
        "NFE2L2 activation or KEAP1 loss vs pathway-wild-type",
        "NRF2 pathway activated",
        "NRF2 pathway wild-type",
        "DepMap NFE2L2 driver alterations or KEAP1 driver/deep-loss alterations",
        "Composite genotypes",
        mutations_by_gene.get("NFE2L2", set()) | loss_sets["keap1_loss"],
    )
    add_analysis(
        "wnt_pathway",
        "APC loss or CTNNB1 activation vs pathway-wild-type",
        "WNT pathway activated",
        "WNT pathway wild-type",
        "DepMap explicit APC or CTNNB1 driver alterations",
        "Composite genotypes",
        mutations_by_gene.get("APC", set()) | mutations_by_gene.get("CTNNB1", set()),
    )

    mutation_sets["alk"] = mutation_sets["alk"] | fusion_sets["alk"]
    fusion_definitions = [
        ("bcr_abl", "BCR-ABL1 fusion", fusion_sets["bcr_abl"]),
        ("alk", "ALK alteration", mutation_sets["alk"]),
        ("ret_fusion", "RET fusion", fusion_sets["ret_fusion"]),
        ("ros1_fusion", "ROS1 fusion", fusion_sets["ros1_fusion"]),
        ("ews_fusion", "EWSR1 fusion", fusion_sets["ews_fusion"]),
        ("kmt2a_fusion", "KMT2A fusion", fusion_sets["kmt2a_fusion"]),
        ("ntrk_fusion", "NTRK1/2/3 fusion", fusion_sets["ntrk_fusion"]),
        ("pax_foxo1_fusion", "PAX3/7-FOXO1 fusion", fusion_sets["pax_foxo1_fusion"]),
        ("fgfr_fusion", "FGFR1/2/3 fusion", fusion_sets["fgfr_fusion"]),
    ]
    for analysis_id, name, positives in fusion_definitions:
        add_analysis(
            analysis_id,
            f"{name} vs alteration-negative",
            name,
            f"{name} negative",
            "DepMap OmicsFusionFiltered.csv" if analysis_id != "alk" else "DepMap ALK driver mutations plus ALK fusions",
            "Oncogenic fusions",
            positives,
        )

    for rank, gene in enumerate(top_driver_genes, start=1):
        key = gene.lower()
        add_analysis(
            key,
            f"{gene} mutant vs non-mutant",
            f"{gene} mutant",
            f"{gene} non-mutant",
            "DepMap OmicsSomaticMutations.csv explicit driver annotations",
            "Top driver mutations by prevalence",
            mutation_sets[key],
            prevalence_rank=rank,
        )
    for key, label, genes in ADDITIONAL_MUTATION_ANALYSES:
        add_analysis(
            key,
            label,
            f"{'/'.join(genes)} mutant",
            f"{'/'.join(genes)} non-mutant",
            "DepMap OmicsSomaticMutations.csv explicit driver annotations",
            "Additional mutation groups",
            mutation_sets[key],
        )

    model_to_ccle = {model_id: info["ccle_name"] for model_id, info in model_info.items() if info["ccle_name"]}
    crispr_groups = {d["id"]: set(d["positive_model_ids"]) for d in analysis_defs}
    rnai_groups = {
        d["id"]: {model_to_ccle[m] for m in d["positive_model_ids"] if m in model_to_ccle}
        for d in analysis_defs
    }
    crispr_eligible = {
        d["id"]: set(d["eligible_model_ids"])
        for d in analysis_defs
        if d.get("eligible_model_ids") is not None
    }
    rnai_eligible = {
        d["id"]: {model_to_ccle[m] for m in d["eligible_model_ids"] if m in model_to_ccle}
        for d in analysis_defs
        if d.get("eligible_model_ids") is not None
    }

    crispr, crispr_included = row_differentials(local_paths["crispr"], crispr_groups, crispr_eligible)
    rnai, rnai_included = column_differentials(local_paths["rnai"], rnai_groups, rnai_eligible)

    analyses = []
    analysis_counts = {}
    for definition in analysis_defs:
        analysis_id = definition["id"]
        included_models = {
            "crispr": crispr_included[analysis_id],
            "rnai": rnai_included[analysis_id],
        }
        analysis_data = {
            "datasets": {
                "crispr": compact_rows(crispr[analysis_id]),
                "rnai": compact_rows(rnai[analysis_id]),
            },
            "included_models": included_models,
        }
        with open(os.path.join(ANALYSIS_DATA_DIR, f"{analysis_id}.json"), "w") as f:
            json.dump(analysis_data, f, separators=(",", ":"))
        analyses.append(
            {
                "id": analysis_id,
                "label": definition["label"],
                "positive_label": definition["positive_label"],
                "negative_label": definition["negative_label"],
                "source": definition["source"],
                "category": definition["category"],
                "effect_metric": "hedges_g",
                "positive_models": definition["positive_models"],
                "data_url": f"/api/dependency-analysis/{analysis_id}",
                **(
                    {"prevalence_rank": definition["prevalence_rank"]}
                    if definition.get("prevalence_rank")
                    else {}
                ),
                "prevalence_total": definition["prevalence_total"],
                "prevalence_denominator": definition["prevalence_denominator"],
            }
        )
        analysis_counts[analysis_id] = included_models

    active_analysis_files = {f"{analysis['id']}.json" for analysis in analyses}
    for filename in os.listdir(ANALYSIS_DATA_DIR):
        if filename.endswith(".json") and filename not in active_analysis_files:
            os.remove(os.path.join(ANALYSIS_DATA_DIR, filename))

    payload = {
        "generated_at": datetime.now(tz=timezone.utc).replace(microsecond=0).isoformat(),
        "comparison": (
            "score = Hedges' g for positive-group versus negative-group dependency; "
            "negative values indicate greater essentiality in the positive group"
        ),
        "sources": {
            "crispr": {
                "release": CRISPR_RELEASE,
                "file": "CRISPRGeneEffect.csv",
                "release_date": files["crispr"]["release_date"],
            },
            "rnai": {
                "release": RNAI_RELEASE,
                "file": "D2_combined_gene_dep_scores.csv",
                "release_date": files["rnai"]["release_date"],
            },
            "models": {
                "release": CRISPR_RELEASE,
                "file": "Model.csv",
                "release_date": files["model"]["release_date"],
            },
        },
        "analyses": analyses,
    }

    with open(SUMMARY_JSON, "w") as f:
        json.dump(payload, f, separators=(",", ":"))
    print(f"Wrote {SUMMARY_JSON}")
    for analysis in analyses:
        included_models = analysis_counts[analysis["id"]]
        print(
            f"{analysis['label']}: CRISPR {len(included_models['crispr']['positive'])}+/"
            f"{included_models['crispr']['negative_n']}-; RNAi "
            f"{len(included_models['rnai']['positive'])}+/"
            f"{included_models['rnai']['negative_n']}-"
        )


if __name__ == "__main__":
    main()
