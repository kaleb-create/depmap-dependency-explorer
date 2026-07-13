import csv
import json
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

DOWNLOAD_INDEX_URL = "https://depmap.org/portal/api/download/files"
CELLOSAURUS_SEARCH_URL = "https://api.cellosaurus.org/search/cell-line"
CRISPR_RELEASE = "DepMap Public 26Q1"
RNAI_RELEASE = "DEMETER2 Data v6"
HPV_TRANSFORMANT_QUERY = "transformant:papillomavirus"

MUTATION_ANALYSES = [
    ("kras", "KRAS mutant vs non-mutant", ["KRAS"]),
    ("egfr", "EGFR mutant vs non-mutant", ["EGFR"]),
    ("pi3k", "PI3K mutant vs non-mutant", ["PIK3CA", "PIK3R1"]),
    ("braf", "BRAF mutant vs non-mutant", ["BRAF"]),
    ("nras", "NRAS mutant vs non-mutant", ["NRAS"]),
    ("hras", "HRAS mutant vs non-mutant", ["HRAS"]),
    ("tp53", "TP53 mutant vs non-mutant", ["TP53"]),
    ("pten", "PTEN mutant vs non-mutant", ["PTEN"]),
    ("apc", "APC mutant vs non-mutant", ["APC"]),
    ("ctnnb1", "CTNNB1 mutant vs non-mutant", ["CTNNB1"]),
]


def ensure_dirs() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(STATIC_DATA_DIR, exist_ok=True)


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
            row.get("VepImpact") in {"HIGH", "MODERATE"}
            or row.get("Hotspot") == "True"
            or row.get("HessDriver") == "True"
            or row.get("OncogeneHighImpact") == "True"
            or row.get("TumorSuppressorHighImpact") == "True"
            or row.get("LikelyLoF") == "True"
        )
    )


def mutation_groups(mutation_path: str) -> dict[str, set[str]]:
    by_gene = {gene: set() for _, _, genes in MUTATION_ANALYSES for gene in genes}
    by_gene["ALK"] = set()
    with open(mutation_path, newline="") as f:
        for row in csv.DictReader(f):
            gene = row.get("HugoSymbol")
            if gene not in by_gene or not is_driver_mutation(row):
                continue
            by_gene[gene].add(row["ModelID"])
    groups = {
        key: set().union(*(by_gene[gene] for gene in genes))
        for key, _, genes in MUTATION_ANALYSES
    }
    groups["alk"] = by_gene["ALK"]
    return groups


def fusion_groups(fusion_path: str) -> dict[str, set[str]]:
    groups = {
        "bcr_abl": set(),
        "alk": set(),
        "ret_fusion": set(),
        "ros1_fusion": set(),
    }
    with open(fusion_path, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("IsDefaultEntryForModel") != "Yes":
                continue
            fusion = " ".join([row.get("CanonicalFusionName", ""), row.get("Gene1", ""), row.get("Gene2", "")]).upper()
            if "BCR" in fusion and "ABL1" in fusion:
                groups["bcr_abl"].add(row["ModelID"])
            if "ALK" in fusion:
                groups["alk"].add(row["ModelID"])
            if "RET" in fusion:
                groups["ret_fusion"].add(row["ModelID"])
            if "ROS1" in fusion:
                groups["ros1_fusion"].add(row["ModelID"])
    return groups


def row_differentials(matrix_path: str, analyses: dict[str, set[str]]) -> tuple[dict[str, list[dict[str, object]]], dict[str, dict[str, object]]]:
    sums = {key: None for key in analyses}
    counts = {key: None for key in analyses}
    neg_sums = {key: None for key in analyses}
    neg_counts = {key: None for key in analyses}
    positive_rows = {key: [] for key in analyses}
    negative_counts_by_analysis = {key: 0 for key in analyses}

    with open(matrix_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        labels = header[1:]
        for key in analyses:
            sums[key] = [0.0] * len(labels)
            counts[key] = [0] * len(labels)
            neg_sums[key] = [0.0] * len(labels)
            neg_counts[key] = [0] * len(labels)

        for row in reader:
            if not row:
                continue
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

            for key, positives in analyses.items():
                is_positive = model_id in positives
                if is_positive:
                    positive_rows[key].append(model_id)
                else:
                    negative_counts_by_analysis[key] += 1
                target_sums = sums[key] if is_positive else neg_sums[key]
                target_counts = counts[key] if is_positive else neg_counts[key]
                for i, value in enumerate(parsed):
                    if value is None:
                        continue
                    target_sums[i] += value
                    target_counts[i] += 1

    datasets: dict[str, list[dict[str, object]]] = {}
    for key in analyses:
        rows = []
        for i, label in enumerate(labels):
            pos_n = counts[key][i]
            neg_n = neg_counts[key][i]
            if pos_n == 0 or neg_n == 0:
                continue
            positive_average = sums[key][i] / pos_n
            negative_average = neg_sums[key][i] / neg_n
            gene, entrez = parse_gene(label)
            rows.append(
                {
                    "gene": gene,
                    "entrez": entrez,
                    "label": label,
                    "score": round(positive_average - negative_average, 6),
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
            "negative_n": negative_counts_by_analysis[key],
        }
        for key in analyses
    }
    return datasets, included


def column_differentials(matrix_path: str, analyses: dict[str, set[str]]) -> tuple[dict[str, list[dict[str, object]]], dict[str, dict[str, object]]]:
    datasets: dict[str, list[dict[str, object]]] = {key: [] for key in analyses}
    included: dict[str, dict[str, object]] = {}

    with open(matrix_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        indexes: dict[str, tuple[list[int], list[int], list[str]]] = {}
        for key, positives in analyses.items():
            pos_indexes = []
            neg_indexes = []
            pos_columns = []
            for i, column in enumerate(header[1:], start=1):
                if column in positives:
                    pos_indexes.append(i)
                    pos_columns.append(column)
                else:
                    neg_indexes.append(i)
            indexes[key] = (pos_indexes, neg_indexes, pos_columns)
            included[key] = {"positive": sorted(pos_columns), "negative_n": len(neg_indexes)}

        for row in reader:
            gene, entrez = parse_gene(row[0])
            for key, (pos_indexes, neg_indexes, _) in indexes.items():
                positive_values = values_at(row, pos_indexes)
                negative_values = values_at(row, neg_indexes)
                if not positive_values or not negative_values:
                    continue
                positive_average = sum(positive_values) / len(positive_values)
                negative_average = sum(negative_values) / len(negative_values)
                datasets[key].append(
                    {
                        "gene": gene,
                        "entrez": entrez,
                        "label": row[0],
                        "score": round(positive_average - negative_average, 6),
                        "positive_average": round(positive_average, 6),
                        "negative_average": round(negative_average, 6),
                        "positive_n": len(positive_values),
                        "negative_n": len(negative_values),
                    }
                )

    for rows in datasets.values():
        add_ranks(rows)
    return datasets, included


def values_at(row: list[str], indexes: list[int]) -> list[float]:
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


def add_ranks(rows: list[dict[str, object]]) -> None:
    rows.sort(key=lambda x: (float(x["score"]), str(x["gene"])))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank


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

    mutation_sets = mutation_groups(local_paths["mutations"])
    fusion_sets = fusion_groups(local_paths["fusions"])
    mutation_sets["bcr_abl"] = fusion_sets["bcr_abl"]
    mutation_sets["alk"] = mutation_sets["alk"] | fusion_sets["alk"]
    mutation_sets["ret_fusion"] = fusion_sets["ret_fusion"]
    mutation_sets["ros1_fusion"] = fusion_sets["ros1_fusion"]

    analysis_defs = [
        {
            "id": "hpv",
            "label": "HPV+ vs HPV-",
            "positive_label": "HPV+",
            "negative_label": "HPV-",
            "source": "Cellosaurus transformant-list",
            "positive_model_ids": hpv_models,
            "positive_models": analysis_models(model_info, hpv_models, hpv_extra),
        },
        {
            "id": "bcr_abl",
            "label": "BCR-ABL fusion vs fusion-negative",
            "positive_label": "BCR-ABL fusion",
            "negative_label": "BCR-ABL fusion-negative",
            "source": "DepMap OmicsFusionFiltered.csv",
            "positive_model_ids": mutation_sets["bcr_abl"],
            "positive_models": analysis_models(model_info, mutation_sets["bcr_abl"]),
        },
        {
            "id": "alk",
            "label": "ALK altered vs unaltered",
            "positive_label": "ALK altered",
            "negative_label": "ALK unaltered",
            "source": "ALK mutations plus ALK fusions",
            "positive_model_ids": mutation_sets["alk"],
            "positive_models": analysis_models(model_info, mutation_sets["alk"]),
        },
    ]
    for key, label, genes in MUTATION_ANALYSES:
        analysis_defs.append(
            {
                "id": key,
                "label": label,
                "positive_label": f"{'/'.join(genes)} mutant",
                "negative_label": f"{'/'.join(genes)} non-mutant",
                "source": "DepMap OmicsSomaticMutations.csv protein-altering/driver variants",
                "positive_model_ids": mutation_sets[key],
                "positive_models": analysis_models(model_info, mutation_sets[key]),
            }
        )
    for key, label in [("ret_fusion", "RET fusion vs fusion-negative"), ("ros1_fusion", "ROS1 fusion vs fusion-negative")]:
        analysis_defs.append(
            {
                "id": key,
                "label": label,
                "positive_label": label.split(" vs ")[0],
                "negative_label": label.split(" vs ")[1],
                "source": "DepMap OmicsFusionFiltered.csv",
                "positive_model_ids": mutation_sets[key],
                "positive_models": analysis_models(model_info, mutation_sets[key]),
            }
        )

    model_to_ccle = {model_id: info["ccle_name"] for model_id, info in model_info.items() if info["ccle_name"]}
    crispr_groups = {d["id"]: set(d["positive_model_ids"]) for d in analysis_defs}
    rnai_groups = {
        d["id"]: {model_to_ccle[m] for m in d["positive_model_ids"] if m in model_to_ccle}
        for d in analysis_defs
    }

    crispr, crispr_included = row_differentials(local_paths["crispr"], crispr_groups)
    rnai, rnai_included = column_differentials(local_paths["rnai"], rnai_groups)

    analyses = []
    for definition in analysis_defs:
        analysis_id = definition["id"]
        analyses.append(
            {
                "id": analysis_id,
                "label": definition["label"],
                "positive_label": definition["positive_label"],
                "negative_label": definition["negative_label"],
                "source": definition["source"],
                "positive_models": definition["positive_models"],
                "datasets": {
                    "crispr": compact_rows(crispr[analysis_id]),
                    "rnai": compact_rows(rnai[analysis_id]),
                },
                "included_models": {
                    "crispr": crispr_included[analysis_id],
                    "rnai": rnai_included[analysis_id],
                },
            }
        )

    payload = {
        "generated_at": datetime.now(tz=timezone.utc).replace(microsecond=0).isoformat(),
        "comparison": "score = positive-group average dependency - negative-group average dependency",
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
        print(
            f"{analysis['label']}: CRISPR {len(analysis['included_models']['crispr']['positive'])}+/"
            f"{analysis['included_models']['crispr']['negative_n']}-; RNAi "
            f"{len(analysis['included_models']['rnai']['positive'])}+/"
            f"{analysis['included_models']['rnai']['negative_n']}-"
        )


if __name__ == "__main__":
    main()
