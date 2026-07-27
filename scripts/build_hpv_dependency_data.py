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

DOWNLOAD_INDEX_URL = "https://depmap.org/portal/api/download/files"
CELLOSAURUS_SEARCH_URL = "https://api.cellosaurus.org/search/cell-line"
CRISPR_RELEASE = "DepMap Public 26Q1"
RNAI_RELEASE = "DEMETER2 Data v6"
HPV_TRANSFORMANT_QUERY = "transformant:papillomavirus"
MTAP_COPY_NUMBER_URL = "https://ndownloader.figshare.com/files/34008428"
MTAP_COPY_NUMBER_FILENAME = "CCLE_gene_cn_22Q1.csv"
MTAP_DELETION_THRESHOLD = 0.4
TOP_DRIVER_COUNT = 20
MIN_GROUP_VALUES = 3
MIN_GROUP_COVERAGE = 0.5

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


def mutation_groups(mutation_path: str) -> tuple[dict[str, set[str]], list[str]]:
    by_gene: dict[str, set[str]] = {}
    with open(mutation_path, newline="") as f:
        for row in csv.DictReader(f):
            gene = row.get("HugoSymbol")
            if not gene or not is_driver_mutation(row):
                continue
            by_gene.setdefault(gene, set()).add(row["ModelID"])
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
    return groups, top_driver_genes


def copy_number_deletion_group(copy_number_path: str, gene: str, threshold: float) -> tuple[set[str], set[str]]:
    deleted: set[str] = set()
    measured: set[str] = set()
    with open(copy_number_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        gene_index = next(
            (i for i, label in enumerate(header) if parse_gene(label)[0] == gene),
            None,
        )
        if gene_index is None:
            raise RuntimeError(f"Missing {gene} in {os.path.basename(copy_number_path)}")
        for row in reader:
            if gene_index >= len(row) or row[gene_index] in {"", "NA"}:
                continue
            model_id = row[0]
            measured.add(model_id)
            if float(row[gene_index]) < threshold:
                deleted.add(model_id)
    return deleted, measured


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


def row_differentials(
    matrix_path: str,
    analyses: dict[str, set[str]],
    eligible_models: dict[str, set[str]] | None = None,
) -> tuple[dict[str, list[dict[str, object]]], dict[str, dict[str, object]]]:
    eligible_models = eligible_models or {}
    sums = {key: None for key in analyses}
    counts = {key: None for key in analyses}
    positive_rows = {key: [] for key in analyses}
    eligible_sums = {key: None for key in eligible_models}
    eligible_counts = {key: None for key in eligible_models}
    eligible_rows = {key: 0 for key in eligible_models}

    with open(matrix_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        labels = header[1:]
        total_sums = [0.0] * len(labels)
        total_counts = [0] * len(labels)
        total_rows = 0
        for key in analyses:
            sums[key] = [0.0] * len(labels)
            counts[key] = [0] * len(labels)
        for key in eligible_models:
            eligible_sums[key] = [0.0] * len(labels)
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
                        eligible_counts[key][i] += 1
                if model_id not in positives:
                    continue
                positive_rows[key].append(model_id)
                for i, value in enumerate(parsed):
                    if value is None:
                        continue
                    sums[key][i] += value
                    counts[key][i] += 1

    datasets: dict[str, list[dict[str, object]]] = {}
    for key in analyses:
        rows = []
        analysis_row_count = eligible_rows[key] if key in eligible_models else total_rows
        analysis_sums = eligible_sums[key] if key in eligible_models else total_sums
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
            negative_average = (analysis_sums[i] - sums[key][i]) / neg_n
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
            all_count = len(all_values)
            for key, (pos_indexes, neg_indexes, _) in indexes.items():
                positive_values = values_at(row, pos_indexes)
                if key in eligible_models:
                    negative_values = values_at(row, neg_indexes)
                    negative_sum = sum(negative_values)
                    negative_count = len(negative_values)
                else:
                    positive_sum = sum(positive_values)
                    negative_sum = all_sum - positive_sum
                    negative_count = all_count - len(positive_values)
                min_positive_n, min_negative_n = minimums[key]
                if len(positive_values) < min_positive_n or negative_count < min_negative_n:
                    continue
                positive_average = sum(positive_values) / len(positive_values)
                negative_average = negative_sum / negative_count
                datasets[key].append(
                    {
                        "gene": gene,
                        "entrez": entrez,
                        "label": row[0],
                        "score": round(positive_average - negative_average, 6),
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


def add_ranks(rows: list[dict[str, object]]) -> None:
    rows.sort(key=lambda x: (float(x["score"]), str(x["gene"])))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank


def min_values_for_group(group_size: int) -> int:
    if group_size <= 0:
        return 1
    return min(group_size, max(MIN_GROUP_VALUES, math.ceil(group_size * MIN_GROUP_COVERAGE)))


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
    local_paths["copy_number_22q1"] = os.path.join(DATA_DIR, MTAP_COPY_NUMBER_FILENAME)
    download(MTAP_COPY_NUMBER_URL, local_paths["copy_number_22q1"])

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

    mutation_sets, top_driver_genes = mutation_groups(local_paths["mutations"])
    fusion_sets = fusion_groups(local_paths["fusions"])
    mutation_sets["bcr_abl"] = fusion_sets["bcr_abl"]
    mutation_sets["alk"] = mutation_sets["alk"] | fusion_sets["alk"]
    mutation_sets["ret_fusion"] = fusion_sets["ret_fusion"]
    mutation_sets["ros1_fusion"] = fusion_sets["ros1_fusion"]
    mtap_deleted, mtap_measured = copy_number_deletion_group(
        local_paths["copy_number_22q1"],
        "MTAP",
        MTAP_DELETION_THRESHOLD,
    )
    mutation_sets["mtap_del"] = mtap_deleted

    analysis_defs = [
        {
            "id": "hpv",
            "label": "HPV+ vs HPV-",
            "positive_label": "HPV+",
            "negative_label": "HPV-",
            "source": "Cellosaurus transformant-list",
            "category": "Featured analyses",
            "positive_model_ids": hpv_models,
            "positive_models": analysis_models(model_info, hpv_models, hpv_extra),
        },
        {
            "id": "bcr_abl",
            "label": "BCR-ABL fusion vs fusion-negative",
            "positive_label": "BCR-ABL fusion",
            "negative_label": "BCR-ABL fusion-negative",
            "source": "DepMap OmicsFusionFiltered.csv",
            "category": "Featured analyses",
            "positive_model_ids": mutation_sets["bcr_abl"],
            "positive_models": analysis_models(model_info, mutation_sets["bcr_abl"]),
        },
        {
            "id": "alk",
            "label": "ALK altered vs unaltered",
            "positive_label": "ALK altered",
            "negative_label": "ALK unaltered",
            "source": "ALK mutations plus ALK fusions",
            "category": "Featured analyses",
            "positive_model_ids": mutation_sets["alk"],
            "positive_models": analysis_models(model_info, mutation_sets["alk"]),
        },
        {
            "id": "mtap_del",
            "label": "MTAP deleted vs intact",
            "positive_label": "MTAP deleted",
            "negative_label": "MTAP intact",
            "source": "DepMap Public 22Q1 CCLE_gene_cn.csv; MTAP linear copy number < 0.4",
            "category": "Featured analyses",
            "prevalence_total": len(mtap_measured & set(model_info)),
            "prevalence_denominator": "DepMap Public 22Q1 models with MTAP copy-number data",
            "positive_model_ids": mutation_sets["mtap_del"],
            "positive_models": analysis_models(model_info, mutation_sets["mtap_del"]),
        },
    ]
    for rank, gene in enumerate(top_driver_genes, start=1):
        key = gene.lower()
        analysis_defs.append(
            {
                "id": key,
                "label": f"{gene} mutant vs non-mutant",
                "positive_label": f"{gene} mutant",
                "negative_label": f"{gene} non-mutant",
                "source": "DepMap OmicsSomaticMutations.csv explicit driver annotations",
                "category": "Top 20 mutations by prevalence",
                "prevalence_rank": rank,
                "positive_model_ids": mutation_sets[key],
                "positive_models": analysis_models(model_info, mutation_sets[key]),
            }
        )
    for key, label, genes in ADDITIONAL_MUTATION_ANALYSES:
        analysis_defs.append(
            {
                "id": key,
                "label": label,
                "positive_label": f"{'/'.join(genes)} mutant",
                "negative_label": f"{'/'.join(genes)} non-mutant",
                "source": "DepMap OmicsSomaticMutations.csv explicit driver annotations",
                "category": "Additional key stratifiers",
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
                "category": "Additional key stratifiers",
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
    crispr_eligible = {"mtap_del": mtap_measured}
    rnai_eligible = {
        "mtap_del": {model_to_ccle[m] for m in mtap_measured if m in model_to_ccle}
    }

    crispr, crispr_included = row_differentials(local_paths["crispr"], crispr_groups, crispr_eligible)
    rnai, rnai_included = column_differentials(local_paths["rnai"], rnai_groups, rnai_eligible)

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
                "category": definition["category"],
                "positive_models": definition["positive_models"],
                **(
                    {"prevalence_rank": definition["prevalence_rank"]}
                    if definition.get("prevalence_rank")
                    else {}
                ),
                **(
                    {
                        "prevalence_total": definition["prevalence_total"],
                        "prevalence_denominator": definition["prevalence_denominator"],
                    }
                    if definition.get("prevalence_total")
                    else {}
                ),
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
