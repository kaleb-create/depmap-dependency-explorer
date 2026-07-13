import csv
import json
import os
import urllib.request
from typing import Any

from scripts.build_hpv_dependency_data import (
    compact_rows,
    parse_gene,
)


BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.environ.get("DEPMAP_DATA_DIR", os.path.join(BASE_DIR, "data"))
CRISPR_PATH = os.path.join(DATA_DIR, "CRISPRGeneEffect.csv")
RNAI_PATH = os.path.join(DATA_DIR, "D2_combined_gene_dep_scores.csv")
MODEL_PATH = os.path.join(DATA_DIR, "Model.csv")

OPENAI_API_URL = "https://api.openai.com/v1/responses"
OPENAI_MODEL = os.environ.get("OPENAI_MODEL", "gpt-5.5")


def load_model_rows() -> list[dict[str, str]]:
    if not os.path.exists(MODEL_PATH):
        raise FileNotFoundError(f"Missing DepMap Model.csv at {MODEL_PATH}. Set DEPMAP_DATA_DIR to the raw DepMap CSV directory.")

    with open(MODEL_PATH, newline="") as f:
        return list(csv.DictReader(f))


def model_matches(model: dict[str, str], spec: dict[str, Any]) -> bool:
    field_checks = [
        ("OncotreeLineage", spec.get("lineages") or []),
        ("OncotreePrimaryDisease", spec.get("diseases") or []),
        ("OncotreeSubtype", spec.get("subtypes") or []),
    ]
    has_structured_filter = any(values for _, values in field_checks)

    if has_structured_filter:
        matched = False
        for field, values in field_checks:
            actual = (model.get(field) or "").lower()
            if any(str(value).lower() == actual for value in values):
                matched = True
                break
        if not matched:
            return False

    haystack = " ".join(
        [
            model.get("CellLineName") or "",
            model.get("StrippedCellLineName") or "",
            model.get("OncotreeLineage") or "",
            model.get("OncotreePrimaryDisease") or "",
            model.get("OncotreeSubtype") or "",
            model.get("ModelSubtypeFeatures") or "",
        ]
    ).lower()

    include_terms = [str(term).lower() for term in spec.get("include_terms") or [] if str(term).strip()]
    exclude_terms = [str(term).lower() for term in spec.get("exclude_terms") or [] if str(term).strip()]

    if include_terms and not any(term in haystack for term in include_terms):
        return False
    if exclude_terms and any(term in haystack for term in exclude_terms):
        return False
    return True


def add_ranks(rows: list[dict[str, object]]) -> None:
    rows.sort(key=lambda x: (float(x["score"]), str(x["gene"])))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank


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


def row_pair_differential(matrix_path: str, positive_ids: set[str], negative_ids: set[str]) -> tuple[list[dict[str, object]], dict[str, Any]]:
    with open(matrix_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        labels = header[1:]
        positive_sums = [0.0] * len(labels)
        positive_counts = [0] * len(labels)
        negative_sums = [0.0] * len(labels)
        negative_counts = [0] * len(labels)
        positive_rows = []
        negative_rows = []

        for row in reader:
            if not row:
                continue
            if row[0] in positive_ids:
                sums = positive_sums
                counts = positive_counts
                positive_rows.append(row[0])
            elif row[0] in negative_ids:
                sums = negative_sums
                counts = negative_counts
                negative_rows.append(row[0])
            else:
                continue

            for i, value in enumerate(row[1:]):
                if value == "" or value.upper() == "NA":
                    continue
                try:
                    sums[i] += float(value)
                    counts[i] += 1
                except ValueError:
                    continue

    rows = []
    for i, label in enumerate(labels):
        if not positive_counts[i] or not negative_counts[i]:
            continue
        positive_average = positive_sums[i] / positive_counts[i]
        negative_average = negative_sums[i] / negative_counts[i]
        gene, _ = parse_gene(label)
        rows.append(
            {
                "gene": gene,
                "score": round(positive_average - negative_average, 6),
                "positive_average": round(positive_average, 6),
                "negative_average": round(negative_average, 6),
                "positive_n": positive_counts[i],
                "negative_n": negative_counts[i],
            }
        )
    add_ranks(rows)
    return rows, {"positive": sorted(positive_rows), "negative_n": len(negative_rows)}


def column_pair_differential(
    matrix_path: str,
    positive_columns: set[str],
    negative_columns: set[str],
) -> tuple[list[dict[str, object]], dict[str, Any]]:
    with open(matrix_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        positive_indexes = []
        negative_indexes = []
        positive_used = []
        negative_used = []
        for i, column in enumerate(header[1:], start=1):
            if column in positive_columns:
                positive_indexes.append(i)
                positive_used.append(column)
            elif column in negative_columns:
                negative_indexes.append(i)
                negative_used.append(column)

        rows = []
        for row in reader:
            positive_values = values_at(row, positive_indexes)
            negative_values = values_at(row, negative_indexes)
            if not positive_values or not negative_values:
                continue
            positive_average = sum(positive_values) / len(positive_values)
            negative_average = sum(negative_values) / len(negative_values)
            gene, _ = parse_gene(row[0])
            rows.append(
                {
                    "gene": gene,
                    "score": round(positive_average - negative_average, 6),
                    "positive_average": round(positive_average, 6),
                    "negative_average": round(negative_average, 6),
                    "positive_n": len(positive_values),
                    "negative_n": len(negative_values),
                }
            )
    add_ranks(rows)
    return rows, {"positive": sorted(positive_used), "negative_n": len(negative_used)}


def choose_models(model_rows: list[dict[str, str]], spec: dict[str, Any]) -> set[str]:
    return {row["ModelID"] for row in model_rows if model_matches(row, spec)}


def available_taxonomy(model_rows: list[dict[str, str]]) -> dict[str, list[str]]:
    taxonomy: dict[str, list[str]] = {}
    for key, column in {
        "lineages": "OncotreeLineage",
        "diseases": "OncotreePrimaryDisease",
        "subtypes": "OncotreeSubtype",
    }.items():
        values = sorted({row[column] for row in model_rows if row.get(column)})
        taxonomy[key] = values[:500]
    return taxonomy


def fallback_spec(prompt: str) -> dict[str, Any]:
    lowered = prompt.lower()
    left, _, right = lowered.partition(" vs ")
    if not right:
        left, _, right = lowered.partition(" versus ")
    if not right:
        left, _, right = lowered.partition(" compared to ")
    right = right or "all other"
    left = left.replace("cell lines", "").replace("cells", "").strip()
    right = right.replace("cell lines", "").replace("cells", "").strip()
    return {
        "label": f"{left.title()} vs {right.title()}",
        "positive_label": left.title() or "Positive",
        "negative_label": right.title() or "Negative",
        "positive": {"include_terms": [left], "exclude_terms": []},
        "negative": {"include_terms": [] if "all other" in right else [right], "exclude_terms": [left]},
        "quality": {
            "summary": "Created with a local keyword fallback because OPENAI_API_KEY is not configured.",
            "strengths": ["Uses DepMap model metadata directly."],
            "weaknesses": ["Keyword matching may miss aliases or include broader disease labels than intended."],
            "recommended_checks": ["Set OPENAI_API_KEY and run the quality review before using the cohort analytically."],
        },
    }


def request_openai_spec(prompt: str, taxonomy: dict[str, list[str]]) -> dict[str, Any]:
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        return fallback_spec(prompt)

    schema = {
        "type": "object",
        "additionalProperties": False,
        "required": ["label", "positive_label", "negative_label", "positive", "negative", "quality"],
        "properties": {
            "label": {"type": "string"},
            "positive_label": {"type": "string"},
            "negative_label": {"type": "string"},
            "positive": {"$ref": "#/$defs/group"},
            "negative": {"$ref": "#/$defs/group"},
            "quality": {
                "type": "object",
                "additionalProperties": False,
                "required": ["summary", "strengths", "weaknesses", "recommended_checks"],
                "properties": {
                    "summary": {"type": "string"},
                    "strengths": {"type": "array", "items": {"type": "string"}},
                    "weaknesses": {"type": "array", "items": {"type": "string"}},
                    "recommended_checks": {"type": "array", "items": {"type": "string"}},
                },
            },
        },
        "$defs": {
            "group": {
                "type": "object",
                "additionalProperties": False,
                "required": ["lineages", "diseases", "subtypes", "include_terms", "exclude_terms"],
                "properties": {
                    "lineages": {"type": "array", "items": {"type": "string"}},
                    "diseases": {"type": "array", "items": {"type": "string"}},
                    "subtypes": {"type": "array", "items": {"type": "string"}},
                    "include_terms": {"type": "array", "items": {"type": "string"}},
                    "exclude_terms": {"type": "array", "items": {"type": "string"}},
                },
            }
        },
    }
    body = {
        "model": OPENAI_MODEL,
        "input": [
            {
                "role": "system",
                "content": (
                    "Convert a plain-English cancer cell-line contrast into strict filters over DepMap "
                    "Model.csv metadata. Prefer exact Oncotree lineages, diseases, and subtypes from the "
                    "provided taxonomy. Use include_terms only for aliases not represented exactly. "
                    "Negative should represent the requested comparator, or all other models by excluding "
                    "the positive group if the prompt implies one-vs-rest."
                ),
            },
            {
                "role": "user",
                "content": json.dumps({"prompt": prompt, "available_taxonomy": taxonomy}),
            },
        ],
        "text": {"format": {"type": "json_schema", "name": "depmap_stratifier", "schema": schema, "strict": True}},
    }
    req = urllib.request.Request(
        OPENAI_API_URL,
        data=json.dumps(body).encode("utf-8"),
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=45) as response:
        payload = json.loads(response.read().decode("utf-8"))

    for item in payload.get("output", []):
        for content in item.get("content", []):
            if content.get("type") in {"output_text", "text"} and content.get("text"):
                return json.loads(content["text"])
    raise RuntimeError("OpenAI returned no structured text.")


def model_payload(model_rows: list[dict[str, str]], model_ids: set[str], quality_by_id: dict[str, str] | None = None) -> list[dict[str, object]]:
    quality_by_id = quality_by_id or {}
    rows = []
    for row in model_rows:
        if row["ModelID"] not in model_ids:
            continue
        item = {
            "model_id": row["ModelID"],
            "cell_line": row["CellLineName"] or row["StrippedCellLineName"],
            "stripped_name": row["StrippedCellLineName"],
            "ccle_name": row["CCLEName"],
            "cellosaurus": row["RRID"],
            "lineage": row["OncotreeLineage"],
            "disease": row["OncotreePrimaryDisease"],
        }
        if row["ModelID"] in quality_by_id:
            item["grouping_note"] = quality_by_id[row["ModelID"]]
        rows.append(item)
    return sorted(rows, key=lambda x: str(x["cell_line"]))


def compute_stratifier(prompt: str) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    model_rows = load_model_rows()
    spec = request_openai_spec(prompt, available_taxonomy(model_rows))
    positive_ids = choose_models(model_rows, spec["positive"])
    negative_ids = choose_models(model_rows, spec["negative"])
    if not negative_ids:
        negative_ids = {row["ModelID"] for row in model_rows} - positive_ids

    if len(positive_ids) < 2 or len(negative_ids) < 2:
        raise ValueError(
            f"Could only match {len(positive_ids)} positive models and {len(negative_ids)} negative models."
        )

    model_to_ccle = {row["ModelID"]: row["CCLEName"] for row in model_rows if row["CCLEName"]}
    crispr, crispr_included = row_pair_differential(CRISPR_PATH, positive_ids, negative_ids)
    rnai, rnai_included = column_pair_differential(
        RNAI_PATH,
        {model_to_ccle[m] for m in positive_ids if m in model_to_ccle},
        {model_to_ccle[m] for m in negative_ids if m in model_to_ccle},
    )

    analysis = {
        "id": "custom-draft",
        "label": spec["label"],
        "positive_label": spec["positive_label"],
        "negative_label": spec["negative_label"],
        "source": f'Plain-English prompt: "{prompt}"',
        "positive_models": model_payload(model_rows, positive_ids),
        "negative_models": model_payload(model_rows, negative_ids),
        "datasets": {
            "crispr": compact_rows(crispr),
            "rnai": compact_rows(rnai),
        },
        "included_models": {
            "crispr": crispr_included,
            "rnai": rnai_included,
        },
    }
    source = {"prompt": prompt, "spec": spec, "positive_model_ids": sorted(positive_ids), "negative_model_ids": sorted(negative_ids)}
    quality = spec["quality"]
    return analysis, source, quality
