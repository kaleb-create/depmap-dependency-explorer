import csv
import hashlib
import ipaddress
import json
import math
import os
import socket
import urllib.parse
import urllib.request
from typing import Any

from scripts.build_hpv_dependency_data import (
    compact_rows,
    hedges_g_from_moments,
    parse_gene,
)


BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.environ.get("DEPMAP_DATA_DIR", os.path.join(BASE_DIR, "data"))
CRISPR_PATH = os.path.join(DATA_DIR, "CRISPRGeneEffect.csv")
RNAI_PATH = os.path.join(DATA_DIR, "D2_combined_gene_dep_scores.csv")
MODEL_PATH = os.path.join(DATA_DIR, "Model.csv")

OPENAI_API_URL = "https://api.openai.com/v1/responses"
OPENAI_MODEL = os.environ.get("OPENAI_MODEL", "gpt-5.5")
STRATIFIER_DATASET_DIR = os.environ.get(
    "STRATIFIER_DATASET_DIR", os.path.join(DATA_DIR, "stratifier_sources")
)
MAX_DATASET_BYTES = int(os.environ.get("MAX_STRATIFIER_DATASET_BYTES", 25 * 1024 * 1024))
MAX_DATASET_ROWS = int(os.environ.get("MAX_STRATIFIER_DATASET_ROWS", 250_000))
MIN_GROUP_VALUES = 3
MIN_GROUP_COVERAGE = 0.5

IDENTIFIER_COLUMNS = (
    "ModelID",
    "DepMap_ID",
    "DepMapID",
    "model_id",
    "CCLEName",
    "CCLE_Name",
    "cell_line_name",
    "CellLineName",
    "StrippedCellLineName",
    "RRID",
    "COSMICID",
)


def normalized(value: object) -> str:
    return "".join(character for character in str(value or "").casefold() if character.isalnum())


def validate_public_url(url: str) -> None:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme not in {"http", "https"} or not parsed.hostname:
        raise ValueError("Dataset download URL must be a public HTTP or HTTPS address.")

    try:
        addresses = socket.getaddrinfo(parsed.hostname, parsed.port or 443, type=socket.SOCK_STREAM)
    except socket.gaierror as exc:
        raise ValueError(f"Could not resolve dataset host {parsed.hostname}.") from exc

    for address in addresses:
        ip = ipaddress.ip_address(address[4][0])
        if not ip.is_global:
            raise ValueError("Dataset download URL resolves to a non-public network address.")


class PublicRedirectHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        validate_public_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


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
        positive_sum_squares = [0.0] * len(labels)
        positive_counts = [0] * len(labels)
        negative_sums = [0.0] * len(labels)
        negative_sum_squares = [0.0] * len(labels)
        negative_counts = [0] * len(labels)
        positive_rows = []
        negative_rows = []

        for row in reader:
            if not row:
                continue
            if row[0] in positive_ids:
                sums = positive_sums
                sum_squares = positive_sum_squares
                counts = positive_counts
                positive_rows.append(row[0])
            elif row[0] in negative_ids:
                sums = negative_sums
                sum_squares = negative_sum_squares
                counts = negative_counts
                negative_rows.append(row[0])
            else:
                continue

            for i, value in enumerate(row[1:]):
                if value == "" or value.upper() == "NA":
                    continue
                try:
                    parsed = float(value)
                    sums[i] += parsed
                    sum_squares[i] += parsed * parsed
                    counts[i] += 1
                except ValueError:
                    continue

    rows = []
    min_positive_n = min_values_for_group(len(positive_rows))
    min_negative_n = min_values_for_group(len(negative_rows))
    for i, label in enumerate(labels):
        if positive_counts[i] < min_positive_n or negative_counts[i] < min_negative_n:
            continue
        positive_average = positive_sums[i] / positive_counts[i]
        negative_average = negative_sums[i] / negative_counts[i]
        raw_difference = positive_average - negative_average
        effect_size = hedges_g_from_moments(
            positive_sums[i],
            positive_sum_squares[i],
            positive_counts[i],
            negative_sums[i],
            negative_sum_squares[i],
            negative_counts[i],
        )
        if effect_size is None:
            continue
        gene, _ = parse_gene(label)
        rows.append(
            {
                "gene": gene,
                "score": round(effect_size, 6),
                "raw_difference": round(raw_difference, 6),
                "positive_average": round(positive_average, 6),
                "negative_average": round(negative_average, 6),
                "positive_n": positive_counts[i],
                "negative_n": negative_counts[i],
            }
        )
    add_ranks(rows)
    return rows, {
        "positive": sorted(positive_rows),
        "negative_n": len(negative_rows),
        "min_positive_n": min_positive_n,
        "min_negative_n": min_negative_n,
    }


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

        min_positive_n = min_values_for_group(len(positive_indexes))
        min_negative_n = min_values_for_group(len(negative_indexes))
        rows = []
        for row in reader:
            positive_values = values_at(row, positive_indexes)
            negative_values = values_at(row, negative_indexes)
            if len(positive_values) < min_positive_n or len(negative_values) < min_negative_n:
                continue
            positive_sum = sum(positive_values)
            positive_sum_squares = sum(value * value for value in positive_values)
            negative_sum = sum(negative_values)
            negative_sum_squares = sum(value * value for value in negative_values)
            positive_average = positive_sum / len(positive_values)
            negative_average = negative_sum / len(negative_values)
            raw_difference = positive_average - negative_average
            effect_size = hedges_g_from_moments(
                positive_sum,
                positive_sum_squares,
                len(positive_values),
                negative_sum,
                negative_sum_squares,
                len(negative_values),
            )
            if effect_size is None:
                continue
            gene, _ = parse_gene(row[0])
            rows.append(
                {
                    "gene": gene,
                    "score": round(effect_size, 6),
                    "raw_difference": round(raw_difference, 6),
                    "positive_average": round(positive_average, 6),
                    "negative_average": round(negative_average, 6),
                    "positive_n": len(positive_values),
                    "negative_n": len(negative_values),
                }
            )
    add_ranks(rows)
    return rows, {
        "positive": sorted(positive_used),
        "negative_n": len(negative_used),
        "min_positive_n": min_positive_n,
        "min_negative_n": min_negative_n,
    }


def min_values_for_group(group_size: int) -> int:
    return max(MIN_GROUP_VALUES, math.ceil(max(0, group_size) * MIN_GROUP_COVERAGE))


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


def dataset_extension(dataset_format: str) -> str:
    return {"csv": ".csv", "tsv": ".tsv", "json": ".json"}.get(dataset_format, ".data")


def retrieve_dataset(dataset: dict[str, Any]) -> dict[str, Any]:
    dataset_format = dataset.get("format") or "metadata_only"
    download_url = (dataset.get("download_url") or "").strip()
    result: dict[str, Any] = {
        "status": "local_metadata" if dataset_format == "metadata_only" else "not_downloaded",
        "format": dataset_format,
        "download_url": download_url,
    }
    if dataset_format == "metadata_only":
        return result
    if dataset_format not in {"csv", "tsv", "json"}:
        raise ValueError(f"Unsupported dataset format: {dataset_format}.")
    if not download_url:
        raise ValueError("The selected external dataset did not include a direct download URL.")

    validate_public_url(download_url)
    request = urllib.request.Request(
        download_url,
        headers={"User-Agent": "DepMap-Dependency-Explorer/1.0"},
    )
    opener = urllib.request.build_opener(PublicRedirectHandler())
    with opener.open(request, timeout=45) as response:
        content_length = int(response.headers.get("Content-Length") or 0)
        if content_length > MAX_DATASET_BYTES:
            raise ValueError(
                f"Selected dataset is larger than the {MAX_DATASET_BYTES // (1024 * 1024)} MB download limit."
            )
        payload = response.read(MAX_DATASET_BYTES + 1)
        final_url = response.geturl()
        content_type = response.headers.get("Content-Type", "").split(";", 1)[0]

    if len(payload) > MAX_DATASET_BYTES:
        raise ValueError(
            f"Selected dataset is larger than the {MAX_DATASET_BYTES // (1024 * 1024)} MB download limit."
        )
    if payload.lstrip().lower().startswith((b"<!doctype html", b"<html")):
        raise ValueError("The dataset download URL returned a web page instead of a data file.")

    os.makedirs(STRATIFIER_DATASET_DIR, exist_ok=True)
    digest = hashlib.sha256(payload).hexdigest()
    path = os.path.join(STRATIFIER_DATASET_DIR, f"{digest}{dataset_extension(dataset_format)}")
    if not os.path.exists(path):
        with open(path, "wb") as f:
            f.write(payload)

    result.update(
        {
            "status": "downloaded",
            "final_url": final_url,
            "content_type": content_type,
            "bytes": len(payload),
            "sha256": digest,
            "cached_path": path,
        }
    )
    return result


def load_external_rows(path: str, dataset_format: str) -> list[dict[str, Any]]:
    if dataset_format in {"csv", "tsv"}:
        with open(path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="," if dataset_format == "csv" else "\t")
            rows = []
            for row in reader:
                rows.append(dict(row))
                if len(rows) >= MAX_DATASET_ROWS:
                    raise ValueError(f"Dataset exceeds the {MAX_DATASET_ROWS:,}-row processing limit.")
            return rows

    with open(path, encoding="utf-8-sig") as f:
        payload = json.load(f)
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict):
        rows = next(
            (payload[key] for key in ("data", "records", "rows") if isinstance(payload.get(key), list)),
            None,
        )
    else:
        rows = None
    if rows is None or any(not isinstance(row, dict) for row in rows):
        raise ValueError("JSON dataset must be an array of objects or contain a data, records, or rows array.")
    if len(rows) > MAX_DATASET_ROWS:
        raise ValueError(f"Dataset exceeds the {MAX_DATASET_ROWS:,}-row processing limit.")
    return rows


def find_column(rows: list[dict[str, Any]], requested: str, candidates: tuple[str, ...] = ()) -> str:
    if not rows:
        raise ValueError("The selected dataset contains no rows.")
    columns = list(rows[0])
    by_normalized = {normalized(column): column for column in columns}
    for candidate in (requested, *candidates):
        if normalized(candidate) in by_normalized:
            return by_normalized[normalized(candidate)]
    raise ValueError(f"Could not find required column '{requested}' in the selected dataset.")


def model_identifier_lookup(model_rows: list[dict[str, str]]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    fields = (
        "ModelID",
        "ModelIDAlias",
        "CCLEName",
        "CellLineName",
        "StrippedCellLineName",
        "RRID",
        "COSMICID",
        "SangerModelID",
    )
    for row in model_rows:
        model_id = row["ModelID"]
        for field in fields:
            raw = row.get(field) or ""
            values = raw.replace(";", ",").split(",") if field == "ModelIDAlias" else [raw]
            for value in values:
                key = normalized(value)
                if key:
                    lookup[key] = model_id
    return lookup


def map_external_dataset(
    model_rows: list[dict[str, str]], dataset: dict[str, Any], retrieval: dict[str, Any]
) -> tuple[set[str], set[str], dict[str, Any]]:
    rows = load_external_rows(retrieval["cached_path"], retrieval["format"])
    identifier_column = find_column(rows, dataset.get("identifier_column") or "", IDENTIFIER_COLUMNS)
    group_column = find_column(rows, dataset.get("group_column") or "")
    positive_values = {normalized(value) for value in dataset.get("positive_values") or [] if normalized(value)}
    negative_values = {normalized(value) for value in dataset.get("negative_values") or [] if normalized(value)}
    if not positive_values:
        raise ValueError("The selected dataset did not define any positive-group values.")

    lookup = model_identifier_lookup(model_rows)
    positive_ids: set[str] = set()
    negative_ids: set[str] = set()
    mapped_rows = 0
    unmatched_identifiers: set[str] = set()
    for row in rows:
        raw_identifier = row.get(identifier_column)
        model_id = lookup.get(normalized(raw_identifier))
        if not model_id:
            if raw_identifier:
                unmatched_identifiers.add(str(raw_identifier))
            continue
        mapped_rows += 1
        group_value = normalized(row.get(group_column))
        if group_value in positive_values:
            positive_ids.add(model_id)
        if group_value in negative_values:
            negative_ids.add(model_id)

    if not negative_values:
        negative_ids = {row["ModelID"] for row in model_rows} - positive_ids
    negative_ids -= positive_ids
    mapping = {
        "status": "mapped",
        "rows": len(rows),
        "mapped_rows": mapped_rows,
        "unmatched_identifier_count": len(unmatched_identifiers),
        "identifier_column": identifier_column,
        "group_column": group_column,
        "positive_models": len(positive_ids),
        "negative_models": len(negative_ids),
    }
    return positive_ids, negative_ids, mapping


def response_sources(payload: dict[str, Any]) -> list[dict[str, str]]:
    sources: dict[str, str] = {}
    for item in payload.get("output", []):
        action = item.get("action") or {}
        for source in action.get("sources") or []:
            if source.get("url"):
                sources[source["url"]] = source.get("title") or source["url"]
        for content in item.get("content") or []:
            for annotation in content.get("annotations") or []:
                if annotation.get("type") == "url_citation" and annotation.get("url"):
                    sources[annotation["url"]] = annotation.get("title") or annotation["url"]
    return [{"title": title, "url": url} for url, title in sources.items()]


def request_openai_spec(prompt: str, taxonomy: dict[str, list[str]]) -> dict[str, Any]:
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY is not configured on this server.")

    schema = {
        "type": "object",
        "additionalProperties": False,
        "required": ["label", "positive_label", "negative_label", "positive", "negative", "dataset", "quality"],
        "properties": {
            "label": {"type": "string"},
            "positive_label": {"type": "string"},
            "negative_label": {"type": "string"},
            "positive": {"$ref": "#/$defs/group"},
            "negative": {"$ref": "#/$defs/group"},
            "dataset": {
                "type": "object",
                "additionalProperties": False,
                "required": [
                    "name",
                    "provider",
                    "source_url",
                    "download_url",
                    "format",
                    "identifier_column",
                    "group_column",
                    "positive_values",
                    "negative_values",
                    "why_selected",
                    "mapping_notes"
                ],
                "properties": {
                    "name": {"type": "string"},
                    "provider": {"type": "string"},
                    "source_url": {"type": "string"},
                    "download_url": {"type": "string"},
                    "format": {"type": "string", "enum": ["csv", "tsv", "json", "metadata_only"]},
                    "identifier_column": {"type": "string"},
                    "group_column": {"type": "string"},
                    "positive_values": {"type": "array", "items": {"type": "string"}},
                    "negative_values": {"type": "array", "items": {"type": "string"}},
                    "why_selected": {"type": "string"},
                    "mapping_notes": {"type": "string"}
                }
            },
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
                    "Build a reproducible stratifier for a DepMap dependency analysis. Search the web for "
                    "the most authoritative, current dataset that defines both sides of the requested cell-line "
                    "contrast and can be mapped deterministically to DepMap ModelID, CCLEName, cell-line name, "
                    "Cellosaurus RRID, or COSMIC ID. Prefer primary maintained sources such as DepMap, "
                    "Cellosaurus, NCI, cBioPortal, COSMIC, or a peer-reviewed data repository. For a downloadable "
                    "CSV, TSV, or JSON table, provide a direct data-file URL, the exact identifier and grouping "
                    "columns, and exact positive/negative values. Leave negative_values empty only for a true "
                    "positive-vs-all-other comparison. Use format metadata_only only when the locally installed "
                    "DepMap Model.csv taxonomy is itself the best cohort-defining dataset; then provide strict "
                    "filters over the supplied taxonomy. For external datasets, positive and negative metadata "
                    "filters should express requested lineage, disease, or subtype restrictions shared with the "
                    "molecular grouping, not mutation or fusion terms absent from Model.csv. Do not invent URLs, "
                    "columns, values, identifiers, or "
                    "dataset capabilities. Critically state weaknesses and validation checks."
                ),
            },
            {
                "role": "user",
                "content": json.dumps({"prompt": prompt, "available_taxonomy": taxonomy}),
            },
        ],
        "tools": [{"type": "web_search"}],
        "tool_choice": "auto",
        "include": ["web_search_call.action.sources"],
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
                spec = json.loads(content["text"])
                spec["_web_sources"] = response_sources(payload)
                return spec
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
    dataset = spec["dataset"]
    retrieval = retrieve_dataset(dataset)
    if retrieval["status"] == "downloaded":
        positive_ids, negative_ids, mapping = map_external_dataset(model_rows, dataset, retrieval)
        positive_scope = choose_models(model_rows, spec["positive"])
        negative_scope = choose_models(model_rows, spec["negative"])
        if positive_scope:
            positive_ids &= positive_scope
        if negative_scope:
            if dataset.get("negative_values"):
                negative_ids &= negative_scope
            else:
                negative_ids = negative_scope - positive_ids
        mapping["positive_models"] = len(positive_ids)
        mapping["negative_models"] = len(negative_ids)
        mapping["metadata_scope_applied"] = bool(positive_scope or negative_scope)
        cohort_method = "external_dataset"
    else:
        positive_ids = choose_models(model_rows, spec["positive"])
        negative_ids = choose_models(model_rows, spec["negative"])
        if not negative_ids:
            negative_ids = {row["ModelID"] for row in model_rows} - positive_ids
        negative_ids -= positive_ids
        mapping = {
            "status": "local_metadata",
            "rows": len(model_rows),
            "mapped_rows": len(model_rows),
            "identifier_column": "ModelID",
            "group_column": dataset.get("group_column") or "DepMap model metadata",
            "positive_models": len(positive_ids),
            "negative_models": len(negative_ids),
        }
        cohort_method = "depmap_metadata"

    if len(positive_ids) < MIN_GROUP_VALUES or len(negative_ids) < MIN_GROUP_VALUES:
        raise ValueError(
            f"The selected dataset mapped only {len(positive_ids)} positive models and "
            f"{len(negative_ids)} negative models; at least {MIN_GROUP_VALUES} are required on each side."
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
        "source": f'{dataset["name"]} ({dataset["provider"]})',
        "effect_metric": "hedges_g",
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
    retrieval_for_source = dict(retrieval)
    if retrieval_for_source.get("cached_path"):
        retrieval_for_source["cached_file"] = os.path.basename(retrieval_for_source.pop("cached_path"))
    source = {
        "prompt": prompt,
        "dataset": dataset,
        "retrieval": retrieval_for_source,
        "mapping": mapping,
        "cohort_method": cohort_method,
        "web_sources": spec.pop("_web_sources", []),
        "spec": spec,
        "positive_model_ids": sorted(positive_ids),
        "negative_model_ids": sorted(negative_ids),
    }
    quality = spec["quality"]
    quality["dataset_summary"] = (
        f'{dataset["name"]} mapped {mapping["positive_models"]} positive and '
        f'{mapping["negative_models"]} negative DepMap models.'
    )
    return analysis, source, quality
