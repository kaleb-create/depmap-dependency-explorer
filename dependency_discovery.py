"""High-throughput discovery of selective DepMap dependency stratifiers."""

from __future__ import annotations

import csv
import heapq
import json
import math
import os
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class Stratifier:
    """A positive cohort defined independently of the dependency matrix."""

    key: str
    label: str
    model_ids: frozenset[str]
    category: str = "Uncategorized"
    components: tuple[str, ...] = ()
    aliases: tuple[str, ...] = ()
    combinable: bool = True


@dataclass(frozen=True)
class DiscoveryConfig:
    """Conservative defaults for the rare-stratifier discovery tier."""

    min_prevalence: float = 0.01
    max_prevalence: float = 0.95
    min_screened_positive: int = 15
    min_screened_negative: int = 50
    min_gene_coverage: float = 0.80
    dependency_threshold: float = -0.5
    max_positive_mean: float = -0.5
    min_negative_mean: float = -0.2
    min_dependency_gap: float = 0.4
    max_hedges_g: float = -1.5
    min_positive_dependency_rate: float = 0.60
    max_negative_dependency_rate: float = 0.15
    min_dependency_rate_delta: float = 0.45
    min_gap_lower_bound: float = 0.20
    confidence_z: float = 1.96
    include_pairwise: bool = True
    max_pair_evaluations: int = 2_000_000
    max_pair_stratifiers: int = 10_000
    max_parent_overlap: float = 0.95
    feature_batch_size: int = 256
    gene_batch_size: int = 256
    top_k: int = 500

    def validate(self) -> None:
        if not 0 < self.min_prevalence < self.max_prevalence < 1:
            raise ValueError("Prevalence bounds must satisfy 0 < min < max < 1")
        if not 0 < self.min_gene_coverage <= 1:
            raise ValueError("min_gene_coverage must be in (0, 1]")
        if self.min_screened_positive < 2 or self.min_screened_negative < 2:
            raise ValueError("Both screened cohort minimums must be at least 2")
        if self.top_k < 1:
            raise ValueError("top_k must be positive")


def discover_dependency_targets(
    gene_effect_path: str,
    stratifiers: Iterable[Stratifier],
    population_model_ids: Sequence[str],
    config: DiscoveryConfig | None = None,
) -> dict[str, object]:
    """Scan many stratifier-gene pairs with bounded matrix operations.

    Cohorts are defined without looking at dependency values. Identical cohorts
    are collapsed, pairwise intersections are support-pruned, and sufficient
    statistics are computed with bounded matrix multiplications. The returned
    score uses the lower confidence bound of the dependency gap, which shrinks
    noisy small-cohort hits without imposing a blunt 5% prevalence floor.
    """

    config = config or DiscoveryConfig()
    config.validate()
    population = tuple(dict.fromkeys(population_model_ids))
    if not population:
        raise ValueError("population_model_ids cannot be empty")
    population_set = set(population)

    prepared = _prepare_stratifiers(stratifiers, population_set, config)
    if config.include_pairwise:
        prepared.extend(
            _frequent_pairwise_stratifiers(prepared, population, config)
        )
        prepared = _deduplicate_stratifiers(prepared)

    screened_models, gene_labels, effects = _load_gene_effect_matrix(gene_effect_path)
    screened_index = {model_id: i for i, model_id in enumerate(screened_models)}
    screened_set = set(screened_index)
    prepared = [
        item
        for item in prepared
        if len(item.model_ids & screened_set) >= config.min_screened_positive
        and len(screened_set - item.model_ids) >= config.min_screened_negative
    ]
    single_count = sum(1 for item in prepared if len(item.components) <= 1)
    if not prepared:
        return {
            "hits": [],
            "metadata": {
                "population_models": len(population),
                "screened_models": len(screened_models),
                "single_stratifiers": single_count,
                "tested_stratifiers": 0,
                "genes": len(gene_labels),
                "tested_pairs": 0,
            },
        }

    cohort_matrix = np.zeros((len(prepared), len(screened_models)), dtype=np.float32)
    for feature_index, item in enumerate(prepared):
        indexes = [screened_index[m] for m in item.model_ids if m in screened_index]
        cohort_matrix[feature_index, indexes] = 1.0

    group_sizes = cohort_matrix.sum(axis=1).astype(np.int32)
    negative_group_sizes = len(screened_models) - group_sizes
    prevalence = np.asarray(
        [len(item.model_ids) / len(population) for item in prepared],
        dtype=np.float64,
    )
    heap: list[tuple[float, int, dict[str, object]]] = []
    sequence = 0

    for gene_start in range(0, effects.shape[1], config.gene_batch_size):
        gene_stop = min(gene_start + config.gene_batch_size, effects.shape[1])
        values = effects[:, gene_start:gene_stop]
        finite = np.isfinite(values)
        valid = finite.astype(np.float32)
        clean = np.where(finite, values, 0.0).astype(np.float32, copy=False)
        squares = clean * clean
        dependent = (finite & (values <= config.dependency_threshold)).astype(np.float32)
        total_count = valid.sum(axis=0)
        total_sum = clean.sum(axis=0)
        total_sum_squares = squares.sum(axis=0)
        total_dependent = dependent.sum(axis=0)

        for feature_start in range(0, len(prepared), config.feature_batch_size):
            feature_stop = min(feature_start + config.feature_batch_size, len(prepared))
            cohort = cohort_matrix[feature_start:feature_stop]
            positive_count = cohort @ valid
            positive_sum = cohort @ clean
            positive_sum_squares = cohort @ squares
            positive_dependent = cohort @ dependent
            negative_count = total_count - positive_count
            negative_sum = total_sum - positive_sum
            negative_sum_squares = total_sum_squares - positive_sum_squares
            negative_dependent = total_dependent - positive_dependent

            block_group_sizes = group_sizes[feature_start:feature_stop, None]
            block_negative_sizes = negative_group_sizes[feature_start:feature_stop, None]
            minimum_positive = np.ceil(
                block_group_sizes * config.min_gene_coverage
            )
            minimum_negative = np.ceil(
                block_negative_sizes * config.min_gene_coverage
            )
            sufficient = (
                (positive_count >= minimum_positive)
                & (negative_count >= minimum_negative)
                & (positive_count >= 2)
                & (negative_count >= 2)
            )

            with np.errstate(divide="ignore", invalid="ignore"):
                positive_mean = positive_sum / positive_count
                negative_mean = negative_sum / negative_count
                gap = negative_mean - positive_mean
                positive_variance = (
                    positive_sum_squares - positive_sum * positive_sum / positive_count
                ) / (positive_count - 1)
                negative_variance = (
                    negative_sum_squares - negative_sum * negative_sum / negative_count
                ) / (negative_count - 1)
                positive_variance = np.maximum(positive_variance, 0.0)
                negative_variance = np.maximum(negative_variance, 0.0)
                degrees_freedom = positive_count + negative_count - 2
                pooled_variance = (
                    (positive_count - 1) * positive_variance
                    + (negative_count - 1) * negative_variance
                ) / degrees_freedom
                correction = 1.0 - 3.0 / (4.0 * degrees_freedom - 1.0)
                hedges_g = correction * (positive_mean - negative_mean) / np.sqrt(
                    pooled_variance
                )
                standard_error = np.sqrt(
                    positive_variance / positive_count
                    + negative_variance / negative_count
                )
                gap_lower_bound = gap - config.confidence_z * standard_error
                positive_dependency_rate = positive_dependent / positive_count
                negative_dependency_rate = negative_dependent / negative_count
                dependency_rate_delta = (
                    positive_dependency_rate - negative_dependency_rate
                )

            passing = (
                sufficient
                & np.isfinite(hedges_g)
                & (positive_mean <= config.max_positive_mean)
                & (negative_mean >= config.min_negative_mean)
                & (gap >= config.min_dependency_gap)
                & (hedges_g <= config.max_hedges_g)
                & (positive_dependency_rate >= config.min_positive_dependency_rate)
                & (negative_dependency_rate <= config.max_negative_dependency_rate)
                & (dependency_rate_delta >= config.min_dependency_rate_delta)
                & (gap_lower_bound >= config.min_gap_lower_bound)
            )

            for local_feature, local_gene in np.argwhere(passing):
                feature_index = feature_start + int(local_feature)
                gene_index = gene_start + int(local_gene)
                n_effective = (
                    positive_count[local_feature, local_gene]
                    * negative_count[local_feature, local_gene]
                    / (
                        positive_count[local_feature, local_gene]
                        + negative_count[local_feature, local_gene]
                    )
                )
                stability_weight = math.sqrt(float(n_effective) / (float(n_effective) + 20.0))
                score = float(
                    gap_lower_bound[local_feature, local_gene]
                    * dependency_rate_delta[local_feature, local_gene]
                    * -hedges_g[local_feature, local_gene]
                    * stability_weight
                )
                item = prepared[feature_index]
                hit = {
                    "stratifier_id": item.key,
                    "stratifier": item.label,
                    "category": item.category,
                    "components": list(item.components),
                    "aliases": list(item.aliases),
                    "prevalence": round(float(prevalence[feature_index]), 6),
                    "positive_models": int(group_sizes[feature_index]),
                    "negative_models": int(negative_group_sizes[feature_index]),
                    "gene": _gene_symbol(gene_labels[gene_index]),
                    "gene_label": gene_labels[gene_index],
                    "positive_mean": round(float(positive_mean[local_feature, local_gene]), 6),
                    "negative_mean": round(float(negative_mean[local_feature, local_gene]), 6),
                    "dependency_gap": round(float(gap[local_feature, local_gene]), 6),
                    "gap_lower_bound": round(
                        float(gap_lower_bound[local_feature, local_gene]), 6
                    ),
                    "hedges_g": round(float(hedges_g[local_feature, local_gene]), 6),
                    "positive_dependency_rate": round(
                        float(positive_dependency_rate[local_feature, local_gene]), 6
                    ),
                    "negative_dependency_rate": round(
                        float(negative_dependency_rate[local_feature, local_gene]), 6
                    ),
                    "positive_n": int(positive_count[local_feature, local_gene]),
                    "negative_n": int(negative_count[local_feature, local_gene]),
                    "discovery_score": round(score, 6),
                }
                sequence += 1
                entry = (score, sequence, hit)
                if len(heap) < config.top_k:
                    heapq.heappush(heap, entry)
                elif score > heap[0][0]:
                    heapq.heapreplace(heap, entry)

    hits = [entry[2] for entry in sorted(heap, reverse=True)]
    return {
        "hits": hits,
        "metadata": {
            "population_models": len(population),
            "screened_models": len(screened_models),
            "single_stratifiers": single_count,
            "tested_stratifiers": len(prepared),
            "pairwise_stratifiers": len(prepared) - single_count,
            "genes": len(gene_labels),
            "tested_pairs": len(prepared) * len(gene_labels),
            "returned_hits": len(hits),
            "minimum_prevalence": config.min_prevalence,
        },
    }


def stratifiers_from_summary(summary_path: str) -> list[Stratifier]:
    """Load all currently curated explorer cohorts as discovery features."""

    with open(summary_path) as handle:
        payload = json.load(handle)
    features: list[Stratifier] = []
    for analysis in payload.get("analyses", []):
        model_ids = frozenset(
            model["model_id"]
            for model in analysis.get("positive_models", [])
            if model.get("model_id")
        )
        if model_ids:
            features.append(
                Stratifier(
                    key=str(analysis["id"]),
                    label=str(analysis["label"]),
                    model_ids=model_ids,
                    category=str(analysis.get("category") or "Curated"),
                    components=(str(analysis["id"]),),
                )
            )
    return features


def model_population(model_path: str) -> list[str]:
    """Return the complete model universe used for prevalence calculations."""

    with open(model_path, newline="") as handle:
        return [row["ModelID"] for row in csv.DictReader(handle) if row.get("ModelID")]


def _prepare_stratifiers(
    stratifiers: Iterable[Stratifier],
    population: set[str],
    config: DiscoveryConfig,
) -> list[Stratifier]:
    minimum = math.ceil(len(population) * config.min_prevalence)
    maximum = math.floor(len(population) * config.max_prevalence)
    normalized = []
    for item in stratifiers:
        models = frozenset(item.model_ids & population)
        if minimum <= len(models) <= maximum:
            normalized.append(
                Stratifier(
                    key=item.key,
                    label=item.label,
                    model_ids=models,
                    category=item.category,
                    components=item.components or (item.key,),
                    aliases=item.aliases,
                    combinable=item.combinable,
                )
            )
    return _deduplicate_stratifiers(normalized)


def _deduplicate_stratifiers(stratifiers: Iterable[Stratifier]) -> list[Stratifier]:
    by_cohort: dict[frozenset[str], list[Stratifier]] = {}
    for item in stratifiers:
        by_cohort.setdefault(item.model_ids, []).append(item)
    deduplicated = []
    for cohort, matches in by_cohort.items():
        matches.sort(key=lambda item: (len(item.components), item.key))
        canonical = matches[0]
        aliases = tuple(
            dict.fromkeys(
                alias
                for item in matches
                for alias in (item.key, *item.aliases)
                if alias != canonical.key
            )
        )
        deduplicated.append(
            Stratifier(
                key=canonical.key,
                label=canonical.label,
                model_ids=cohort,
                category=canonical.category,
                components=canonical.components,
                aliases=aliases,
                combinable=any(item.combinable for item in matches),
            )
        )
    deduplicated.sort(key=lambda item: item.key)
    return deduplicated


def _frequent_pairwise_stratifiers(
    stratifiers: Sequence[Stratifier],
    population: Sequence[str],
    config: DiscoveryConfig,
) -> list[Stratifier]:
    """Generate frequent intersections with bitsets and Apriori support pruning."""

    model_index = {model_id: i for i, model_id in enumerate(population)}
    minimum = math.ceil(len(population) * config.min_prevalence)
    maximum = math.floor(len(population) * config.max_prevalence)
    candidates = [item for item in stratifiers if item.combinable]
    bitsets = [_models_to_bits(item.model_ids, model_index) for item in candidates]
    seen = {_models_to_bits(item.model_ids, model_index) for item in stratifiers}
    combinations = []
    evaluations = 0
    stop = False
    for left_index, left in enumerate(candidates):
        if stop:
            break
        left_bits = bitsets[left_index]
        left_size = left_bits.bit_count()
        for right_index in range(left_index + 1, len(candidates)):
            evaluations += 1
            if evaluations > config.max_pair_evaluations:
                stop = True
                break
            right = candidates[right_index]
            cohort_bits = left_bits & bitsets[right_index]
            size = cohort_bits.bit_count()
            if size < minimum or size > maximum or cohort_bits in seen:
                continue
            if size / min(left_size, bitsets[right_index].bit_count()) > config.max_parent_overlap:
                continue
            seen.add(cohort_bits)
            components = tuple(dict.fromkeys((*left.components, *right.components)))
            combinations.append(
                Stratifier(
                    key="and:" + "&".join(components),
                    label=f"{left.label} AND {right.label}",
                    model_ids=_bits_to_models(cohort_bits, population),
                    category="Discovered combinations",
                    components=components,
                )
            )
            if len(combinations) >= config.max_pair_stratifiers:
                stop = True
                break
    return combinations


def _models_to_bits(models: Iterable[str], model_index: dict[str, int]) -> int:
    bits = 0
    for model in models:
        index = model_index.get(model)
        if index is not None:
            bits |= 1 << index
    return bits


def _bits_to_models(bits: int, population: Sequence[str]) -> frozenset[str]:
    return frozenset(
        model for index, model in enumerate(population) if bits & (1 << index)
    )


def _load_gene_effect_matrix(path: str) -> tuple[list[str], list[str], np.ndarray]:
    model_ids = []
    rows = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        gene_labels = header[1:]
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Malformed gene-effect row {row_number} in {path}")
            model_ids.append(row[0])
            rows.append(
                np.fromiter(
                    (_float_or_nan(value) for value in row[1:]),
                    dtype=np.float32,
                    count=len(gene_labels),
                )
            )
    if len(model_ids) != len(set(model_ids)):
        raise ValueError(f"Duplicate model identifiers in {path}")
    matrix = np.vstack(rows) if rows else np.empty((0, len(gene_labels)), dtype=np.float32)
    return model_ids, gene_labels, matrix


def _float_or_nan(value: str) -> float:
    if not value or value.upper() == "NA":
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def _gene_symbol(label: str) -> str:
    return label.rsplit(" (", 1)[0]
