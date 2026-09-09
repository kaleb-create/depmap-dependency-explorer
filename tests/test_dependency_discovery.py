import csv
import os
import tempfile
import unittest

from dependency_discovery import DiscoveryConfig, Stratifier, discover_dependency_targets


class DependencyDiscoveryTests(unittest.TestCase):
    def test_finds_selective_dependency_and_support_pruned_combination(self):
        population = [f"M{i}" for i in range(20)]
        feature_a = Stratifier(
            key="a",
            label="Feature A",
            model_ids=frozenset(population[:8]),
            components=("a",),
        )
        feature_b = Stratifier(
            key="b",
            label="Feature B",
            model_ids=frozenset(population[4:12]),
            components=("b",),
        )
        values = []
        for index in range(20):
            if 4 <= index < 8:
                selective = -1.0 + (index - 5.5) * 0.02
            else:
                selective = ((index % 3) - 1) * 0.01
            panessential = -1.0 + ((index % 5) - 2) * 0.01
            values.append([selective, panessential])

        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "effects.csv")
            with open(path, "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["ModelID", "SELECTIVE (1)", "PANESSENTIAL (2)"])
                for model, row in zip(population, values):
                    writer.writerow([model, *row])

            result = discover_dependency_targets(
                path,
                [feature_a, feature_b],
                population,
                DiscoveryConfig(
                    min_prevalence=0.10,
                    min_screened_positive=3,
                    min_screened_negative=3,
                    min_gene_coverage=1.0,
                    max_positive_mean=-0.8,
                    min_negative_mean=-0.1,
                    min_dependency_gap=0.8,
                    max_hedges_g=-1.0,
                    min_positive_dependency_rate=0.9,
                    max_negative_dependency_rate=0.1,
                    min_dependency_rate_delta=0.8,
                    min_gap_lower_bound=0.5,
                    max_parent_overlap=0.9,
                    top_k=10,
                ),
            )

        self.assertEqual(result["metadata"]["pairwise_stratifiers"], 1)
        self.assertEqual(len(result["hits"]), 1)
        hit = result["hits"][0]
        self.assertEqual(hit["gene"], "SELECTIVE")
        self.assertEqual(hit["components"], ["a", "b"])
        self.assertEqual(hit["positive_models"], 4)
        self.assertEqual(hit["positive_dependency_rate"], 1.0)
        self.assertEqual(hit["negative_dependency_rate"], 0.0)


if __name__ == "__main__":
    unittest.main()
