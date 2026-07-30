import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

import qc_dashboard_visuals as visuals


class VisualFormattingTests(unittest.TestCase):
    def test_target_colors_are_stable_and_controls_are_grey(self):
        self.assertEqual(visuals.target_color("IgG", is_control=True), "#8F8D87")
        self.assertEqual(visuals.target_color("CTCF"), "#2F78D1")
        self.assertEqual(visuals.target_color("GATA1"), "#1FAE7A")
        self.assertEqual(visuals.target_color("RUNX1"), "#F06432")
        self.assertEqual(visuals.target_color("NEW_TF"), visuals.target_color("NEW_TF"))

    def test_html_numbers_use_three_significant_digits_without_rounding_exact_fields(self):
        self.assertEqual(visuals.format_significant(17_432_100), "17.4M")
        self.assertEqual(visuals.format_significant(0.012345), "0.0123")
        self.assertEqual(visuals.format_significant(0.000012345), "1.23e-05")
        self.assertEqual(visuals.format_significant(1250, exact=True), "1250")

    def test_html_numbers_promote_compact_suffixes_after_rounding(self):
        self.assertEqual(visuals.format_significant(999.9), "1.00K")
        self.assertEqual(visuals.format_significant(999_500), "1.00M")
        self.assertEqual(visuals.format_significant(999_500_000), "1.00B")


class DistributionTransformTests(unittest.TestCase):
    def test_shared_250bp_bins_include_zero_bins_and_normalize_each_sample(self):
        result = visuals.bin_weighted_series({
            "A": [(0, 1), (249, 1), (250, 2), (501, 1)],
            "B": [(260, 4)],
        })
        self.assertEqual(
            [(row["bin_start"], row["bin_end"]) for row in result["A"]],
            [(0, 250), (250, 500), (500, 750)],
        )
        self.assertEqual(
            [row["count"] for row in result["B"]],
            [0, 4, 0],
        )
        self.assertAlmostEqual(sum(row["percent"] for row in result["A"]), 100.0)

    def test_fragment_histogram_ecdf_retains_zero_and_reaches_one_hundred_percent(self):
        rows = visuals.histogram_ecdf([(0, 2), (1, 1), (10, 1)])
        self.assertEqual(rows[0], {"value": 0, "count": 2, "cumulative_percent": 50.0})
        self.assertEqual(rows[-1]["cumulative_percent"], 100.0)

    def test_distribution_transforms_reject_invalid_inputs_and_combine_duplicate_values(self):
        with self.assertRaisesRegex(ValueError, "bin_size"):
            visuals.bin_weighted_series({"A": [(0, 1)]}, bin_size=0)
        with self.assertRaisesRegex(ValueError, "position"):
            visuals.bin_weighted_series({"A": [(-1, 1)]})
        self.assertEqual(
            visuals.histogram_ecdf([(2, 1), (2, 2), (3, 0)]),
            [{"value": 2, "count": 3, "cumulative_percent": 100.0}],
        )


class EndpointLabelPackingTests(unittest.TestCase):
    def test_endpoint_labels_are_separated_and_constrained_to_the_plot(self):
        packed = visuals.pack_endpoint_labels(
            [("first", 9.0), ("second", 9.1), ("third", 9.2)],
            lower=0.0,
            upper=10.0,
            minimum_gap=1.0,
        )
        self.assertEqual(packed, {"first": 8.0, "second": 9.0, "third": 10.0})

    def test_endpoint_label_packing_rejects_an_impossible_layout(self):
        with self.assertRaisesRegex(ValueError, "impossible"):
            visuals.pack_endpoint_labels(
                [("first", 0.0), ("second", 1.0), ("third", 2.0)],
                lower=0.0,
                upper=1.0,
                minimum_gap=1.0,
            )
