import sys
import unittest
from pathlib import Path
import re


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


class BarAndScatterRenderingTests(unittest.TestCase):
    def test_bar_panel_labels_axes_values_targets_and_exact_tooltips(self):
        panel = visuals.render_bar_panel(
            "Mapped reads",
            [{
                "sample_id": "701_CTCF",
                "assay_target": "CTCF",
                "is_control": False,
                "mapped_percent": 82.9,
            }],
            value_key="mapped_percent",
            axis_label="%",
        )

        self.assertIn('aria-label="Mapped reads (%)"', panel)
        self.assertIn('class="axis-title axis-title-y"', panel)
        self.assertGreaterEqual(panel.count('class="axis-tick-label axis-tick-y"'), 5)
        self.assertIn('class="bar-value"', panel)
        self.assertIn(">82.9%<", panel)
        self.assertIn("<title>701_CTCF: 82.9</title>", panel)
        self.assertIn('data-assay-target="CTCF"', panel)
        self.assertIn('fill="#2F78D1"', panel)

    def test_log10_bar_panel_has_explicit_power_ticks_and_na_for_zero(self):
        panel = visuals.render_bar_panel(
            "Usable fragments after filtering",
            [
                {
                    "sample_id": "large",
                    "assay_target": "RUNX1",
                    "is_control": False,
                    "fragments": 1_000_000,
                },
                {
                    "sample_id": "zero",
                    "assay_target": "IgG",
                    "is_control": True,
                    "fragments": 0,
                },
            ],
            value_key="fragments",
            axis_label="Fragments",
            log10_axis=True,
        )

        self.assertIn("Fragments (log10)", panel)
        self.assertRegex(panel, r">10\^-?\d+(?:\.\d+)?<")
        self.assertIn('data-sample-id="zero"', panel)
        self.assertIn(">NA<", panel)

    def test_log10_bar_panel_keeps_the_smallest_positive_bar_visible(self):
        """Using the smallest observation as baseline would erase a real sample."""
        panel = visuals.render_bar_panel(
            "Usable fragments",
            [
                {
                    "sample_id": "small",
                    "assay_target": "CTCF",
                    "is_control": False,
                    "fragments": 1,
                },
                {
                    "sample_id": "large",
                    "assay_target": "CTCF",
                    "is_control": False,
                    "fragments": 100,
                },
            ],
            value_key="fragments",
            axis_label="Fragments",
            log10_axis=True,
        )

        small = re.search(
            r'class="bar" data-sample-id="small"[^>]*height="([0-9.]+)"',
            panel,
        )
        self.assertIsNotNone(small)
        self.assertGreater(float(small.group(1)), 0)

    def test_range_panel_draws_full_range_iqr_and_median(self):
        panel = visuals.render_range_panel(
            "Peak width median and range",
            [{
                "sample_id": "701_CTCF",
                "assay_target": "CTCF",
                "is_control": False,
                "minimum": 100,
                "q25": 200,
                "median": 300,
                "q75": 400,
                "maximum": 700,
            }],
            minimum_key="minimum",
            q25_key="q25",
            median_key="median",
            q75_key="q75",
            maximum_key="maximum",
            axis_label="Peak width (bp)",
        )

        self.assertIn('aria-label="Peak width median and range (Peak width (bp))"', panel)
        self.assertIn('class="range-min-max"', panel)
        self.assertIn('class="range-iqr"', panel)
        self.assertIn('class="range-median"', panel)
        self.assertIn("<title>701_CTCF: min 100; Q25 200; median 300; Q75 400; max 700</title>", panel)

    def test_scatter_panel_labels_every_point_and_both_axes(self):
        panel = visuals.render_scatter_panel(
            "Peak count vs usable fragments",
            [{
                "sample_id": "701_CTCF",
                "assay_target": "CTCF",
                "is_control": False,
                "usable": 6.17,
                "peaks": 17.2,
            }],
            x_key="usable",
            y_key="peaks",
            x_axis_label="Usable fragments (millions)",
            y_axis_label="Peaks (thousands)",
        )

        self.assertIn('class="axis-title axis-title-x"', panel)
        self.assertIn('class="axis-title axis-title-y"', panel)
        self.assertIn('class="scatter-point"', panel)
        self.assertIn('class="scatter-leader"', panel)
        self.assertIn('class="scatter-label"', panel)
        self.assertIn(">701_CTCF<", panel)
        self.assertIn("<title>701_CTCF: x 6.17; y 17.2</title>", panel)

    def test_empty_panel_is_accessibly_named_and_explicit(self):
        panel = visuals.render_bar_panel(
            "Mapped reads", [], value_key="mapped_percent", axis_label="%"
        )

        self.assertIn('role="status"', panel)
        self.assertIn('aria-label="Mapped reads (%)"', panel)
        self.assertIn("No numeric data available", panel)

    def test_short_labels_stay_horizontal_and_long_labels_rotate(self):
        short = visuals.render_bar_panel(
            "Metric",
            [{"sample_id": "S1", "assay_target": "CTCF", "is_control": False, "value": 1}],
            value_key="value",
            axis_label="Units",
        )
        long = visuals.render_bar_panel(
            "Metric",
            [{
                "sample_id": "long_sample_identifier_replicate_alpha",
                "assay_target": "CTCF",
                "is_control": False,
                "value": 1,
            }],
            value_key="value",
            axis_label="Units",
        )

        self.assertNotIn('class="sample-label rotated"', short)
        self.assertIn('class="sample-label rotated"', long)

    def test_one_hundred_samples_expand_within_horizontal_scroll(self):
        panel = visuals.render_bar_panel(
            "Mapped reads",
            [{
                "sample_id": f"sample_{index:03d}",
                "assay_target": "CTCF",
                "is_control": False,
                "value": index + 1,
            } for index in range(100)],
            value_key="value",
            axis_label="%",
        )

        width = re.search(r'viewBox="0 0 ([0-9.]+) ', panel)
        self.assertIsNotNone(width)
        self.assertGreater(float(width.group(1)), 640)
        self.assertEqual(panel.count('class="bar"'), 100)
        self.assertIn('class="panel-scroll"', panel)

    def test_panel_grid_has_responsive_named_group(self):
        grid = visuals.render_panel_grid(
            ["<article>one</article>", "<article>two</article>"],
            aria_label="Sequencing QC panels",
        )

        self.assertIn('class="panel-grid"', grid)
        self.assertIn('aria-label="Sequencing QC panels"', grid)
        self.assertEqual(grid.count("<article>"), 2)
