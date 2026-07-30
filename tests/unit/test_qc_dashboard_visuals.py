import sys
import unittest
from pathlib import Path
import re
import math


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


class DistributionRenderingTests(unittest.TestCase):
    def metadata(self):
        return {
            "S1": {
                "assay_target": "CTCF",
                "is_control": False,
            },
            "S2": {
                "assay_target": "CTCF",
                "is_control": False,
            },
            "I1": {
                "assay_target": "IgG",
                "is_control": True,
            },
        }

    def test_binned_distribution_uses_midpoints_percent_and_target_styles(self):
        panel = visuals.render_binned_distribution(
            "Insert-size distribution",
            {
                "S1": [
                    {"bin_start": 0, "bin_end": 250, "count": 3, "percent": 75.0},
                    {"bin_start": 250, "bin_end": 500, "count": 1, "percent": 25.0},
                ],
                "S2": [
                    {"bin_start": 0, "bin_end": 250, "count": 1, "percent": 25.0},
                    {"bin_start": 250, "bin_end": 500, "count": 3, "percent": 75.0},
                ],
            },
            self.metadata(),
            x_axis_label="Insert size (bp)",
        )

        self.assertRegex(
            panel,
            r'class="axis-title axis-title-y"[^>]*>Percent of read pairs</text>',
        )
        self.assertIn('data-bin-midpoint="125"', panel)
        self.assertIn('data-bin-midpoint="375"', panel)
        last_midpoint = re.search(
            r'data-bin-midpoint="375"[^>]*(?:cx|x)="([0-9.]+)"',
            panel,
        )
        self.assertIsNotNone(last_midpoint)
        self.assertLess(float(last_midpoint.group(1)), 535.0)
        self.assertEqual(panel.count('data-assay-target="CTCF"'), 2)
        self.assertEqual(
            len(re.findall(r'class="distribution-trace"[^>]*stroke="#2F78D1"', panel)),
            2,
        )
        dashes = re.findall(
            r'data-sample-id="S[12]"[^>]*stroke-dasharray="([^"]+)"',
            panel,
        )
        self.assertEqual(len(set(dashes)), 2)
        self.assertEqual(panel.count('class="endpoint-label"'), 2)
        self.assertIn("<title>S1: [0, 250) bp; 75%; count 3</title>", panel)

    def test_endpoint_labels_are_packed_inside_plot_and_get_leaders(self):
        panel = visuals.render_binned_distribution(
            "Peak-width distribution",
            {
                sample_id: [
                    {"bin_start": 0, "bin_end": 250, "count": 1, "percent": 50.0},
                    {"bin_start": 250, "bin_end": 500, "count": 1, "percent": 50.0},
                ]
                for sample_id in ("S1", "S2", "I1")
            },
            self.metadata(),
            x_axis_label="Peak width (bp)",
        )

        positions = [
            float(value)
            for value in re.findall(r'class="endpoint-label"[^>]* y="([0-9.]+)"', panel)
        ]
        self.assertEqual(len(positions), 3)
        self.assertTrue(all(36.0 <= value <= 236.0 for value in positions))
        self.assertGreaterEqual(
            min(b - a for a, b in zip(sorted(positions), sorted(positions)[1:])),
            13.0,
        )
        self.assertIn('class="endpoint-leader"', panel)

    def test_fragments_per_peak_ecdf_has_separate_zero_and_log_ticks(self):
        panel = visuals.render_ecdf(
            "Fragments per peak",
            {
                "S1": [
                    {"value": 0, "count": 2, "cumulative_percent": 40.0},
                    {"value": 1, "count": 1, "cumulative_percent": 60.0},
                    {"value": 10, "count": 1, "cumulative_percent": 80.0},
                    {"value": 100, "count": 1, "cumulative_percent": 100.0},
                ],
            },
            self.metadata(),
            x_axis_label="Fragments per peak",
            zero_origin=True,
        )

        self.assertIn('data-zero-origin="true"', panel)
        for label in ("0", "1", "10", "100"):
            self.assertRegex(panel, rf'class="axis-tick-label axis-tick-x"[^>]*>{label}<')
        self.assertIn("Cumulative percent of peaks", panel)
        self.assertIn("<title>S1: 0 fragments; 40% cumulative; count 2</title>", panel)

    def test_profile_chart_uses_target_color_tss_reference_and_direct_labels(self):
        panel = visuals.render_profile_chart(
            "TSS profiles",
            {
                "S1": [(-10, 1.0), (0, 4.0), (10, 2.0)],
                "S2": [(-10, 2.0), (0, 5.0), (10, 2.0)],
            },
            self.metadata(),
            x_axis_label="Position relative to TSS (bp)",
            y_axis_label="Mean coverage (RPKM)",
        )

        self.assertEqual(
            len(re.findall(r'class="profile-trace"[^>]*stroke="#2F78D1"', panel)),
            2,
        )
        self.assertIn('class="zero-reference"', panel)
        self.assertIn("Mean coverage (RPKM)", panel)
        self.assertEqual(panel.count('class="endpoint-label"'), 2)
        self.assertIn(">S1<", panel)
        self.assertIn(">S2<", panel)


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

    def test_nonempty_all_missing_bar_and_range_inputs_are_unavailable(self):
        """Staged identities without values must not create a fake 0–1 chart."""
        bar = visuals.render_bar_panel(
            "Mapped reads",
            [
                {"sample_id": "S1", "mapped_percent": None},
                {"sample_id": "S2", "mapped_percent": float("nan")},
            ],
            value_key="mapped_percent",
            axis_label="%",
        )
        range_panel = visuals.render_range_panel(
            "Peak width median and range",
            [{"sample_id": "S1", "minimum": None, "q25": None,
              "median": None, "q75": None, "maximum": None}],
            minimum_key="minimum",
            q25_key="q25",
            median_key="median",
            q75_key="q75",
            maximum_key="maximum",
            axis_label="Peak width (bp)",
        )

        for panel in (bar, range_panel):
            self.assertIn('role="status"', panel)
            self.assertIn("No numeric data available", panel)
            self.assertNotIn("<svg", panel)

    def test_short_labels_stay_horizontal_and_long_labels_fit_inside_viewbox(self):
        short = visuals.render_bar_panel(
            "Metric",
            [{"sample_id": "S1", "assay_target": "CTCF", "is_control": False, "value": 1}],
            value_key="value",
            axis_label="Units",
        )
        long_id = "long_sample_identifier_" + "replicate_alpha_" * 8
        long = visuals.render_bar_panel(
            "Metric",
            [{
                "sample_id": long_id,
                "assay_target": "CTCF",
                "is_control": False,
                "value": 1,
            }],
            value_key="value",
            axis_label="Units",
        )

        self.assertNotIn('class="sample-label rotated"', short)
        self.assertIn('class="sample-label rotated"', long)
        viewbox = re.search(
            r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', long
        )
        label = re.search(
            r'class="sample-label rotated" x="([0-9.]+)" y="([0-9.]+)"',
            long,
        )
        self.assertIsNotNone(viewbox)
        self.assertIsNotNone(label)
        projected_text = len(long_id) * 7.0 / math.sqrt(2.0)
        self.assertGreaterEqual(
            float(viewbox.group(1)),
            float(label.group(1)) + projected_text + 10.0,
        )
        self.assertGreaterEqual(
            float(viewbox.group(2)),
            float(label.group(2)) + projected_text + 10.0,
        )

    def test_long_range_labels_also_fit_inside_viewbox(self):
        long_id = "range_sample_" + "very_long_identifier_" * 8
        panel = visuals.render_range_panel(
            "Peak width",
            [{
                "sample_id": long_id,
                "minimum": 10,
                "q25": 20,
                "median": 30,
                "q75": 40,
                "maximum": 50,
            }],
            minimum_key="minimum",
            q25_key="q25",
            median_key="median",
            q75_key="q75",
            maximum_key="maximum",
            axis_label="Peak width (bp)",
        )

        viewbox = re.search(
            r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', panel
        )
        label = re.search(
            r'class="sample-label rotated" x="([0-9.]+)" y="([0-9.]+)"',
            panel,
        )
        self.assertIsNotNone(viewbox)
        self.assertIsNotNone(label)
        projected_text = len(long_id) * 7.0 / math.sqrt(2.0)
        self.assertGreaterEqual(
            float(viewbox.group(1)),
            float(label.group(1)) + projected_text + 10.0,
        )
        self.assertGreaterEqual(
            float(viewbox.group(2)),
            float(label.group(2)) + projected_text + 10.0,
        )

    def test_scatter_long_direct_label_has_intrinsic_scrollable_rail(self):
        long_id = "scatter_sample_" + "very_long_identifier_" * 8
        panel = visuals.render_scatter_panel(
            "Peak count vs usable fragments",
            [{
                "sample_id": long_id,
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

        svg = re.search(
            r'<svg class="panel-chart scatter-chart" width="([0-9.]+)" '
            r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"',
            panel,
        )
        label_x = re.search(
            r'class="scatter-label" x="([0-9.]+)"', panel
        )
        self.assertIsNotNone(svg)
        self.assertIsNotNone(label_x)
        self.assertEqual(float(svg.group(1)), float(svg.group(2)))
        self.assertGreaterEqual(
            float(svg.group(2)),
            float(label_x.group(1)) + len(long_id) * 7.0 + 16.0,
        )
        self.assertIn('class="panel-scroll"', panel)

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
