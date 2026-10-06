"""Focused offline tests for ILAMB discovery, coordinates and missing data."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import netCDF4
import numpy as np
import requests
from shapely.geometry import Polygon, box
import xarray as xr

from ilamb_regions import BENCHMARKS, DownloadLinks, colour_normalization, download, main, plot_models, plt, subset_field
from ilamb_regions import area_weighted_mean, cell_areas_for_field, normalize_mean_biases, plot_mean_biases
from ilamb_regions import aggregate_bias_scores, plot_aggregate_scores
from ilamb_regions import regional_lai_et_ratio, plot_lai_et_ratios
from ilamb_regions import regional_gpp_et_ratio, plot_gpp_et_ratios


class DownloadTests(unittest.TestCase):
    def test_connection_failure_retries_then_caches(self):
        session = MagicMock()
        response = MagicMock()
        response.__enter__.return_value = response
        response.iter_content.return_value = [b"netcdf data"]
        session.get.side_effect = [requests.ConnectionError("disconnected"), response]
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "model.nc"
            self.assertEqual(download(session, "https://example.org/model.nc", target), target)
            self.assertEqual(target.read_bytes(), b"netcdf data")
            self.assertEqual(session.get.call_count, 2)
            self.assertFalse(target.with_suffix(".nc.part").exists())
            download(session, "https://example.org/model.nc", target)
            self.assertEqual(session.get.call_count, 2)

    def test_connection_retries_are_bounded(self):
        session = MagicMock()
        session.get.side_effect = requests.ConnectionError("disconnected")
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "model.nc"
            with self.assertRaises(requests.ConnectionError):
                download(session, "https://example.org/model.nc", target)
            self.assertEqual(session.get.call_count, 3)
            self.assertFalse(target.exists())
            self.assertFalse(target.with_suffix(".nc.part").exists())


class RegionTests(unittest.TestCase):
    def make_dataset(self, path):
        with netCDF4.Dataset(path, "w") as dataset:
            group = dataset.createGroup("MeanState")
            group.createDimension("lat_", 3)
            group.createDimension("lon_", 3)
            latitude = group.createVariable("lat_", "f8", ("lat_",))
            longitude = group.createVariable("lon_", "f8", ("lon_",))
            latitude[:] = [1, 0, -1]
            longitude[:] = [289, 290, 291]
            field = group.createVariable("bias_map_of_lai", "f8", ("lat_", "lon_"))
            field.units = "1"
            field[:] = np.ma.array(np.arange(9).reshape(3, 3),
                                   mask=[[True, False, False],
                                         [False, False, False],
                                         [False, False, False]])
            self.assertNotIn("_FillValue", field.ncattrs())

    def test_implicit_fill_and_coordinate_normalization(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.nc"
            self.make_dataset(path)
            data = subset_field(path, "bias_map_of_lai", box(-72, -2, -68, 2))
            np.testing.assert_equal(data.lon.values, [-71, -70, -69])
            np.testing.assert_equal(data.lat.values, [-1, 0, 1])
            self.assertTrue(np.isnan(data.sel(lat=1, lon=-71).item()))
            self.assertEqual(int(np.isfinite(data.values).sum()), 8)
            self.assertEqual(data.sel(lat=-1, lon=-69).item(), 8)

    def test_polygon_hole_masks_grid_centres(self):
        geometry = Polygon([(-72, -2), (-68, -2), (-68, 2), (-72, 2)],
                           holes=[[(-70.5, -0.5), (-69.5, -0.5),
                                   (-69.5, 0.5), (-70.5, 0.5)]])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.nc"
            self.make_dataset(path)
            data = subset_field(path, "bias_map_of_lai", geometry)
            self.assertTrue(np.isnan(data.sel(lat=0, lon=-70).item()))
            self.assertEqual(int(np.isfinite(data.values).sum()), 7)

    def test_missing_field_fails_with_available_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.nc"
            self.make_dataset(path)
            with self.assertRaisesRegex(ValueError, "available.*bias_map_of_lai"):
                subset_field(path, "not_a_field", box(-72, -2, -68, 2))

    def test_region_without_data_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.nc"
            self.make_dataset(path)
            with self.assertRaisesRegex(ValueError, "no grid centres"):
                subset_field(path, "bias_map_of_lai", box(10, 10, 20, 20))

    def test_discovery_excludes_benchmark_and_images(self):
        parser = DownloadLinks()
        parser.feed('<a href="MODIS_Benchmark.nc">reference</a>'
                    '<a href="MODIS_CLM-FATES.nc">model</a>'
                    '<a href="CLM-FATES_global_bias.png">image</a>')
        self.assertEqual(parser.links, ["MODIS_CLM-FATES.nc"])

    def test_discovery_uses_selected_source(self):
        for source in ("MODIS", "GLEAMv3.3a", "FLUXCOM"):
            with self.subTest(source=source):
                parser = DownloadLinks(source)
                parser.feed(f'<a href="{source}_Benchmark.nc">reference</a>'
                            f'<a href="{source}_CLM-FATES.nc">model</a>'
                            '<a href="Other_CLM-FATES.nc">other</a>')
                self.assertEqual(parser.links, [f"{source}_CLM-FATES.nc"])


class ColourScaleTests(unittest.TestCase):
    def test_cli_lai_defaults_and_overrides(self):
        for options, expected in [
            ([], (3.0, None, None)),
            (["--limit", "2"], (2.0, None, None)),
            (["--vmin", "-2", "--vmax", "4"], (None, -2.0, 4.0)),
            (["--field", "biasscore_map_of_lai"], (None, None, None)),
            (["--benchmark", "modis-et"], (None, None, None)),
            (["--benchmark", "avh15c1-lai"], (3.0, None, None)),
        ]:
            with self.subTest(options=options), tempfile.TemporaryDirectory() as directory:
                with patch("ilamb_regions.colour_normalization", wraps=colour_normalization) as normalize:
                    with patch("ilamb_regions.discover_models", side_effect=RuntimeError("stop before downloads")):
                        with self.assertRaisesRegex(RuntimeError, "stop before downloads"):
                            main(["--offline", "--cache", directory, "--output-dir", directory] + options)
                self.assertEqual(normalize.call_args.args[2:], expected)

    def test_bias_limits_ignore_extreme_model(self):
        norm = colour_normalization(np.array([-4, 0, 27]), "bias_map_of_lai",
                                    vmin=-3, vmax=3)
        self.assertEqual((norm.vmin, norm.vmax), (-3, 3))
        self.assertEqual(norm(0), 0.5)

    def test_asymmetric_bias_limits_keep_zero_centred(self):
        norm = colour_normalization(np.array([-4, 0, 27]), "bias_map_of_lai",
                                    vmin=-2, vmax=4)
        self.assertEqual((norm.vmin, norm.vmax), (-2, 4))
        self.assertEqual(norm(0), 0.5)

    def test_score_limits_apply_to_selected_field(self):
        norm = colour_normalization(np.array([0, 0.5, 1]), "biasscore_map_of_lai",
                                    vmin=0.2, vmax=0.8)
        self.assertEqual((norm.vmin, norm.vmax), (0.2, 0.8))

    def test_existing_defaults_and_symmetric_limit(self):
        values = np.array([-4, 0, 27])
        for field, kwargs, expected in [
            ("bias_map_of_lai", {}, (-27, 27)),
            ("bias_map_of_lai", {"limit": 3}, (-3, 3)),
            ("biasscore_map_of_lai", {}, (0, 1)),
        ]:
            with self.subTest(field=field, kwargs=kwargs):
                norm = colour_normalization(values, field, **kwargs)
                self.assertEqual((norm.vmin, norm.vmax), expected)

    def test_invalid_limits_fail(self):
        for kwargs in [
            {"vmin": -3}, {"vmax": 3}, {"vmin": 3, "vmax": -3},
            {"vmin": -3, "vmax": -3}, {"vmin": float("nan"), "vmax": 3},
            {"vmin": -3, "vmax": float("inf")}, {"vmin": 1, "vmax": 3},
            {"limit": 3, "vmin": -3, "vmax": 3}, {"limit": float("nan")},
        ]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                colour_normalization(np.array([-4, 27]), "bias_map_of_lai", **kwargs)


class MeanBiasTests(unittest.TestCase):
    def test_area_weights_and_missing_values(self):
        data = xr.DataArray([[0.0, 10.0], [0.0, 10.0]], dims=("lat", "lon"),
                            coords={"lat": [-1, 1], "lon": [-71, -69]})
        areas = xr.full_like(data, 1.0)
        areas.loc[dict(lon=-69)] = 3.0
        self.assertAlmostEqual(area_weighted_mean(data, areas), 7.5)
        data.loc[dict(lat=1, lon=-69)] = np.nan
        self.assertAlmostEqual(area_weighted_mean(data, areas), 6.0)

    def test_cell_areas_use_bounds_and_normalize_coordinates(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.nc"
            with netCDF4.Dataset(path, "w") as dataset:
                group = dataset.createGroup("MeanState")
                for name, length in [("lat_", 2), ("lon_", 2), ("nb", 2)]:
                    group.createDimension(name, length)
                for name, values, bounds in [
                    ("lat_", [30, -30], [[0, 60], [-60, 0]]),
                    ("lon_", [289, 291], [[288.5, 289.5], [289.5, 292.5]]),
                ]:
                    coordinate = group.createVariable(name, "f8", (name,))
                    coordinate[:] = values
                    coordinate.bounds = name + "bounds"
                    group.createVariable(coordinate.bounds, "f8", (name, "nb"))[:] = bounds
                group.createVariable("bias_map_of_lai", "f8", ("lat_", "lon_"))
            areas = cell_areas_for_field(path, "bias_map_of_lai")
            np.testing.assert_equal(areas.lat.values, [-30, 30])
            np.testing.assert_equal(areas.lon.values, [-71, -69])
            self.assertAlmostEqual(float(areas.sel(lat=30, lon=-69) / areas.sel(lat=30, lon=-71)), 3)
            expected = 6371000.0 ** 2 * np.sin(np.deg2rad(60)) * np.deg2rad(1)
            self.assertAlmostEqual(float(areas.sel(lat=30, lon=-71)) / expected, 1)

    def test_normalization_is_signed_and_separate_per_benchmark(self):
        rows = [{"benchmark": "lai", "mean_bias": -2.0},
                {"benchmark": "lai", "mean_bias": 4.0},
                {"benchmark": "et", "mean_bias": 0.0},
                {"benchmark": "et", "mean_bias": 0.0}]
        normalized = normalize_mean_biases(rows)
        self.assertEqual([row["normalized_bias"] for row in normalized], [-0.5, 1, 0, 0])
        self.assertEqual([row["normalization_scale"] for row in normalized], [4, 4, 0, 0])
        self.assertNotIn("normalized_bias", rows[0])

    def test_mean_plot_has_three_points_without_excluded_references_and_bold_fates(self):
        models = ["CLM-FATES", "CLM6.0"]
        rows = [{"model": model, "benchmark": benchmark, "normalized_bias": 0.25}
                for model in models for benchmark in BENCHMARKS]
        with patch("ilamb_regions.plt.close"), patch("matplotlib.figure.Figure.savefig"):
            plot_mean_biases(rows, models, "Amazon basin", Path("unused"))
            figure = plt.gcf()
        try:
            axis = figure.axes[0]
            self.assertEqual(len(axis.collections), 3)
            self.assertFalse(any(reference in label for reference in ("GLEAM", "WECANN", "AVH15C1")
                                 for label in axis.get_legend_handles_labels()[1]))
            self.assertTrue(all(len(collection.get_offsets()) == 2 for collection in axis.collections))
            self.assertEqual([label.get_fontweight() for label in axis.get_xticklabels()], ["bold", "normal"])
            self.assertEqual(axis.get_ylabel(), "Normalized area-weighted mean bias")
        finally:
            plt.close(figure)


class LaiEtRatioTests(unittest.TestCase):
    def make_fields(self):
        lai = xr.DataArray([[1.0, 9.0], [1.0, 9.0]], dims=("lat", "lon"),
                           coords={"lat": [-1, 1], "lon": [-71, -69]}, attrs={"units": "1"})
        et = lai.copy(data=np.array([[1.0, 3.0], [1.0, 3.0]]))
        et.attrs = {"units": "mm d-1"}
        areas = xr.full_like(lai, 1.0)
        areas.loc[dict(lon=-69)] = 3.0
        return lai, et, areas

    def test_ratio_of_area_weighted_means_not_mean_of_ratios(self):
        lai, et, areas = self.make_fields()
        result = regional_lai_et_ratio(lai, et, areas, areas)
        self.assertEqual(result["mean_lai"], 7)
        self.assertEqual(result["mean_et_mm_day"], 2.5)
        self.assertAlmostEqual(result["lai_et_ratio_day_mm"], 2.8)
        self.assertEqual(result["valid_cells"], 4)
        self.assertEqual(result["valid_area_m2"], 8)

    def test_joint_valid_mask_is_used_for_both_means(self):
        lai, et, areas = self.make_fields()
        lai.loc[dict(lat=-1, lon=-69)] = np.nan
        et.loc[dict(lat=1, lon=-69)] = np.nan
        result = regional_lai_et_ratio(lai, et, areas, areas)
        self.assertEqual(result["mean_lai"], 1)
        self.assertEqual(result["mean_et_mm_day"], 1)
        self.assertEqual(result["valid_cells"], 2)

    def test_nonpositive_denominator_and_mismatched_grids_fail(self):
        lai, et, areas = self.make_fields()
        with self.assertRaisesRegex(ValueError, "positive regional mean ET"):
            regional_lai_et_ratio(lai, et.copy(data=np.zeros_like(et.values)), areas, areas)
        with self.assertRaises(ValueError):
            regional_lai_et_ratio(lai, et.assign_coords(lon=[-72, -70]), areas, areas)
        with self.assertRaises(ValueError):
            regional_lai_et_ratio(lai, et, areas, areas * 2)

    def test_plot_has_observation_line_models_on_y_and_bold_fates(self):
        highlighted = ["CLM-FATES", "CLM6.0", "JULES", "EDv3", "LPX-Bern", "LPJmL",
                   "JSBACH", "CLASSIC", "VISIT", "LPJ-EOSIM", "DLEM"]
        rows = [{"model": model, "lai_et_ratio_day_mm": 1.1}
            for model in highlighted + ["VISIT-UT", "E3SM"]]
        observed = {"lai_et_ratio_day_mm": 1.2}
        with patch("ilamb_regions.plt.close"), patch("matplotlib.figure.Figure.savefig"):
            plot_lai_et_ratios(rows, observed, "Amazon basin", Path("unused"))
            figure = plt.gcf()
        try:
            axis = figure.axes[0]
            self.assertEqual(axis.get_ylabel(), "Model")
            np.testing.assert_equal(axis.lines[0].get_xdata(), [1.2, 1.2])
            self.assertEqual(axis.get_yticklabels()[0].get_fontweight(), "bold")
            for label in axis.get_yticklabels():
                self.assertEqual(label.get_color(), "red" if label.get_text() in highlighted else "black")
            self.assertTrue(axis.yaxis_inverted())
        finally:
            plt.close(figure)


class GppEtRatioTests(unittest.TestCase):
    def test_gpp_ratio_uses_means_and_common_finite_cells(self):
        gpp = xr.DataArray([[2.0, 12.0], [2.0, 12.0]], dims=("lat", "lon"),
                           coords={"lat": [-1, 1], "lon": [-71, -69]}, attrs={"units": "g m-2 d-1"})
        et = gpp.copy(data=np.array([[1.0, 3.0], [1.0, 3.0]]))
        et.attrs = {"units": "mm d-1"}
        areas = xr.full_like(gpp, 1.0)
        areas.loc[dict(lon=-69)] = 3.0
        result = regional_gpp_et_ratio(gpp, et, areas, areas)
        self.assertEqual(result["mean_gpp_g_m2_day"], 9.5)
        self.assertEqual(result["mean_et_mm_day"], 2.5)
        self.assertAlmostEqual(result["gpp_et_ratio_g_m2_mm"], 3.8)
        gpp.loc[dict(lat=-1, lon=-69)] = np.nan
        et.loc[dict(lat=1, lon=-69)] = np.nan
        self.assertEqual(regional_gpp_et_ratio(gpp, et, areas, areas)["gpp_et_ratio_g_m2_mm"], 2)

    def test_gpp_ratio_plot_labels_and_red_models(self):
        rows = [{"model": "CLM-FATES", "gpp_et_ratio_g_m2_mm": 2.1},
                {"model": "VISIT-UT", "gpp_et_ratio_g_m2_mm": 2.8}]
        with patch("ilamb_regions.plt.close"), patch("matplotlib.figure.Figure.savefig"):
            plot_gpp_et_ratios(rows, {"gpp_et_ratio_g_m2_mm": 2.4}, "Amazon basin", Path("unused"))
            figure = plt.gcf()
        try:
            axis = figure.axes[0]
            np.testing.assert_equal(axis.lines[0].get_xdata(), [2.4, 2.4])
            self.assertIn("FLUXCOM GPP / MODIS ET", axis.lines[0].get_label())
            self.assertIn("mean GPP", axis.get_xlabel())
            self.assertEqual(axis.get_yticklabels()[0].get_color(), "red")
            self.assertEqual(axis.get_yticklabels()[0].get_fontweight(), "bold")
            self.assertEqual(axis.get_yticklabels()[1].get_color(), "black")
        finally:
            plt.close(figure)


class AggregateScoreTests(unittest.TestCase):
    def make_rows(self, model, values):
        return [{"model": model, "benchmark": key, "normalized_bias": value}
                for key, value in zip(("modis-et", "gleam-et", "modis-lai", "fluxcom-gpp", "avh15c1-lai", "wecann-gpp"), values)]

    def test_both_bias_signs_penalize_without_category_cancellation(self):
        rows = self.make_rows("positive", [0.5] * 6)
        rows += self.make_rows("negative", [-0.5] * 6)
        rows += self.make_rows("mixed", [0.5, -0.5, 0.5, -0.5, -0.5, 0.5])
        scores = aggregate_bias_scores(rows, ["positive", "negative", "mixed"])
        self.assertTrue(all(row["aggregate_score"] == 0.5 for row in scores))

    def test_et_references_do_not_double_weight_et(self):
        rows = self.make_rows("example", [1, 1, 0, 0, 0, 0])
        score = aggregate_bias_scores(rows, ["example"])[0]
        self.assertEqual(score["et_score"], 0)
        self.assertAlmostEqual(score["aggregate_score"], 2 / 3)

    def test_gleam_does_not_affect_scores_or_export(self):
        baseline = aggregate_bias_scores(self.make_rows("example", [0.2, 0, 0, 0, 0, 0]), ["example"])
        changed = aggregate_bias_scores(self.make_rows("example", [0.2, 1, 0, 0, 0, 0]), ["example"])
        self.assertEqual(baseline, changed)
        self.assertEqual(baseline[0]["et_score"], 0.8)
        self.assertNotIn("gleam_et_score", baseline[0])

    def test_wecann_does_not_affect_scores_or_export(self):
        baseline = aggregate_bias_scores(self.make_rows("example", [0, 0, 0, 0.2, 0, 0]), ["example"])
        changed = aggregate_bias_scores(self.make_rows("example", [0, 0, 0, 0.2, 0, 1]), ["example"])
        self.assertEqual(baseline, changed)
        self.assertEqual(baseline[0]["gpp_score"], 0.8)
        self.assertNotIn("wecann_gpp_score", baseline[0])

    def test_zero_bias_is_one_and_models_are_ranked(self):
        rows = self.make_rows("worst", [1, -1, 1, -1, -1, 1])
        rows += self.make_rows("best", [0] * 6)
        scores = aggregate_bias_scores(rows, ["worst", "best"])
        self.assertEqual([row["model"] for row in scores], ["best", "worst"])
        self.assertEqual([row["aggregate_score"] for row in scores], [1, 0])

    def test_avh15c1_does_not_affect_lai_scores_or_export(self):
        baseline = aggregate_bias_scores(self.make_rows("example", [0, 0, 0.2, 0, 0, 0]), ["example"])
        changed = aggregate_bias_scores(self.make_rows("example", [0, 0, 0.2, 0, 1, 0]), ["example"])
        self.assertEqual(baseline, changed)
        self.assertEqual(baseline[0]["lai_score"], 0.8)
        self.assertNotIn("avh15c1_lai_score", baseline[0])

    def test_missing_or_invalid_benchmark_scores_fail(self):
        with self.assertRaises(ValueError):
            aggregate_bias_scores(self.make_rows("model", [0, 0, 0]), ["model"])
        with self.assertRaises(ValueError):
            aggregate_bias_scores(self.make_rows("model", [0, 0, 0, 2, 0, 0]), ["model"])

    def test_requested_aggregate_names_are_red_and_others_black(self):
        highlighted = {"CLM-FATES", "CLM6.0", "JULES", "EDv3", "LPX-Bern", "LPJmL",
                       "JSBACH", "CLASSIC", "VISIT", "LPJ-EOSIM", "DLEM"}
        models = sorted(highlighted) + ["VISIT-UT", "ORCHIDEEv4", "E3SM"]
        rows = [row for model in models for row in self.make_rows(model, [0] * 6)]
        scores = aggregate_bias_scores(rows, models)
        with patch("ilamb_regions.plt.close"), patch("matplotlib.figure.Figure.savefig"):
            plot_aggregate_scores(scores, "Amazon basin", Path("unused"))
            figure = plt.gcf()
        try:
            labels = {label.get_text(): label for label in figure.axes[0].get_yticklabels()}
            for model in models:
                self.assertEqual(labels[model].get_color(), "red" if model in highlighted else "black")
            self.assertEqual(labels["CLM-FATES"].get_fontweight(), "bold")
        finally:
            plt.close(figure)

    def test_plot_models_are_on_y_axis_and_fates_is_bold(self):
        scores = aggregate_bias_scores(self.make_rows("CLM-FATES", [0.2, -0.2, 0.3, -0.4, 0.2, -0.3]), ["CLM-FATES"])
        with patch("ilamb_regions.plt.close"), patch("matplotlib.figure.Figure.savefig"):
            plot_aggregate_scores(scores, "Amazon basin", Path("unused"))
            figure = plt.gcf()
        try:
            axis = figure.axes[0]
            self.assertEqual(axis.get_ylabel(), "Model")
            self.assertEqual(axis.get_yticklabels()[0].get_text(), "CLM-FATES")
            self.assertEqual(axis.get_yticklabels()[0].get_fontweight(), "bold")
            self.assertTrue(axis.yaxis_inverted())
            self.assertEqual(len(axis.collections), 5)
            for collection in axis.collections[1:]:
                np.testing.assert_equal(collection.get_offsets()[:, 1], [0])
            self.assertFalse(any(reference in label for reference in ("GLEAM", "WECANN", "AVH15C1")
                                 for label in axis.get_legend_handles_labels()[1]))
        finally:
            plt.close(figure)


class PlotStyleTests(unittest.TestCase):
    def test_et_and_gpp_labels_use_reference_and_native_units(self):
        for source, variable, unit in [("MODIS", "ET", "mm d-1"),
                                       ("GLEAMv3.3a", "ET", "mm d-1"),
                                       ("FLUXCOM", "GPP", "g m-2 d-1")]:
            with self.subTest(source=source, variable=variable):
                field = f"bias_map_of_{variable.lower()}"
                data = xr.DataArray([[0.1, 0.3], [0.5, 0.7]], dims=("lat", "lon"),
                                    coords={"lat": [-1, 1], "lon": [-71, -69]},
                                    attrs={"units": unit})
                with patch("ilamb_regions.plt.close"), patch("matplotlib.figure.Figure.savefig"):
                    plot_models([("CLM-FATES", data)], box(-72, -2, -68, 2), field,
                                "Amazon basin", Path("unused"), source=source, variable=variable)
                    figure = plt.gcf()
                try:
                    self.assertIn(f"{source} {variable}", figure._suptitle.get_text())
                    self.assertEqual(figure.axes[-1].get_xlabel(),
                                     f"{variable} bias: model minus {source} [{unit}]")
                    self.assertEqual(figure.axes[0].title.get_fontweight(), "bold")
                finally:
                    plt.close(figure)

    def test_only_fates_title_is_bold_for_bias_and_score(self):
        data = xr.DataArray([[0.1, 0.3], [0.5, 0.7]], dims=("lat", "lon"),
                            coords={"lat": [-1, 1], "lon": [-71, -69]},
                            attrs={"units": "1"})
        results = [(model, data) for model in ("CLM-FATES", "CLM6.0", "E3SM")]
        for field in ("bias_map_of_lai", "biasscore_map_of_lai"):
            with self.subTest(field=field):
                with patch("ilamb_regions.plt.close"), patch("matplotlib.figure.Figure.savefig"):
                    plot_models(results, box(-72, -2, -68, 2), field, "Amazon basin",
                                Path("unused"))
                    figure = plt.gcf()
                try:
                    titles = {axis.get_title(): axis.title.get_fontweight()
                              for axis in figure.axes if axis.get_title()}
                    self.assertEqual(titles, {"CLM-FATES": "bold", "CLM6.0": "normal",
                                              "E3SM": "normal"})
                finally:
                    plt.close(figure)


if __name__ == "__main__":
    unittest.main()