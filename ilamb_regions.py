"""Download ILAMB spatial results and compare models within a watershed."""

import argparse
import csv
import getpass
import io
import json
import math
import os
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urljoin, urlparse
import zipfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import netCDF4
import numpy as np
import requests
import shapefile
from shapely.geometry import mapping, shape
from shapely.ops import unary_union
from shapely import contains_xy
import xarray as xr


DEFAULT_URL = (
    "https://gws-access.jasmin.ac.uk/public/landsurf_rdg/pmcguire/"
    "ILAMB_output_TRENDYv15/GCP_270926a4/"
)
BENCHMARK_PATH = "EcosystemandCarbonCycle/LeafAreaIndex/MODIS/"
BENCHMARKS = {
    "modis-lai": {"path": BENCHMARK_PATH, "source": "MODIS", "variable": "LAI",
                  "field": "bias_map_of_lai", "slug": "modis_lai"},
    "modis-et": {"path": "HydrologyCycle/Evapotranspiration/MODIS/", "source": "MODIS",
                 "variable": "ET", "field": "bias_map_of_et", "slug": "modis_et"},
    "gleam-et": {"path": "HydrologyCycle/Evapotranspiration/GLEAMv3.3a/",
                 "source": "GLEAMv3.3a", "variable": "ET", "field": "bias_map_of_et",
                 "slug": "gleam_et"},
    "fluxcom-gpp": {"path": "EcosystemandCarbonCycle/GrossPrimaryProductivity/FLUXCOM/",
                    "source": "FLUXCOM", "variable": "GPP", "field": "bias_map_of_gpp",
                    "slug": "fluxcom_gpp"},
    "wecann-gpp": {"path": "EcosystemandCarbonCycle/GrossPrimaryProductivity/WECANN/",
                   "source": "WECANN", "variable": "GPP", "field": "bias_map_of_gpp",
                   "slug": "wecann_gpp"},
    "avh15c1-lai": {"path": "EcosystemandCarbonCycle/LeafAreaIndex/AVH15C1/",
                    "source": "AVH15C1", "variable": "LAI", "field": "bias_map_of_lai",
                    "slug": "avh15c1_lai"},
}
COMPARISON_BENCHMARKS = tuple(key for key in BENCHMARKS
                             if key not in ("gleam-et", "wecann-gpp", "avh15c1-lai"))
BASIN_URL = (
    "https://data.hydrosheds.org/file/hydrobasins/standard/"
    "hybas_sa_lev04_v1c.zip"
)
AMAZON_MAIN_BAS = 6040007000
AGGREGATE_RED_MODELS = {
    "CLM-FATES", "CLM6.0", "JULES", "EDv3", "LPX-Bern", "LPJmL",
    "JSBACH", "CLASSIC", "VISIT", "LPJ-EOSIM", "DLEM",
}
ROOT = Path(__file__).resolve().parent
ET_RATIO_SPECS = {
    "LAI": {"benchmark": "modis-lai", "field": "timeint_of_lai", "native_units": "1",
        "mean_key": "mean_lai", "ratio_key": "lai_et_ratio_day_mm",
        "units": "day/mm", "reference_label": "MODIS"},
    "GPP": {"benchmark": "fluxcom-gpp", "field": "timeint_of_gpp", "native_units": "g m-2 d-1",
        "mean_key": "mean_gpp_g_m2_day", "ratio_key": "gpp_et_ratio_g_m2_mm",
        "units": "g/m2/mm", "reference_label": "FLUXCOM GPP / MODIS ET"},
}


class DownloadLinks(HTMLParser):
    def __init__(self, source="MODIS"):
        super().__init__()
        self.source = source
        self.links = []

    def handle_starttag(self, tag, attrs):
        href = dict(attrs).get("href", "")
        if tag == "a" and href.startswith(self.source + "_") and href.endswith(".nc"):
            if href != self.source + "_Benchmark.nc":
                self.links.append(href)


def download(session, url, target, offline=False):
    if target.exists():
        return target
    if offline:
        raise FileNotFoundError(f"Not cached: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    for attempt in range(3):
        try:
            with session.get(url, timeout=(30, 180), stream=True) as response:
                response.raise_for_status()
                with partial.open("wb") as handle:
                    for chunk in response.iter_content(1024 * 1024):
                        handle.write(chunk)
            partial.replace(target)
            return target
        except (requests.ConnectionError, requests.Timeout):
            if attempt == 2:
                raise
            print(f"Retrying {target.name} after connection failure ({attempt + 1}/2)", flush=True)
        finally:
            partial.unlink(missing_ok=True)


def discover_models(session, benchmark_url, cache, offline=False, source="MODIS"):
    page = source + ".html"
    index = download(session, urljoin(benchmark_url, page), cache / page, offline)
    parser = DownloadLinks(source)
    parser.feed(index.read_text())
    links = list(dict.fromkeys(parser.links))
    if not links:
        raise ValueError(f"No model NetCDF download links found in {page}")
    return links


def amazon_geometry(cache, offline=False):
    archive_path = cache / "hybas_sa_lev04_v1c.zip"
    with requests.Session() as public_session:
        download(public_session, BASIN_URL, archive_path, offline)
    stem = "hybas_sa_lev04_v1c"
    with zipfile.ZipFile(archive_path) as archive:
        reader = shapefile.Reader(
            shp=io.BytesIO(archive.read(stem + ".shp")),
            dbf=io.BytesIO(archive.read(stem + ".dbf")),
        )
        members = [shape(item.shape.__geo_interface__)
                   for item in reader.iterShapeRecords()
                   if item.record["MAIN_BAS"] == AMAZON_MAIN_BAS]
    if not members:
        raise ValueError("Amazon MAIN_BAS not found in HydroBASINS archive")
    return unary_union(members)


def read_geometry(path):
    content = json.loads(path.read_text())
    if content["type"] == "FeatureCollection":
        geometry = unary_union([shape(item["geometry"]) for item in content["features"]])
    elif content["type"] == "Feature":
        geometry = shape(content["geometry"])
    else:
        geometry = shape(content)
    if geometry.is_empty or not geometry.is_valid:
        raise ValueError("Region must be a valid, nonempty WGS84 polygon")
    if geometry.geom_type not in ("Polygon", "MultiPolygon"):
        raise ValueError("Region must contain polygons")
    return geometry


def subset_field(path, field, geometry):
    with xr.open_dataset(path, group="MeanState") as dataset:
        if field not in dataset:
            available = [name for name in dataset.data_vars if "map" in name]
            raise ValueError(f"{path.name}: missing {field}; available: {available}")
        data = dataset[field].load().squeeze(drop=True)
    with netCDF4.Dataset(path) as dataset:
        native = dataset.groups["MeanState"].variables[field][:]
        data = data.copy(data=np.ma.filled(native, np.nan).squeeze())
    latitude = [dim for dim in data.dims if dim.startswith("lat")]
    longitude = [dim for dim in data.dims if dim.startswith("lon")]
    if data.ndim != 2 or len(latitude) != 1 or len(longitude) != 1:
        raise ValueError(f"Expected a lat/lon map, got {data.dims}")
    data = data.rename({latitude[0]: "lat", longitude[0]: "lon"})
    data = data.assign_coords(lon=(data.lon + 180) % 360 - 180)
    data = data.sortby("lat").sortby("lon").transpose("lat", "lon")
    west, south, east, north = geometry.bounds
    data = data.sel(lat=slice(south, north), lon=slice(west, east))
    if not data.size:
        raise ValueError(f"{path.name}: no grid centres in region bounds")
    longitudes, latitudes = np.meshgrid(data.lon.values, data.lat.values)
    mask = contains_xy(geometry, longitudes, latitudes)
    data = data.where(mask)
    if not np.isfinite(data.values).any():
        raise ValueError(f"{path.name}: no valid {field} values inside region")
    data.attrs.update(region_mask="Grid centre inside watershed polygon")
    return data


def colour_normalization(values, field, limit=None, vmin=None, vmax=None):
    signed = field.startswith(("bias_", "shift_"))
    if (vmin is None) != (vmax is None):
        raise ValueError("Supply --vmin and --vmax together")
    if limit is not None and vmin is not None:
        raise ValueError("Use either --limit or --vmin/--vmax, not both")
    if limit is not None and (not np.isfinite(limit) or limit <= 0):
        raise ValueError("--limit must be finite and positive")
    if vmin is not None:
        if not np.isfinite([vmin, vmax]).all() or vmin >= vmax:
            raise ValueError("Colour limits must be finite, with --vmin < --vmax")
        if signed:
            if not vmin < 0 < vmax:
                raise ValueError("Bias and shift colour limits must straddle zero")
            return TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        return Normalize(vmin=vmin, vmax=vmax)
    if signed:
        maximum = limit if limit is not None else float(np.max(np.abs(values)))
        if maximum <= 0:
            maximum = 1.0
        return TwoSlopeNorm(vmin=-maximum, vcenter=0, vmax=maximum)
    if "score" in field:
        return Normalize(vmin=0, vmax=1)
    return Normalize(vmin=float(values.min()), vmax=float(values.max()))


def plot_models(results, geometry, field, region_name, output, columns=5, limit=None,
                vmin=None, vmax=None, source="MODIS", variable="LAI"):
    units = {data.attrs.get("units", "1") for _, data in results}
    if len(units) != 1:
        raise ValueError(f"Inconsistent units across models: {units}")
    unit = units.pop()
    values = np.concatenate([data.values[np.isfinite(data.values)] for _, data in results])
    signed = field.startswith(("bias_", "shift_"))
    norm = colour_normalization(values, field, limit, vmin, vmax)
    cmap = "RdBu_r" if signed else "viridis"
    below = bool(np.any(values < norm.vmin))
    above = bool(np.any(values > norm.vmax))
    extension = "both" if below and above else "min" if below else "max" if above else "neither"
    columns = min(columns, len(results))
    rows = math.ceil(len(results) / columns)
    figure, axes = plt.subplots(rows, columns, figsize=(3.5 * columns, 3.1 * rows),
                                squeeze=False)
    west, south, east, north = geometry.bounds
    polygons = list(geometry.geoms) if geometry.geom_type == "MultiPolygon" else [geometry]
    for axis, (model, data) in zip(axes.flat, results):
        image = axis.pcolormesh(data.lon, data.lat, data.values, cmap=cmap,
                               norm=norm, shading="auto", rasterized=True)
        for polygon in polygons:
            longitude, latitude = polygon.exterior.xy
            axis.plot(longitude, latitude, color="0.25", linewidth=0.5)
            for interior in polygon.interiors:
                longitude, latitude = interior.xy
                axis.plot(longitude, latitude, color="0.25", linewidth=0.3)
        axis.set(xlim=(west - 0.7, east + 0.7), ylim=(south - 0.7, north + 0.7),
                 title=model, facecolor="#f0f0f0")
        axis.title.set_fontweight("bold" if model == "CLM-FATES" else "normal")
        axis.set_aspect(1 / np.cos(np.deg2rad((south + north) / 2)))
        axis.tick_params(labelsize=8)
        axis.grid(alpha=0.15, linewidth=0.5)
        axis.set_xlabel("Longitude", fontsize=8)
        axis.set_ylabel("Latitude", fontsize=8)
    for axis in list(axes.flat)[len(results):]:
        axis.set_visible(False)
    metric = field.split("_map_of_", 1)[0].replace("_", " ")
    display_unit = "m2/m2" if variable == "LAI" and unit == "1" and "score" not in field else unit
    if field.startswith("bias_map_"):
        label = f"{variable} bias: model minus {source} [{display_unit}]"
    else:
        label = f"{variable} {metric} [{display_unit}]"
    figure.suptitle(f"{region_name} | {source} {variable} | {len(results)} TRENDYv15 models",
                   fontsize=16)
    figure.subplots_adjust(left=0.055, right=0.97, top=0.93, bottom=0.12,
                           hspace=0.38, wspace=0.3)
    colorbar_axis = figure.add_axes((0.25, 0.045, 0.5, 0.015))
    colorbar = figure.colorbar(image, cax=colorbar_axis, orientation="horizontal",
                              extend=extension)
    colorbar.set_label(label)
    for suffix in (".png", ".pdf"):
        figure.savefig(output.with_suffix(suffix), dpi=180, facecolor="white")
    plt.close(figure)
    return float(norm.vmin), float(norm.vmax)


def cell_areas_for_field(path, field):
    with netCDF4.Dataset(path) as dataset:
        group = dataset.groups["MeanState"]
        dimensions = group.variables[field].dimensions
        latitude_name = next(name for name in dimensions if name.startswith("lat"))
        longitude_name = next(name for name in dimensions if name.startswith("lon"))
        latitude = group.variables[latitude_name]
        longitude = group.variables[longitude_name]
        latitude_bounds = np.asarray(group.variables[latitude.bounds][:])
        longitude_bounds = np.asarray(group.variables[longitude.bounds][:])
        latitude_values = np.asarray(latitude[:])
        longitude_values = np.asarray(longitude[:])
    latitude_factor = np.abs(np.diff(np.sin(np.deg2rad(latitude_bounds)), axis=1)).ravel()
    longitude_width = np.abs(np.diff(longitude_bounds, axis=1)).ravel()
    longitude_width = np.minimum(longitude_width, 360 - longitude_width)
    areas = 6371000.0 ** 2 * latitude_factor[:, None] * np.deg2rad(longitude_width)[None, :]
    data = xr.DataArray(areas, dims=("lat", "lon"),
                        coords={"lat": latitude_values, "lon": (longitude_values + 180) % 360 - 180})
    return data.sortby("lat").sortby("lon")


def area_weighted_mean(data, areas):
    weights = areas.reindex_like(data)
    if not np.isfinite(weights.values).all() or np.any(weights.values <= 0):
        raise ValueError("Missing or invalid cell areas for the bias grid")
    finite = data.where(np.isfinite(data))
    if not np.isfinite(finite.values).any():
        raise ValueError("No valid bias values for area averaging")
    return float(finite.weighted(weights).mean().item())


def regional_et_ratio(numerator, et, numerator_areas, et_areas, variable):
    spec = ET_RATIO_SPECS[variable]
    if numerator.attrs.get("units") != spec["native_units"] or et.attrs.get("units") != "mm d-1":
        raise ValueError(f"{variable}/ET requires {variable} in {spec['native_units']} and ET in mm/day")
    numerator, et = xr.align(numerator, et, join="exact")
    areas = numerator_areas.reindex_like(numerator)
    other_areas = et_areas.reindex_like(et)
    if not np.allclose(areas.values, other_areas.values, rtol=1e-10, atol=0):
        raise ValueError(f"{variable} and ET cell bounds must describe the same grid")
    common = np.isfinite(numerator) & np.isfinite(et)
    mean_numerator = area_weighted_mean(numerator.where(common), areas)
    mean_et = area_weighted_mean(et.where(common), areas)
    if mean_et <= 0 or mean_numerator < 0:
        raise ValueError(f"{variable}/ET requires positive regional mean ET and nonnegative mean {variable}")
    return {spec["mean_key"]: mean_numerator, "mean_et_mm_day": mean_et,
            spec["ratio_key"]: mean_numerator / mean_et,
            "valid_cells": int(common.sum().item()),
            "valid_area_m2": float(areas.where(common).sum().item())}


def regional_lai_et_ratio(lai, et, lai_areas, et_areas):
    return regional_et_ratio(lai, et, lai_areas, et_areas, "LAI")


def regional_gpp_et_ratio(gpp, et, gpp_areas, et_areas):
    return regional_et_ratio(gpp, et, gpp_areas, et_areas, "GPP")


def plot_et_ratios(rows, observed, region_name, output, variable):
    spec = ET_RATIO_SPECS[variable]
    figure, axis = plt.subplots(figsize=(11, max(5, len(rows) * 0.38 + 2)))
    positions = np.arange(len(rows))
    ratios = [row[spec["ratio_key"]] for row in rows]
    reference = observed[spec["ratio_key"]]
    axis.hlines(positions, np.minimum(ratios, reference), np.maximum(ratios, reference),
                color="0.85", linewidth=1.1, zorder=1)
    axis.scatter(ratios, positions, color="#087e8b", s=50, label="Models", zorder=3)
    axis.axvline(reference, color="#202020", linestyle="--", linewidth=1.5,
                 label=f"{spec['reference_label']} observed {variable}/ET = {reference:.3f} {spec['units']}")
    axis.set(xlim=(0, max(max(ratios), reference) * 1.08),
             ylim=(-0.6, len(rows) - 0.4), ylabel="Model",
             xlabel=f"Ratio of regional mean {variable} to mean ET ({spec['units']})")
    axis.set_yticks(positions, [row["model"] for row in rows])
    axis.invert_yaxis()
    for label in axis.get_yticklabels():
        label.set_fontweight("bold" if label.get_text() == "CLM-FATES" else "normal")
        label.set_color("red" if label.get_text() in AGGREGATE_RED_MODELS else "black")
    axis.grid(axis="x", color="0.9", linewidth=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle(f"{region_name} | {variable} / evapotranspiration", fontsize=16, y=0.98)
    figure.legend(*axis.get_legend_handles_labels(), loc="upper center",
                  bbox_to_anchor=(0.57, 0.945), ncol=2, frameon=False)
    figure.subplots_adjust(left=0.2, right=0.98, bottom=0.08, top=0.89)
    for suffix in (".png", ".pdf"):
        figure.savefig(output.with_suffix(suffix), dpi=180, facecolor="white")
    plt.close(figure)


def plot_lai_et_ratios(rows, observed, region_name, output):
    plot_et_ratios(rows, observed, region_name, output, "LAI")


def plot_gpp_et_ratios(rows, observed, region_name, output):
    plot_et_ratios(rows, observed, region_name, output, "GPP")


def summarize_et_ratios(cache, output_dir, region_name, variable, models=None, offline=False):
    spec = ET_RATIO_SPECS[variable]
    numerator_key = spec["benchmark"]
    slug = "".join(character if character.isalnum() else "_"
                   for character in region_name.lower()).strip("_")
    manifests = {}
    directories = {}
    for key in (numerator_key, "modis-et"):
        benchmark = BENCHMARKS[key]
        stem = output_dir / f"{benchmark['slug']}_{benchmark['field']}_{slug}"
        manifests[key] = json.loads(stem.with_suffix(".json").read_text())
        directory = cache / urlparse(manifests[key]["ilamb_url"]).path.strip("/").replace("/", "_")
        directories[key] = directory if key == "modis-lai" else directory / key
    numerator_manifest = manifests[numerator_key]
    et_manifest = manifests["modis-et"]
    if (numerator_manifest["ilamb_url"] != et_manifest["ilamb_url"]
            or numerator_manifest["region"] != region_name or et_manifest["region"] != region_name):
        raise ValueError(f"{variable} and ET inputs must describe the same run and region")
    models = models or numerator_manifest["models"]
    if set(models) - (set(numerator_manifest["models"]) & set(et_manifest["models"])):
        raise ValueError(f"Selected models must be available for both {variable} and ET")
    geometry = read_geometry(output_dir / numerator_manifest["boundary"])
    references = {key: directory / f"{BENCHMARKS[key]['source']}_Benchmark.nc"
                  for key, directory in directories.items()}
    with requests.Session() as session:
        if not offline and any(not path.exists() for path in references.values()):
            username = os.environ.get("ILAMB_USERNAME") or input("ILAMB username: ")
            password = os.environ.get("ILAMB_PASSWORD") or getpass.getpass("ILAMB password: ")
            session.auth = (username, password)
        for key, path in references.items():
            url = urljoin(numerator_manifest["ilamb_url"], BENCHMARKS[key]["path"] + path.name)
            download(session, url, path, offline)

    def read_pair(numerator_path, et_path):
        numerator = subset_field(numerator_path, spec["field"], geometry)
        et = subset_field(et_path, "timeint_of_et", geometry)
        return regional_et_ratio(numerator, et, cell_areas_for_field(numerator_path, spec["field"]),
                                 cell_areas_for_field(et_path, "timeint_of_et"), variable)

    observed = {"model": spec["reference_label"] + " observations", "kind": "observation",
                **read_pair(references[numerator_key], references["modis-et"])}
    rows = [{"model": model, "kind": "model",
             **read_pair(directories[numerator_key] / f"{BENCHMARKS[numerator_key]['source']}_{model}.nc",
                         directories["modis-et"] / f"MODIS_{model}.nc")} for model in models]
    if variable == "GPP":
        rows.sort(key=lambda row: (abs(row[spec["ratio_key"]]), row["model"].casefold()))
        ordering = "Ascending absolute GPP/ET ratio"
    else:
        rows.sort(key=lambda row: (abs(row[spec["ratio_key"]] - observed[spec["ratio_key"]]),
                                  row["model"].casefold()))
        ordering = "Closest to the observed ratio first"
    output = output_dir / f"{variable.lower()}_et_ratio_{slug}"
    plot_et_ratios(rows, observed, region_name, output, variable)
    with output.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows + [observed])
    metadata = {"region": region_name, "ilamb_url": numerator_manifest["ilamb_url"],
                "models": models, "observed": observed, "units": spec["units"],
                "fields": ["MeanState/" + spec["field"], "MeanState/timeint_of_et"],
                "method": "Ratio of area-weighted regional means, not mean of cell ratios",
                "mask": f"Common finite {variable} and ET cells within the basin, separately for each model or observation pair",
                "weighting": f"Spherical area from native cell bounds; {variable} and ET grids must match within each pair",
                "ordering": ordering,
                "reference_files": {key: str(path) for key, path in references.items()},
                "caveat": "Uses ILAMB-supplied period means without temporal harmonization. Native grids and valid footprints differ between models and observations; ratio agreement does not ensure either field is unbiased."}
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2))
    print(f"Saved {output.with_suffix('.png')} and PDF; observed {variable}/ET = {observed[spec['ratio_key']]:.4f} {spec['units']}", flush=True)
    return rows, observed


def summarize_lai_et_ratios(cache, output_dir, region_name, models=None, offline=False):
    return summarize_et_ratios(cache, output_dir, region_name, "LAI", models, offline)


def summarize_gpp_et_ratios(cache, output_dir, region_name, models=None, offline=False):
    return summarize_et_ratios(cache, output_dir, region_name, "GPP", models, offline)


def normalize_mean_biases(rows):
    scales = {}
    for row in rows:
        value = row["mean_bias"]
        if not np.isfinite(value):
            raise ValueError("Mean biases must be finite")
        scales[row["benchmark"]] = max(scales.get(row["benchmark"], 0.0), abs(value))
    normalized = []
    for row in rows:
        scale = scales[row["benchmark"]]
        normalized.append(dict(row, normalization_scale=scale,
                               normalized_bias=row["mean_bias"] / scale if scale else 0.0))
    return normalized


def plot_mean_biases(rows, models, region_name, output):
    styles = [
        ("modis-et", "ET | MODIS", "#087e8b", "o"),
        ("modis-lai", "LAI | MODIS", "#c7384f", "^"),
        ("fluxcom-gpp", "GPP | FLUXCOM", "#d39117", "s"),
    ]
    figure, axis = plt.subplots(figsize=(16, 6))
    positions = np.arange(len(models))
    for offset, (benchmark, label, color, marker) in zip(np.linspace(-0.3, 0.3, len(styles)), styles):
        values = {row["model"]: row["normalized_bias"] for row in rows if row["benchmark"] == benchmark}
        axis.scatter(positions + offset, [values[model] for model in models], label=label,
                     color=color, marker=marker, s=55, edgecolors="white", linewidths=0.6, zorder=3)
    axis.axhline(0, color="0.25", linewidth=1)
    axis.set(ylim=(-1.08, 1.08), xlim=(-0.7, len(models) - 0.3),
             ylabel="Normalized area-weighted mean bias", xlabel="Model")
    axis.set_xticks(positions, models, rotation=55, ha="right", rotation_mode="anchor")
    for label in axis.get_xticklabels():
        label.set_fontweight("bold" if label.get_text() == "CLM-FATES" else "normal")
    axis.set_yticks(np.linspace(-1, 1, 9))
    axis.grid(axis="y", color="0.88", linewidth=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle(f"{region_name} | normalized mean bias | {len(models)} TRENDYv15 models", fontsize=16)
    figure.legend(*axis.get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.52, 0.93),
                  ncol=3, frameon=False)
    figure.subplots_adjust(left=0.07, right=0.99, bottom=0.27, top=0.82)
    for suffix in (".png", ".pdf"):
        figure.savefig(output.with_suffix(suffix), dpi=180, facecolor="white")
    plt.close(figure)


def aggregate_bias_scores(rows, models):
    scores = []
    for model in models:
        model_rows = [row for row in rows if row["model"] == model
                      and row["benchmark"] not in ("gleam-et", "wecann-gpp", "avh15c1-lai")]
        values = {row["benchmark"]: float(row["normalized_bias"]) for row in model_rows}
        if len(model_rows) != len(COMPARISON_BENCHMARKS) or set(values) != set(COMPARISON_BENCHMARKS):
            raise ValueError(f"Expected all {len(COMPARISON_BENCHMARKS)} comparison biases exactly once for {model}; rerun --summary")
        if any(not np.isfinite(value) or abs(value) > 1 for value in values.values()):
            raise ValueError("Normalized biases must be finite and between -1 and +1")
        category = {key: 1 - abs(value) for key, value in values.items()}
        variable_scores = {
            variable.lower() + "_score": float(np.mean([
                category[key] for key in COMPARISON_BENCHMARKS if BENCHMARKS[key]["variable"] == variable
            ])) for variable in ("ET", "LAI", "GPP")
        }
        scores.append({"model": model, **variable_scores,
                       "aggregate_score": float(np.mean(list(variable_scores.values()))),
                       **{key.replace("-", "_") + "_score": category[key] for key in COMPARISON_BENCHMARKS}})
    return sorted(scores, key=lambda row: (-row["aggregate_score"], row["model"].casefold()))


def plot_aggregate_scores(scores, region_name, output):
    figure, axis = plt.subplots(figsize=(11, max(5, len(scores) * 0.38 + 2)))
    positions = np.arange(len(scores))
    styles = [
        ("modis_et_score", "ET | MODIS", "#087e8b", "o"),
        ("modis_lai_score", "LAI | MODIS", "#c7384f", "^"),
        ("fluxcom_gpp_score", "GPP | FLUXCOM", "#d39117", "s"),
        ("aggregate_score", "Aggregate", "#202020", "D"),
    ]
    axis.hlines(positions, 0, [row["aggregate_score"] for row in scores],
                color="0.87", linewidth=1.2, zorder=1)
    for key, label, color, marker in styles:
        aggregate = key == "aggregate_score"
        axis.scatter([row[key] for row in scores], positions, label=label,
                     facecolors="none" if aggregate else color, marker=marker,
                     s=90 if aggregate else 40, edgecolors=color if aggregate else "white",
                     linewidths=1.3 if aggregate else 0.5, zorder=3)
    axis.set(xlim=(-0.025, 1.025), ylim=(-0.6, len(scores) - 0.4),
             xlabel="Mean-bias score (higher is better)", ylabel="Model")
    axis.set_yticks(positions, [row["model"] for row in scores])
    axis.invert_yaxis()
    for label in axis.get_yticklabels():
        label.set_fontweight("bold" if label.get_text() == "CLM-FATES" else "normal")
        label.set_color("red" if label.get_text() in AGGREGATE_RED_MODELS else "black")
    axis.set_xticks(np.linspace(0, 1, 6))
    axis.grid(axis="x", color="0.9", linewidth=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle(f"{region_name} | aggregate mean-bias scores", fontsize=16, y=0.98)
    figure.legend(*axis.get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.57, 0.95),
                  ncol=3, frameon=False)
    figure.subplots_adjust(left=0.2, right=0.98, bottom=0.08, top=0.9)
    for suffix in (".png", ".pdf"):
        figure.savefig(output.with_suffix(suffix), dpi=180, facecolor="white")
    plt.close(figure)


def summarize_aggregate_scores(output_dir, region_name):
    slug = "".join(character if character.isalnum() else "_"
                   for character in region_name.lower()).strip("_")
    source = output_dir / f"normalized_mean_bias_{slug}"
    metadata = json.loads(source.with_suffix(".json").read_text())
    with source.with_suffix(".csv").open() as handle:
        rows = list(csv.DictReader(handle))
    scores = aggregate_bias_scores(rows, metadata["models"])
    output = output_dir / f"aggregate_mean_bias_scores_{slug}"
    plot_aggregate_scores(scores, region_name, output)
    with output.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scores[0]))
        writer.writeheader()
        writer.writerows(scores)
    metadata.update(source_summary=source.with_suffix(".csv").name,
                    category_score="1 - abs(normalized regional mean bias)",
                    variable_references={variable: [BENCHMARKS[key]["source"] for key in COMPARISON_BENCHMARKS
                                                    if BENCHMARKS[key]["variable"] == variable]
                                         for variable in ("ET", "LAI", "GPP")},
                    et_score="MODIS ET score; GLEAM excluded",
                    lai_score="MODIS LAI score; AVH15C1 excluded",
                    gpp_score="FLUXCOM GPP score; WECANN excluded",
                    aggregate_score="Equal-weight mean of ET, LAI and GPP scores",
                    ordering="Highest aggregate score first",
                    caveat="Custom relative mean-bias score, not ILAMB skill. Spatial errors can still cancel within a regional mean.")
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2))
    print(f"Saved {output.with_suffix('.png')} and PDF; {len(scores)} model scores in CSV", flush=True)
    return scores


def summarize_biases(cache, output_dir, region_name, models=None):
    slug = "".join(character if character.isalnum() else "_"
                   for character in region_name.lower()).strip("_")
    rows = []
    inputs = []
    reference_run = None
    for key in COMPARISON_BENCHMARKS:
        benchmark = BENCHMARKS[key]
        stem = output_dir / f"{benchmark['slug']}_{benchmark['field']}_{slug}"
        manifest_path = stem.with_suffix(".json")
        manifest = json.loads(manifest_path.read_text())
        if reference_run is None:
            reference_run = manifest["ilamb_url"]
            models = models or manifest["models"]
        if manifest["ilamb_url"] != reference_run or manifest["region"] != region_name:
            raise ValueError("Summary inputs must describe the same ILAMB run and region")
        if set(models) - set(manifest["models"]):
            raise ValueError(f"Missing selected models for {key}")
        inputs.append(str(manifest_path))
        raw_cache = cache / urlparse(manifest["ilamb_url"]).path.strip("/").replace("/", "_")
        if key != "modis-lai":
            raw_cache = raw_cache / key
        for model in models:
            with xr.open_dataset(stem.with_name(stem.name + "_subsets") / f"{model}.nc") as dataset:
                data = dataset[benchmark["field"]].load()
            areas = cell_areas_for_field(raw_cache / f"{benchmark['source']}_{model}.nc", benchmark["field"])
            mean = area_weighted_mean(data, areas)
            valid_area = float(areas.reindex_like(data).where(np.isfinite(data)).sum().item())
            rows.append({"model": model, "benchmark": key, "variable": benchmark["variable"],
                         "reference": benchmark["source"], "units": data.attrs.get("units", "1"),
                         "mean_bias": mean, "valid_cells": int(np.isfinite(data.values).sum()),
                         "valid_area_m2": valid_area})
    rows = normalize_mean_biases(rows)
    output = output_dir / f"normalized_mean_bias_{slug}"
    plot_mean_biases(rows, models, region_name, output)
    with output.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {"region": region_name, "ilamb_url": reference_run, "models": models,
                "inputs": inputs, "normalization": "mean bias / largest absolute model mean bias per benchmark",
                "zero_scale": "All zero biases map to zero", "weighting": "Spherical cell area from native bounds",
                "mask": "Saved basin grid-centre mask; average over valid cells only",
                "caveat": "Signed mean errors can cancel; this is not an ILAMB skill score. Scales depend on the plotted models."}
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2))
    print(f"Saved {output.with_suffix('.png')} and PDF; {len(rows)} means in CSV", flush=True)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--cache", type=Path, default=ROOT / "cache")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--models", nargs="+", help="Default: all discovered models")
    parser.add_argument("--benchmark", choices=BENCHMARKS, default="modis-lai",
                        help="Reference dataset and variable; default: modis-lai")
    parser.add_argument("--field", help="MeanState variable; default: selected benchmark's bias map")
    parser.add_argument("--region-geojson", type=Path, help="Optional WGS84 watershed polygon")
    parser.add_argument("--region-name", default="Amazon basin")
    parser.add_argument("--columns", type=int, default=5)
    parser.add_argument("--limit", type=float, help="Symmetric bias colour limit; LAI bias defaults to 3")
    parser.add_argument("--vmin", type=float, help="Minimum colour value for the selected field")
    parser.add_argument("--vmax", type=float, help="Maximum colour value for the selected field")
    parser.add_argument("--summary", action="store_true",
                        help="Compare normalized mean biases from all saved benchmarks; no downloads")
    parser.add_argument("--aggregate", action="store_true",
                        help="Plot category and aggregate scores from the saved mean-bias summary")
    parser.add_argument("--lai-et-ratio", action="store_true",
                        help="Plot model regional mean LAI/ET ratios with the MODIS observed reference")
    parser.add_argument("--gpp-et-ratio", action="store_true",
                        help="Plot model mean GPP/ET ratios with the FLUXCOM GPP / MODIS ET reference")
    parser.add_argument("--offline", action="store_true", help="Use only cached downloads")
    args = parser.parse_args(argv)
    if args.gpp_et_ratio:
        if args.lai_et_ratio or args.summary or args.aggregate:
            parser.error("Run --gpp-et-ratio separately from other comparison options")
        summarize_gpp_et_ratios(args.cache, args.output_dir, args.region_name, args.models, args.offline)
        return
    if args.lai_et_ratio:
        if args.summary or args.aggregate:
            parser.error("Run --lai-et-ratio separately from --summary and --aggregate")
        summarize_lai_et_ratios(args.cache, args.output_dir, args.region_name, args.models, args.offline)
        return
    if args.aggregate:
        if args.models:
            parser.error("Use --summary --models first to normalize the selected models, then --aggregate")
        if args.summary:
            summarize_biases(args.cache, args.output_dir, args.region_name)
        summarize_aggregate_scores(args.output_dir, args.region_name)
        return
    if args.summary:
        summarize_biases(args.cache, args.output_dir, args.region_name, args.models)
        return
    benchmark = BENCHMARKS[args.benchmark]
    args.field = args.field or benchmark["field"]
    if args.columns < 1:
        parser.error("--columns must be positive")
    explicit_limits = any(value is not None for value in (args.limit, args.vmin, args.vmax))
    if (benchmark["variable"] == "LAI" and args.field == "bias_map_of_lai"
            and not explicit_limits):
        args.limit = 3.0
    try:
        colour_normalization(np.array([0.0, 1.0]), args.field, args.limit, args.vmin, args.vmax)
    except ValueError as error:
        parser.error(str(error))
    benchmark_url = urljoin(args.url.rstrip("/") + "/", benchmark["path"])
    cache = args.cache / urlparse(args.url).path.strip("/").replace("/", "_")
    if args.benchmark != "modis-lai":
        cache = cache / args.benchmark
    cache.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with requests.Session() as session:
        if not args.offline:
            username = os.environ.get("ILAMB_USERNAME") or input("ILAMB username: ")
            password = os.environ.get("ILAMB_PASSWORD") or getpass.getpass("ILAMB password: ")
            session.auth = (username, password)
        links = discover_models(session, benchmark_url, cache, args.offline, benchmark["source"])
        available = {link[len(benchmark["source"]) + 1:-3]: link for link in links}
        models = args.models or list(available)
        unknown = set(models) - set(available)
        if unknown:
            parser.error(f"Unknown models: {sorted(unknown)}; available: {list(available)}")
        geometry = (read_geometry(args.region_geojson) if args.region_geojson
                    else amazon_geometry(args.cache, args.offline))
        results = []
        for model in models:
            link = available[model]
            print(f"Loading {model}", flush=True)
            path = download(session, urljoin(benchmark_url, link), cache / link, args.offline)
            data = subset_field(path, args.field, geometry)
            data.attrs.update(model=model, source_url=urljoin(benchmark_url, link),
                              region=args.region_name, benchmark=args.benchmark)
            results.append((model, data))
    slug = "".join(character if character.isalnum() else "_"
                   for character in args.region_name.lower()).strip("_")
    stem = args.output_dir / f"{benchmark['slug']}_{args.field}_{slug}"
    subset_dir = args.output_dir / f"{stem.name}_subsets"
    subset_dir.mkdir(exist_ok=True)
    for model, data in results:
        data.to_netcdf(subset_dir / f"{model}.nc")
    boundary = args.output_dir / f"{slug}.geojson"
    boundary.write_text(json.dumps({"type": "Feature", "properties": {
        "name": args.region_name,
        "source": str(args.region_geojson) if args.region_geojson else BASIN_URL,
        "mask_method": "Grid-cell centre within polygon",
    }, "geometry": mapping(geometry)}))
    colour_min, colour_max = plot_models(
        results, geometry, args.field, args.region_name, stem, args.columns,
        args.limit, args.vmin, args.vmax, benchmark["source"], benchmark["variable"],
    )
    manifest = {"ilamb_url": args.url, "field": args.field, "models": models,
            "benchmark": args.benchmark, "benchmark_path": benchmark["path"],
            "source": benchmark["source"], "variable": benchmark["variable"],
            "units": results[0][1].attrs.get("units", "1"),
                "region": args.region_name, "boundary": boundary.name,
                "mask_method": "Grid-cell centre within polygon", "colour_limit": args.limit,
                "colour_min": colour_min, "colour_max": colour_max}
    stem.with_suffix(".json").write_text(json.dumps(manifest, indent=2))
    print(f"Saved {stem.with_suffix('.png')} and PDF; {len(results)} model subsets", flush=True)


if __name__ == "__main__":
    main()