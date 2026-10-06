# ILAMB TRENDY regional comparisons

Download LAI, evapotranspiration (ET), and gross primary productivity (GPP)
spatial diagnostics from the authenticated TRENDYv15 ILAMB run, mask them to
a watershed, and compare all available models using one shared colour scale.
Amazon bias maps are available for all 24 models in `GCP_270926a4` against
MODIS and AVH15C1 LAI, MODIS and GLEAMv3.3a ET, and FLUXCOM and WECANN GPP.
CLM-FATES panel titles are always bold. WECANN photosynthesis is represented
by the site's gross primary productivity benchmark.
GLEAM, WECANN and AVH15C1 remain available as standalone maps, but are
excluded from the current mean-bias and aggregate comparisons.

## Results

- [Amazon bias comparison, PNG](outputs/modis_lai_bias_map_of_lai_amazon_basin.png)
- [Amazon bias comparison, PDF](outputs/modis_lai_bias_map_of_lai_amazon_basin.pdf)
- [Run metadata](outputs/modis_lai_bias_map_of_lai_amazon_basin.json)
- [MODIS ET bias, PNG](outputs/modis_et_bias_map_of_et_amazon_basin.png)
- [MODIS ET bias, PDF](outputs/modis_et_bias_map_of_et_amazon_basin.pdf)
- [GLEAMv3.3a ET bias, PNG](outputs/gleam_et_bias_map_of_et_amazon_basin.png)
- [GLEAMv3.3a ET bias, PDF](outputs/gleam_et_bias_map_of_et_amazon_basin.pdf)
- [FLUXCOM GPP bias, PNG](outputs/fluxcom_gpp_bias_map_of_gpp_amazon_basin.png)
- [FLUXCOM GPP bias, PDF](outputs/fluxcom_gpp_bias_map_of_gpp_amazon_basin.pdf)
- [WECANN GPP bias, PNG](outputs/wecann_gpp_bias_map_of_gpp_amazon_basin.png)
- [WECANN GPP bias, PDF](outputs/wecann_gpp_bias_map_of_gpp_amazon_basin.pdf)
- [AVH15C1 LAI bias, PNG](outputs/avh15c1_lai_bias_map_of_lai_amazon_basin.png)
- [AVH15C1 LAI bias, PDF](outputs/avh15c1_lai_bias_map_of_lai_amazon_basin.pdf)
- [Normalized mean-bias comparison, PNG](outputs/normalized_mean_bias_amazon_basin.png)
- [Normalized mean-bias comparison, PDF](outputs/normalized_mean_bias_amazon_basin.pdf)
- [Mean biases and normalization scales, CSV](outputs/normalized_mean_bias_amazon_basin.csv)
- [Aggregate mean-bias scores, PNG](outputs/aggregate_mean_bias_scores_amazon_basin.png)
- [Aggregate mean-bias scores, PDF](outputs/aggregate_mean_bias_scores_amazon_basin.pdf)
- [Category and aggregate scores, CSV](outputs/aggregate_mean_bias_scores_amazon_basin.csv)
- [Model LAI/ET ratios with MODIS reference, PNG](outputs/lai_et_ratio_amazon_basin.png)
- [Model LAI/ET ratios with MODIS reference, PDF](outputs/lai_et_ratio_amazon_basin.pdf)
- [Regional LAI, ET and ratios, CSV](outputs/lai_et_ratio_amazon_basin.csv)
- [Model GPP/ET ratios with observed reference, PNG](outputs/gpp_et_ratio_amazon_basin.png)
- [Model GPP/ET ratios with observed reference, PDF](outputs/gpp_et_ratio_amazon_basin.pdf)
- [Regional GPP, ET and ratios, CSV](outputs/gpp_et_ratio_amazon_basin.csv)

The readable comparison uses a symmetric +/-3 m2/m2 colour scale; values
outside that range saturate, as indicated by the colourbar arrows. Original
bias values, including extremes, are preserved in the saved NetCDF data.
Native-grid basin subsets are saved as individual NetCDF files in
`outputs/modis_lai_bias_map_of_lai_amazon_basin_subsets/`.
The ET comparisons use +/-2 mm/day, and GPP uses +/-3 g/m2/day. These are
display limits, not thresholds for model skill. Each benchmark has its own
NetCDF subset directory and JSON metadata, including the chosen scale and units.

## Run

Python 3.10 or newer:

```sh
python -m pip install -r requirements.txt
python ilamb_regions.py --vmin -3 --vmax 3
```

The tool prompts for the username and password, hiding password input.
Alternatively, provide `ILAMB_USERNAME` and `ILAMB_PASSWORD` through your
environment. Credentials are not stored in source, cache or output files.
The HTML index and complete model NetCDF files are cached, so future plots
need neither network access nor credentials:

```sh
python ilamb_regions.py --offline --vmin -3 --vmax 3
python ilamb_regions.py --offline --models CLM-FATES CLM6.0 E3SM --limit 3
python ilamb_regions.py --offline --field biasscore_map_of_lai
python ilamb_regions.py --offline --field rmsescore_map_of_lai
```

Select the reference dataset and variable with `--benchmark`; the default
remains `modis-lai`. Download a benchmark once without `--offline`, then use
these cached commands:

```sh
python ilamb_regions.py --offline --benchmark modis-et --vmin -2 --vmax 2
python ilamb_regions.py --offline --benchmark gleam-et --vmin -2 --vmax 2
python ilamb_regions.py --offline --benchmark fluxcom-gpp --vmin -3 --vmax 3
python ilamb_regions.py --offline --benchmark wecann-gpp --vmin -3 --vmax 3
python ilamb_regions.py --offline --benchmark avh15c1-lai --vmin -3 --vmax 3
```

| Benchmark | Reference | Default field | Bias units |
| --- | --- | --- | --- |
| `modis-lai` | MODIS LAI | `bias_map_of_lai` | m2/m2 |
| `modis-et` | MODIS evapotranspiration | `bias_map_of_et` | mm/day |
| `gleam-et` | GLEAMv3.3a evapotranspiration | `bias_map_of_et` | mm/day |
| `fluxcom-gpp` | FLUXCOM GPP | `bias_map_of_gpp` | g/m2/day |
| `wecann-gpp` | WECANN photosynthesis (GPP) | `bias_map_of_gpp` | g/m2/day |
| `avh15c1-lai` | AVH15C1 LAI | `bias_map_of_lai` | m2/m2 |

Benchmark-specific caches prevent the identically named MODIS LAI and ET
files from colliding. The original LAI cache and output names are preserved.

`--field` selects an existing variable in the NetCDF `MeanState` group;
an unavailable name reports the file's available map variables. Score maps
use a shared 0-1 scale by default, with higher scores indicating better performance.
Set `--vmin` and `--vmax` together to choose independent minimum and maximum
colour values for whichever variable is selected by `--field`:

```sh
python ilamb_regions.py --offline --field bias_map_of_lai --vmin -3 --vmax 3
python ilamb_regions.py --offline --field bias_map_of_lai --vmin -2 --vmax 4
python ilamb_regions.py --offline --field biasscore_map_of_lai --vmin 0.2 --vmax 1
```

All models share the chosen scale. Out-of-range values saturate in the plot,
and colourbar arrows indicate which ends contain saturated values; saved
NetCDF data are unchanged. Bias and shift bounds must straddle zero to keep
the diverging scale centred there. `--limit 3` remains a shortcut for a
symmetric -3 to +3 bias scale; do not combine it with `--vmin`/`--vmax`.
With no explicit limits, LAI bias defaults to -3 to +3 so a routine rerun
does not restore the extreme-model scale. Other bias fields still use the
full range unless limits are supplied. Only the selected comparison is
generated, with no additional full-range figure. The resolved bounds are
recorded in the metadata.

```sh
python ilamb_regions.py --offline
```

Other regions can be supplied as WGS84 longitude/latitude GeoJSON polygons
(local regions not crossing the antimeridian):

```sh
python ilamb_regions.py --offline \
  --region-geojson path/to/watershed.geojson --region-name "Other watershed"
```

Use `--output-dir` for alternative versions without overwriting a previous
plot. `--url` selects another run with the same benchmark layout. These six
benchmarks are supported; arbitrary ILAMB benchmarks are not auto-discovered.

## LAI/ET Ratio

```sh
python ilamb_regions.py --lai-et-ratio --offline
```

This plots models on the Y-axis and their regional mean LAI divided by
regional mean ET on the X-axis. A dashed vertical line shows the MODIS
observed ratio; models are ordered by distance from that line. CLM-FATES is bold.

The inputs are the existing ILAMB `MeanState/timeint_of_lai` and
`MeanState/timeint_of_et` period-mean fields, not bias maps or normalized
scores. MODIS observations come from the two benchmark NetCDF files. On a
first run with missing benchmark files, omit `--offline` to download them
with the usual authentication; later runs use the cache. Existing LAI and ET
map manifests and cached model files are required. `--models` can select a subset.

The ratio is the quotient of two area-weighted basin means, not the mean
of cell-by-cell ratios. Within each model or observation pair, both means
use exactly the same finite LAI/ET cells inside the watershed. Native grid
coordinates and cell bounds must match within each pair. Mean ET must be
positive. With dimensionless LAI and ET in mm/day, the ratio has units day/mm.

The CSV contains 24 model rows and one observation row, including the original
means, ratio, valid-cell count and total valid area. A JSON sidecar records
the observation value, sources and calculation method. This does not change
the bias or aggregate-score diagrams.

These ratios use ILAMB's supplied period means without additional temporal
harmonization. Model and observation grids and valid footprints can differ.
Agreement in LAI/ET alone does not imply agreement in either LAI or ET, and
the ordering is not an overall skill ranking.

## GPP/ET Ratio

```sh
python ilamb_regions.py --gpp-et-ratio --offline
```

This creates a separate diagram in the same layout as LAI/ET: models on the
Y-axis, ratios of area-weighted regional mean GPP to mean ET on the X-axis,
and a dashed observed reference line. GPP/ET models are ordered by increasing
absolute ratio, not by distance from observations; LAI/ET ordering is unchanged.
The same selected model names are red,
with CLM-FATES bold. Existing LAI/ET outputs are not overwritten.

Model means come from `MeanState/timeint_of_gpp` and `timeint_of_et`. The
observed ratio combines FLUXCOM GPP and MODIS ET benchmark means. The shared
ratio calculation uses jointly valid basin cells within each pair and native
cell areas, with matching grids and positive mean ET required. Units are
g/m2/mm, not a normalized bias score. The CSV contains 24 model rows and an
observation row, including both original means and their ratio.

Omit `--offline` on a first run if the FLUXCOM benchmark file is not cached.
As with LAI/ET, these are ratios of supplied ILAMB period means without
additional temporal harmonization; model and observed footprints can differ,
and ratio agreement alone does not establish agreement in GPP or ET.

## Mean-Bias Comparison

```sh
python ilamb_regions.py --summary --offline
```

This uses three saved bias subsets and their original cached NetCDF files;
there are no downloads or authentication prompts. The plot has models on the
X-axis and normalized mean bias on the Y-axis, with separate coloured markers
for MODIS ET, MODIS LAI and FLUXCOM GPP.
CLM-FATES stays bold. The comparison uses the same 24 models for all three
benchmarks, yielding 72 means. GLEAM ET, WECANN GPP and AVH15C1 LAI are excluded.

Each regional mean is weighted by spherical grid-cell area computed from the
native latitude/longitude bounds, not by cell count or latitude alone. The
comparison grids can have unequal cell widths. Only valid, basin-masked cells
contribute, and original biases are used regardless of map colour saturation.

Normalization is separate for each benchmark: divide each model's mean bias
by the largest absolute mean bias among the plotted models for that benchmark.
Values lie between -1 and +1: negative means underestimation, positive means
overestimation, and zero means no net regional bias. If all models have zero
bias in a benchmark, all normalized values are zero. These are normalized
biases, not ILAMB skill scores or percentage biases. The normalization compares
relative ensemble extremes, not a common physical or uncertainty-based scale.

The CSV retains each mean in its native units, the normalization denominator,
the normalized result, the number of valid cells and their total area.
`--models` selects a smaller model set and recomputes normalization for that
set. The JSON sidecar records the method and input manifests. Signed positive
and negative local errors can cancel in a mean, and different models may have
different valid-data coverage; a small net bias alone does not establish high skill.

## Aggregate Scores

```sh
python ilamb_regions.py --aggregate --offline
```

This reads the saved normalized mean-bias summary and creates another diagram,
without replacing the signed-bias plot. Models are on the Y-axis, ranked by
aggregate score, and scores are on the X-axis. Three coloured markers show the
individual reference scores; hollow black diamonds show the aggregate. Every
marker sits on its model's exact horizontal row. CLM-FATES is bold.

Each benchmark score is one minus the absolute normalized mean bias: one
means zero net regional bias, and zero means the largest absolute mean bias
used for that benchmark's normalization. Both signs are penalized equally.
The ET and LAI categories use MODIS only, and GPP uses FLUXCOM only.
GLEAM, WECANN and AVH15C1 contribute neither points nor scores.
The aggregate
then averages ET, LAI and GPP with equal weights; adding references does not
give a variable extra influence. Absolute values are taken before category
aggregation, preventing positive and negative category biases from cancelling.

The CSV records the three category scores, the aggregate and all three
individual reference scores. These are custom relative mean-bias scores, not official
ILAMB skill scores. Spatial positive and negative errors can still cancel
within the regional mean, so the ranking measures net mean bias rather than
overall model skill or mean absolute spatial error.

Run `--summary --aggregate --offline` to recompute both diagrams. For a model
subset, first run `--summary --models ...`, then `--aggregate`; normalization
uses the selected ensemble recorded in the saved summary.

## Scientific Choices

- Bias is ILAMB's existing `MeanState/bias_map_of_lai`: model minus MODIS,
  dimensionless leaf area per ground area (m2/m2). Red is overestimation;
  blue is underestimation. It is raw bias, not the normalized bias skill score.
- ET and GPP likewise use the existing `MeanState/bias_map_of_et` and
  `MeanState/bias_map_of_gpp` diagnostics: model minus the selected reference.
  Native daily flux units are retained, with no conversion to annual totals.
  The MODIS ET benchmark is `HydrologyCycle/Evapotranspiration/MODIS`, not
  the separate MOD16A2 benchmark. GLEAM uses the site's `GLEAMv3.3a` product.
- The benchmark describes a derived monthly MODIS LAI climatology. This
  tool uses the existing ILAMB period-integrated diagnostic; it does not
  reconstruct dates or recompute bias from monthly model output.
- The basin boundary is the union of nine HydroBASINS South America level-4
  polygons with `MAIN_BAS = 6040007000`, about 5.91 million km2. This is the
  connected Amazon watershed, excluding neighboring coastal basins and
  Tocantins. It is not the Amazon biome or a rectangular Amazon region.
- A cell is retained when its centre is inside the basin polygon; this is
  not fractional-area clipping. Missing observations or model data remain
  missing. Grey means outside-basin or missing data.
- Each field's own latitude/longitude coordinates are used, including the
  ILAMB comparison grid (`lat_`, `lon_`), which may differ between models.
  There is no additional regridding or spatial smoothing. Maps use geographic
  longitude/latitude with aspect correction at the region's central latitude.
- Native NetCDF masks are respected, including implicit default fill values
  that have no explicit `_FillValue` attribute in these ILAMB files.

Boundary source: [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins),
Lehner and Grill (2013). Consult HydroSHEDS' terms before redistributing its
boundary data. Cached downloads and generated outputs are git-ignored.

## Checks

```sh
python -m unittest discover -s . -p 'test_*.py' -v
```

These offline tests cover implicit fill values, coordinate normalization,
polygon holes, missing fields, empty regions, model discovery, and explicit
colour bounds (including asymmetric bias scales and invalid limits), as well
as benchmark labels/units, bold CLM-FATES titles, and bounded download retries.
Summary tests also cover unequal-area weighting, missing values, native bounds,
signed normalization and the three-point-per-model plot.
Aggregate tests cover sign-independent penalties, equal variable weights,
zero-bias scores, ranking, row alignment and excluded-reference independence.
Ratio tests cover area-weighted means, joint valid-cell masks, matching grids,
positive denominators and the observed reference line.