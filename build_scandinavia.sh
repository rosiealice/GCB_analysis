#!/usr/bin/env bash
set -euo pipefail

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

output_dir=outputs/scandinavia
common=(--offline
  --region-geojson "$output_dir/scandinavia.geojson"
  --region-name Scandinavia
  --north-limit 73
  --columns 8
  --output-dir "$output_dir")

for benchmark in modis-lai modis-et gleam-et fluxcom-gpp wecann-gpp \
  avh15c1-lai esacci-biomass-change esacci-biomass-stock; do
  case "$benchmark" in
    modis-lai) limits=(--vmin -3 --vmax 3) ;;
    modis-et|gleam-et) limits=(--vmin -2 --vmax 2) ;;
    fluxcom-gpp|wecann-gpp|avh15c1-lai) limits=(--vmin -3 --vmax 3) ;;
    esacci-biomass-change) limits=(--vmin -2 --vmax 2) ;;
    esacci-biomass-stock) limits=() ;;
  esac
  python ilamb_regions.py "${common[@]}" --benchmark "$benchmark" "${limits[@]}"
done

python ilamb_regions.py --offline --summary --region-name Scandinavia --output-dir "$output_dir"
python ilamb_regions.py --offline --aggregate --region-name Scandinavia --output-dir "$output_dir"
python ilamb_regions.py --offline --lai-et-ratio --region-name Scandinavia --output-dir "$output_dir"
python ilamb_regions.py --offline --gpp-et-ratio --region-name Scandinavia --output-dir "$output_dir"
python ilamb_regions.py --offline --lai-gpp-ratio --region-name Scandinavia --output-dir "$output_dir"

marp --pptx --allow-local-files \
  --output "$output_dir/scandinavia_plots.pptx" \
  "$output_dir/scandinavia_plots.md"
marp --pdf --allow-local-files \
  --output "$output_dir/scandinavia_plots.pdf" \
  "$output_dir/scandinavia_plots.md"

destination=/datalake/NS9560K/www/diagnostics/noresm/rosief/GCB_analysis
mkdir -p "$destination"
cp "$output_dir/scandinavia_plots.pptx" \
  "$output_dir/scandinavia_plots.pdf" "$destination"/