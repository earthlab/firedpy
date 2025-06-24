#!/bin/bash

export FIREDPY_ED_USER=erve3705
export FIREDPY_ED_PWD=luxHizp1032\!\!

# Loop through spatial values
for spatial in {4..5}; do
  # Loop through temporal values
  for temporal in {5..11}; do
	  output_file="$(pwd)/fired_conus_ak_2000_to_2025_${spatial}_${temporal}.gpkg"

    if [ -f "$output_file" ]; then
      echo "Skipping spatial=$spatial, temporal=$temporal — $output_file already exists"
      continue
    fi
    echo "Running firedpy with spatial=$spatial and temporal=$temporal"

    python bin/firedpy.py \
      -spatial "$spatial" \
      -temporal "$temporal" \
      -eco_region_level 3 \
      -land_cover_type 1 \
      -shape_type gpkg \
      --tile_choice b \
      --tile_name conus_ak \
      --daily true \
      -start_year 0 \
      -end_year 0 \
      --full_csv false \
      --n_cores 2 \
      --cleanup false
	  mv "$(pwd)/output/outputs/shapefiles/fired_conus_ak_2000_to_2025_events.gpkg "fired_conus_ak_2000_to_2025_${spatial}_${temporal}.gpkg"
    echo "Finished run with spatial=$spatial and temporal=$temporal"
    echo "----------------------------------------"
  done
done
