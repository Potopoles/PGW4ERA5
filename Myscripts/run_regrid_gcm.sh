#!/bin/bash

filename=$1

input_data_dir=/home/dargueso/BDY_DATA/CMIP6/PGW4ERA/Clim_Deltas/Amon
output_data_dir=/home/dargueso/BDY_DATA/CMIP6/PGW4ERA/regridded_on_ERA5

while read line; do
echo "GCM: $line"
gcm_name=$line

python step_02_preproc_deltas.py -i ${input_data_dir}/${gcm_name}/ -o ${output_data_dir}/${gcm_name} -e era5_daily_sfc_20230301.nc   -v tas regridding

done < $filename