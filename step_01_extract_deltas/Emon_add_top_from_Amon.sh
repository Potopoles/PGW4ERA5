#!/bin/bash

# base directory where output should be stored
out_base_dir=/net/o3/hymet_nobackup/heimc/data/pgw/deltas/native

# name of the GCM to extract data for
gcm_name=MPI-ESM1-2-HR

## type of CMIP6 model output (e.g. monthly or daily, etc.)
## to use
# high-resolution monthly data for only very few GCMs
table_ID=Emon

## select variables to add model top
var_names=(ua va ta zg hur)
#var_names=(hur)

## CMIP experiments to use to compute climate deltas
## --> climate delta = future climatology - ERA climatology
# CMIP experiment to use for ERA climatology 
era_climate_experiment=historical
# CMIP experiment to use for future climatology 
future_climate_experiment=ssp585

# iterate over both experiments and climate delta
experiments=($era_climate_experiment $future_climate_experiment delta)


out_dir=$out_base_dir/$table_ID/$gcm_name

for var_name in ${var_names[@]}; do
    echo "#################################################################"
    echo $var_name
    echo "#################################################################"


    # add Amon model top values to Emon
    if [[ "$table_ID" == "Emon" ]]; then
        for experiment in ${experiments[@]}
        do
            echo $experiment

            #mv $out_dir/${var_name}_${experiment}.nc \
            #    $out_dir/Emon_model_bottom_${var_name}_${experiment}.nc
            cdo sellevel,100000,97500,95000,92500,90000,87500,85000,82500,80000,77500,75000,70000,65000,60000,55000,50000,45000,40000,35000,30000,25000,22500,20000,17500,15000,12500,10000 $out_dir/${var_name}_${experiment}.nc \
                $out_dir/Emon_model_bottom_${var_name}_${experiment}.nc

            Amon_out_base_dir=$out_base_dir/Amon
            Amon_out_dir=$Amon_out_base_dir/$gcm_name
            cdo sellevel,7000,5000,3000,2000,1000,500,100 \
                $Amon_out_dir/${var_name}_${experiment}.nc \
                $out_dir/Amon_model_top_${var_name}_${experiment}.nc
            cdo -O merge \
                $out_dir/Emon_model_bottom_${var_name}_${experiment}.nc \
                $out_dir/Amon_model_top_${var_name}_${experiment}.nc \
                $out_dir/${var_name}_${experiment}.nc
        done
    fi


done

