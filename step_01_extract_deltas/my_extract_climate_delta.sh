#!/bin/bash
##############################################################################
## Template script to extract climate deltas from CMIP GCMs.
## Computes climatologies and the climate deltas for specific CMIP 
## experiments and members.
## The script only serves as inspiration for the extraction of the climate
## deltas and it is not generally valid and has to be adjusted 
## for specific use case.

##### IMPORTANT NOTES FOR Emon DATA USERS:
## FIRST
## hur is not available in the Emon output and has to be computed from hus.
## This is not excat because it matters if hur is computed from hus
## on the model levels every time the GCM writes an output or whether it is
## computed on the monthly-aggregated output on pressure levels.
## Therefore, the hur is used from the coarse resolution output group Amon.
## However, the Amon output is vertically interpolated to the higher resolution
## of Emon using information from the computed Emon hur=f(hus,ta)
## data. This helps to do a "better informed" vertical interpolation of the
## coarse Amon hur to the higher resolved Emon grid.
## This steps are done by the convert_Emon_hus_to_hur.py script which is
## called automatically in the current implementation.
## Consequently, Amon must have been exctracted before Emon can be!
## SECOND
## After computing hur for Emon, the delta still has to be computed
## by running the script again only for var_names=(hur).
## (turn off i_exctract_vars but keep i_compute_delta)
## THIRD
## The Emon data does not reach up as far as the Amon data
## Thus, after completing this script, run add_Emon_model_top.sh
## to add the model top data from Amon to the Emon fields.
## FOURTH
## These Emon fixes are still a bit experimental and not very user-friendly.
## Sorry for that!
##############################################################################
Usage="
##############################################################################
## Template script to extract climate deltas from CMIP GCMs.
## Computes climatologies and the climate deltas for specific CMIP 
## experiments and members.
## The script only serves as inspiration for the extraction of the climate
## deltas and it is not generally valid and has to be adjusted 
## for specific use case.

-h to see this message
-s <dataset_name> to get info on a specific data set (TODO)
<> to get deltas for Amon dataset (should be removed later on)
<dataset_name> to get the desired delta from this dataset
    Curently supported are: Amon, Omon, Emon, day, CFday"

while getopts 'h' flag; do
  case "${flag}" in
    h) echo "$Usage" ; exit 0;
  esac
done

# USER SETTINGS
##############################################################################
# base directory where cmip6 data is stored
cmip_data_dir=/home/dargueso/BDY_DATA/CMIP6/PGW4ERA/
# base directory where output should be stored
out_base_dir=/home/dargueso/BDY_DATA/CMIP6/PGW4ERA/Clim_Deltas

filename=$2

# name of the GCM to extract data for

while read line; do

echo "GCM: $line"
gcm_name=$line


## CMIP experiments to use to compute climate deltas
## --> climate delta = future climatology - ERA climatology
# CMIP experiment to use for ERA climatology 
era_climate_experiment=historical
# CMIP experiment to use for future climatology 
future_climate_experiment=ssp585


## type of CMIP6 model output (e.g. monthly or daily, etc.)
## to use
# standard monthly output
#table_ID=Amon
## high-resolution monthly data for only very few GCMs
#table_ID=Emon
## standard daily output
#table_ID=day
## CFMIP daily output
#table_ID=CFday
## ocean monthly output
#table_ID=Omon

##Experimental addition to use command-line args to select type
##Can (should) be done in a foor loop so that it runs for multiple data sets at the same time
if [[ "$1" == "Amon" ]]; then
    table_ID=$1
elif [[ "$1" == "Omon" ]]; then
    table_ID=$1
elif [[ "$1" == "Emon" ]]; then
    table_ID=$1
elif [[ "$1" == "day" ]]; then
    table_ID=$1
elif [[ "$1" == "CFday" ]]; then
    table_ID=$1
elif [ $# -eq 0 ]; then
    echo "Please specify which dataset the deltas should be extracted from. See -h for more info"
    exit 0
else
    txt="This dataset is currently not supported
Check -h for all supported datasets"
    echo "$txt"
    exit 0
fi

## select variables to extract
if [[ "$table_ID" == "Amon" ]]; then
    var_names=(ts psl tas hurs vas uas ps ua va ta hur zg)
elif [[ "$table_ID" == "day" ]]; then
    var_names=(tas hurs ps ua va ta hur zg)
elif [[ "$table_ID" == "Emon" ]]; then
    var_names=(ua va ta hus zg)
elif [[ "$table_ID" == "CFday" ]]; then
    var_names=(ua va ta hur)
elif [[ "$table_ID" == "Omon" ]]; then
    var_names=(tos)
fi
## for Emon, compute hur delta separately after 
## it has been derived from hus
## (turn off i_exctract_vars but keep i_compute_delta)
#var_names=(hur)


# should variables be exctracted for the two climatologies=
i_extract_vars=1
# should climate deltas be computed?
i_compute_delta=1
# for Emon, hur is not available and needs to be approximated
# using the high-resolution hus climatology, as well as the hur climatology
if [[ "$table_ID" == "Emon" ]]; then
    i_convert_hus_to_hur=1
else
    i_convert_hus_to_hur=0
fi

## subdomain for which to extract GCM data
## should be either global (0,360,-90,90)
## or anything larger than ERA5 subdomain
## except for storage and performance reasons, there is no benefit of
## using a subdomain.
box=0,360,-90,90
# subdomain
#box=-74,40,-45,35
#box=-73,37,-42,34

# select appropriate cdo time aggregation command
# depending if input data is monthly or daily.
if [[ "$table_ID" == "day" ]]; then
    cdo_agg_command=ydaymean
else 
    cdo_agg_command=ymonmean
fi

# iterate over both experiments to extract data
experiments=($era_climate_experiment $future_climate_experiment)

##############################################################################

out_dir=$out_base_dir/$table_ID/$gcm_name
echo $out_dir
mkdir -p $out_dir

for var_name in ${var_names[@]}; do
    echo "#################################################################"
    echo $var_name
    echo "#################################################################"

    for experiment in ${experiments[@]}; do
        echo "#######################################"
        echo $experiment
        echo "#######################################"

        if [[ $i_extract_vars == 1 ]]; then

            # data folder hierarchy for CMIP6
            inp_dir=$cmip_data_dir/$experiment/$table_ID/$var_name/$gcm_name

            # start of the CMIP6 file names
            file_name_base=${table_ID}_${gcm_name}_${experiment}

            ## overwrite old data
            #rm $out_dir/${var_name}_${experiment}.nc

            ## compute ERA climatology
            if [[ "$experiment" == "$era_climate_experiment" ]]; then
                # extract full time series
                cdo sellonlatbox,$box \
                    -selyear,1985/2014 \
                    -cat \
                    $inp_dir/${var_name}_${file_name_base}_*.nc \
                    $inp_dir/${var_name}_${file_name_base}_*.nc \
                    $out_dir/${var_name}_${experiment}_full.nc

            ## compute future experiment climatology
            elif [[ "$experiment" == "$future_climate_experiment" ]]; then
                # extract full time series
                cdo sellonlatbox,$box \
                    -selyear,2070/2099 \
                    -cat \
                    $inp_dir/${var_name}_${file_name_base}_*.nc \
                    $out_dir/${var_name}_${experiment}_full.nc
            fi

            # aggregate to yearly monthly/daily means
            # in principal this could be done directly during extraction
            # step above. However, for Emon computation of hur from hus
            # this should be done on the basis of monthly values not
            # with the mean annual cycle. Due to this, the full time
            # series is stored as well.
            cdo $cdo_agg_command \
                $out_dir/${var_name}_${experiment}_full.nc \
                $out_dir/${var_name}_${experiment}.nc
        fi

        ## convert hus to hur if required
        if [[ $i_convert_hus_to_hur == 1 ]]; then
            if [[ "$var_name" == "hus" ]]; then
                echo Convert Emon hus to hur using Amon hur data.

                Amon_out_dir=$out_base_dir/Amon/$gcm_name

                python Emon_convert_hus_to_hur.py  \
                    $out_dir/hus_${experiment}_full.nc \
                    $out_dir/ta_${experiment}_full.nc \
                    $out_dir/hur_${experiment}_full.nc \
                    -a $Amon_out_dir/hur_${experiment}_full.nc
                    
                # aggregate to yearly monthly/daily means
                cdo $cdo_agg_command \
                    $out_dir/hur_${experiment}_full.nc \
                    $out_dir/hur_${experiment}.nc
            fi
        fi

    done

    ## compute delta (future climatology - ERA climatology)
    if [[ $i_compute_delta == 1 ]]; then
        cdo sub $out_dir/${var_name}_$future_climate_experiment.nc \
                $out_dir/${var_name}_$era_climate_experiment.nc \
                $out_dir/${var_name}_delta.nc
    fi

done
done < $filename

# link surface fields from Amon because they are not available
# in Emon
if [[ "$table_ID" == "Emon" ]]; then
    cd $out_dir
    ln -s $out_base_dir/Amon/$gcm_name/*s_*.nc .
    cd -
fi

