#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     PGW for ERA5 preprocessing of climate deltas
authors		    Before 2022: original developments by Roman Brogli
                Since 2022:  upgrade to PGW for ERA5 by Christoph Heim 
"""
##############################################################################
import os, argparse
import xarray as xr
import numpy as np
from pathlib import Path
from functions import filter_data, regrid_lat_lon
from settings import (
    i_debug,
    file_name_bases,
    LON_ERA, LAT_ERA,
    TIME_GCM, PLEV_GCM, LON_GCM, LAT_GCM,
    i_use_xesmf_regridding,
)
##############################################################################

## input arguments
parser = argparse.ArgumentParser(description =
            'PGW for ERA5: Preprocess GCM data before modifying ' +
            'the ERA5 files. The main PGW routine (step_03_apply_to_era.py) ' +
            'requires GCM climate delta files (SCEN-HIST) for ' +
            'ta,hur,ua,va,zg,hurs,tas as well ' +
            'as the GCM HIST climatology file for ps. This script here by ' +
            'default preprocesses both the SCEN-HIST and HIST files ' +
            'for the variables it is run for. Both are thus required for ' +
            'every variable being processed. The script ' +
            'looks for the inputs files using the naming convention ' +
            '${var_name}_${file_name_base}.nc, where the ${file_name_base} ' +
            'for the SCEN-HIST and the HIST files ' +
            'can be set in settings.py, among other things. If this script ' +
            'is used to preprocess daily GCM data, one can run it twice and '+
            'store the intermediate results: once '+
            'for processing_step "smoothing" and once for "regridding". ' +
            'More details are given below.')

# processing step to perform during script execution
parser.add_argument('processing_step', type=str, 
            choices=['smoothing','regridding'],
            help='Possible processing steps are: ' +
            'smoothing: [For daily climate deltas, a smoothing of ' +
            'the annual cycle should be applied. For monthly ' +
            'climate deltas this is not necessary.] ' +
            'regridding: [If the climate deltas are not on the same ' +
            'horizontal grid as ERA5, they can be regridded here. '+
            'WARNING: The default interpolation routine ' +
            '(i_use_xesmf_regridding = 0) assumes a regular '+
            '(thus non-rotated) lat/lon grid for ' +
            'input (GCM data) and output (ERA5 data) grids! ' +
            'If this is not the case for the GCM data, using the ' +
            'xESMF package for regridding may help ' +
            '(i_use_xesmf_regridding = 1). However, such cases have not ' +
            'been tested in detail and may require code adjustments in the ' +
            'function "regrid_lat_lon" in "functions.py".]')

# input directory
parser.add_argument('-i', '--input_dir', type=str,
            help='Directory with input GCM data files (SCEN-HIST, HIST) ' +
            'for selected processing step.')

# output directory
parser.add_argument('-o', '--output_dir', type=str,
            help='Directory where the preprocessed output GCM data files ' +
            'for the selected processing step should be stored.')

# target ERA5 example file to take grid information
parser.add_argument('-e', '--era5_file_path', type=str, default=None,
            help='Path to example ERA5 file ' +
            'from which to take grid information for regridding.')

# variable(s) to process
parser.add_argument('-v', '--var_names', type=str,
            help='Variable names (e.g. tas) to process. Separate ' +
            'multiple variable names with "," (e.g. tas,ta). Default is ' +
            'to process all required variables ta,hur,ua,va,zg,hurs,tas,ps.',
            default='ta,hur,ua,va,zg,hurs,tas,ps')


args = parser.parse_args()
print(args)
##############################################################################

# make sure required input arguments are set.
if args.input_dir is None:
    raise ValueError('Input directory (-i) is required.')
if args.output_dir is None:
    raise ValueError('Output directory (-o) is required.')
if (args.processing_step == 'regridding') and (args.era5_file_path is None):
    raise ValueError('era5_file_path is required for regridding step.')

# create output directory
Path(args.output_dir).mkdir(exist_ok=True, parents=True)

# set up list of variable names
var_names = args.var_names.split(',')
print('Run {} for variable names {}.'.format(
        args.processing_step, var_names))


##############################################################################
# iterate over all variables to preprocess
for var_name in var_names:
    print(var_name)
    if var_name == 'ps':
        clim_periods = ['HIST','SCEN-HIST']
    else:
        clim_periods = ['SCEN-HIST']
    # iterate over the two types of GCM data files
    # (HIST and SCEN-HIST)
    for clim_period in clim_periods:

        var_file_name = file_name_bases[clim_period].format(var_name)

        inp_file = os.path.join(args.input_dir, var_file_name)
        out_file = os.path.join(args.output_dir, var_file_name)


        ## smoothing
        if args.processing_step == 'smoothing':

            filter_data(inp_file, var_name, out_file)

        ## regridding
        elif args.processing_step == 'regridding':
            
            ds_gcm = xr.open_dataset(inp_file)
            ds_era5 = xr.open_dataset(args.era5_file_path)

            ds_gcm = regrid_lat_lon(ds_gcm, ds_era5, var_name,
                                    method='bilinear',
                                    i_use_xesmf=i_use_xesmf_regridding)

            ds_gcm.to_netcdf(out_file)

