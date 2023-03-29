#!/usr/bin/env python
'''
@File    :  rename_coordinates.py
@Time    :  2023/03/20 16:07:30
@Author  :  Daniel Argüeso
@Version :  1.0
@Contact :  d.argueso@uib.es
@License :  (C)Copyright 2022, Daniel Argüeso
@Project :  PGW4ERA5
@Desc    :  None
'''

import xarray as xr
import argparse
from glob import glob
import os

class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


#####################################################################
#####################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description="PURPOSE: Check the completeness of the CMIP6 files for PGW"
    )

    parser.add_argument(
        "-m",
        "--models",
        dest="models",
        help="Optional input list of models",
        type=str,
        nargs="?",
        default=None,
    )

    args = parser.parse_args()
    return args


args = parse_args()
models_str = args.models
models = args.models.split(",")

input_folder = "/home/dargueso/BDY_DATA/CMIP6/PGW4ERA/Clim_Deltas/"
tableID = "Amon"

#####################################################################
#####################################################################


def main():
    for GCM in models:
        print(f"{bcolors.HEADER}Renaming coordinates {GCM}{bcolors.ENDC}")

        filesin = sorted(glob(f"{input_folder}/{tableID}/{GCM}/*nc"))

        for ifile in filesin:

            fin = xr.open_dataset(ifile)
            fout = fin.rename({'longitude':'lon','latitude':'lat'})

            fout.to_netcdf(f"./{os.path.basename(ifile)}")

            os.replace(f"./{os.path.basename(ifile)}",ifile)

###############################################################################
# __main__  scope
###############################################################################

if __name__ == "__main__":
    raise SystemExit(main())

###############################################################################


