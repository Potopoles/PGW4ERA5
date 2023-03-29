#!/usr/bin/env python
"""
@File    :  check_completeness_simulation.py
@Time    :  2023/02/20 18:11:40
@Author  :  Daniel Argüeso
@Version :  1.0
@Contact :  d.argueso@uib.es
@License :  (C)Copyright 2023, Daniel Argüeso
@Project :  EPICC
@Desc    :  None
"""

import os
import argparse
from glob import glob
import pandas as pd
import datetime as dt
import xarray as xr
from tqdm.auto import tqdm
import logging


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


if models_str is None:
    with open("list_CMIP6_historical.txt") as f:
        models_hist = f.read().splitlines()
    with open("list_CMIP6_ssp585.txt") as f:
        models_proj = f.read().splitlines()

    models = list(set(models_hist + models_proj))
    # models = ["GFDL-CM4_r1i1p1f1"]  # , "MPI-ESM1-2-HR_r1i1p1f1", "MRI-ESM2-0_r1i1p1f1"]
else:
    models = args.models.split(",")

logfile = "check_CMIP6_PGW.log"
if os.path.exists(logfile):
    os.remove(logfile)

logging.basicConfig(
    filename=logfile,
    filemode="a",
    format="%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
    level=logging.INFO,
)

logging.info("Checking completeness of CMIP6 files for PGW")

logger = logging.getLogger("checkCMI6PGW")

scenarios = {"historical": [1850, 2014], "ssp585": [2015, 2100]}
variables = [
    "hur",
    "hurs",
    "ps",
    "psl",
    "ta",
    "tas",
    "ts",
    "ua",
    "uas",
    "va",
    "vas",
    "zg",
]


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


def main():
    for GCM in models:
        complete = True
        first = True
        for scen in scenarios.keys():
            print(f"{bcolors.HEADER}Checking {GCM} {scen}{bcolors.ENDC}")
            for var in (pbar := tqdm(variables)):
                pbar.set_description(f"Processing {var}")
                path_in = f"./{scen}/{GCM}/"
                finlesin_var = sorted(glob(f"{path_in}/{var}_*"))
                # print(finlesin_var)
                if finlesin_var:
                    finvar = xr.open_mfdataset(
                        finlesin_var, use_cftime=True, combine="by_coords"
                    )

                    # Checking that the period is covered
                    if (finvar.time[0].dt.year > scenarios[scen][0]) or (
                        finvar.time[-1].dt.year < scenarios[scen][1]
                    ):
                        complete = False
                        if first:
                            logger.info(
                                f"{bcolors.FAIL}{scen} {GCM} is not complete{bcolors.ENDC}"
                            )
                            first = False

                        logger.info(f"{GCM} {scen} {var} is not complete")
                    finvar.close()
                else:
                    complete = False
                    if first:
                        logger.info(
                            f"{bcolors.FAIL}{scen} {GCM} is not complete{bcolors.ENDC}"
                        )
                        first = False
                    logger.info(f"{GCM} {scen} {var} is not complete")

        if complete:
            print(f"{bcolors.OKGREEN} {GCM} is complete{bcolors.ENDC}")
            logger.info(f"{bcolors.OKGREEN} {GCM} is complete{bcolors.ENDC}")
        else:
            print(
                f"{bcolors.FAIL}{GCM} is not complete{bcolors.ENDC}"
                f"{bcolors.BOLD} Check log file {logfile} for details{bcolors.ENDC}"
            )


###############################################################################
# __main__  scope
###############################################################################

if __name__ == "__main__":
    raise SystemExit(main())

###############################################################################
