#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     Auxiliary functions for PGW for ERA5
authors		Before 2022:    original developments by Roman Brogli
                Since 2022:     upgrade to PGW for ERA5 by Christoph Heim 
                2022:           udpates by Jonas Mensch
"""
##############################################################################
import os,math,warnings
import pandas as pd
import xarray as xr
import numpy as np
from pyvista import PolyData
from pyproj import Geod

from numba import njit
from datetime import datetime,timedelta
from constants import CON_RD,CON_G,CON_MW_MD
from settings import (
    i_debug,
    i_use_xesmf_regridding,
    file_name_bases,
    TIME_ERA, LEV_ERA, HLEV_ERA, LON_ERA, LAT_ERA,
    TIME_GCM, PLEV_GCM, LON_GCM, LAT_GCM,
    LAT_GCM_OCEAN, LON_GCM_OCEAN, TIME_GCM_OCEAN,
)
if i_use_xesmf_regridding:
    import xesmf as xe

## TODO DEBUG
#from matplotlib import pyplot as plt
##############################################################################


##############################################################################
##### ARBITRARY FUNCTIONS
##############################################################################
def dt64_to_dt(dt64):
    """
    Converts a numpy datetime64 object to a python datetime object 
    Input:
      date - a np.datetime64 object
    Output:
      DATE - a python datetime object
    source: 
      https://gist.github.com/blaylockbk/1677b446bc741ee2db3e943ab7e4cabd
    """
    timestamp = ((dt64 - np.datetime64('1970-01-01T00:00:00Z'))
                 / np.timedelta64(1, 's'))
    return(datetime.utcfromtimestamp(timestamp))


##############################################################################
##### PHYSICAL COMPUTATIONS
##############################################################################

def specific_humidity_to_vapor_pressure(hus, pa):
    """
    Compute vapor pressure
    https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
    """
    vapp = hus * pa / (CON_MW_MD + 0.378 * hus)
    return(vapp)

def vapor_pressure_to_specific_humidity(vapp, pa):
    """
    Compute vapor pressure
    https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
    """
    hus = CON_MW_MD * vapp / (pa - (1 - CON_MW_MD) * vapp)
    return(hus)

def saturation_vapor_pressure_water_or_ice(pa, ta, water=True):
    """
    Compute saturation vapor pressure over water or ice 
    according to ECMWF IFS documentation (7.93).
    """
    T0 = 273.16     # K
    if water:
        a1 = 611.21 # Pa
        a3 = 17.502 # 1
        a4 = 32.19  # K
    else: # over ice
        a1 = 611.21 # Pa
        a3 = 22.587 # 1
        a4 = -0.7   # K
    vapp_sat_x = a1 * np.exp(a3 * (ta - T0) / (ta - a4))
    return(vapp_sat_x)

def saturation_vapor_pressure_water_and_ice(pa, ta):
    """
    Compute saturation vapor pressure according to ECMWF IFS documentation (7.92).
    """
    T0 = 273.16     # K
    Ti = 250.16     # K
    alpha = xr.full_like(ta, np.nan)
    alpha = xr.where(ta >= T0, 1, alpha)
    alpha = xr.where(ta <= Ti, 0, alpha)
    alpha = xr.where((ta < T0) & (ta > Ti), np.power((ta - Ti)/(T0 - Ti), 2.), alpha)
    vapp_sat = (
        alpha * saturation_vapor_pressure_water_or_ice(pa, ta, water=True) +
        (1 - alpha) * saturation_vapor_pressure_water_or_ice(pa, ta, water=False)
    )
    return(vapp_sat)

def specific_to_relative_humidity(hus, pa, ta):
    """
    Compute relative humidity from specific humidity
    following the implementation in the IFS model.
    """
    hur = (
        specific_humidity_to_vapor_pressure(hus, pa) / 
        saturation_vapor_pressure_water_and_ice(pa, ta)
    ) * 100
    return(hur)

def relative_to_specific_humidity(hur, pa, ta):
    """
    Compute relative humidity from specific humidity
    following the implementation in the IFS model.
    """
    vapp = hur / 100 * saturation_vapor_pressure_water_and_ice(pa, ta)
    hus = vapor_pressure_to_specific_humidity(vapp, pa)
    return(hus)


def integ_geopot(pa_hl, zgs, ta, hus, level1, p_ref):
    """
    Integrate ERA5 geopotential from surfce to a reference pressure
    level p_ref.
    """
    ## take log half-level pressure difference (located at full levels)
    # make sure pressure is not exactly zero because of ln
    pa_hl = pa_hl.where(pa_hl > 0, 0.0001)
    dlnpa = np.log(pa_hl).diff(
                dim=HLEV_ERA, 
                label='lower').rename({HLEV_ERA:LEV_ERA})

    # create geopotential array and fill with surface geopotential
    phi_hl = zgs.expand_dims(dim={HLEV_ERA:level1}).copy()

    # compute virtual temperature
    tav = ta * (1 + 0.61 * hus)

    ## integrate over model half levels
    for l in sorted(tav[LEV_ERA].values, reverse=True):
        # geopotential at full level
        phi_hl.loc[{HLEV_ERA:l}] = (
                phi_hl.sel({HLEV_ERA:l+1}) +
                (CON_RD * tav.sel({LEV_ERA:l}) * dlnpa.sel({LEV_ERA:l}))
        )

            
    phi_hl = phi_hl.transpose(TIME_ERA, HLEV_ERA, LAT_ERA, LON_ERA)

    ## integrate from last half level below reference pressure
    ## up to reference pressure
    # determine level below reference pressure
    p_diff = pa_hl - p_ref
    p_diff = p_diff.where(p_diff >= 0, np.nan)
    try:
        ind_ref_star = p_diff.argmin(dim=HLEV_ERA)
    except ValueError :
        raise ValueError("p_ref locally lies below the surface. Please set a lower reference pressue (p_ref_inp) in settings.py")
    
    hl_ref_star = p_diff[HLEV_ERA].isel({HLEV_ERA:ind_ref_star})
    # get pressure and geopotential of that level
    p_ref_star = pa_hl.sel({HLEV_ERA:hl_ref_star})
    phi_ref_star = phi_hl.sel({HLEV_ERA:hl_ref_star})

    # finally interpolate geopotential to reference
    # pressure level
    phi_ref = (
            phi_ref_star -
            (CON_RD * tav.sel({LEV_ERA:hl_ref_star-1})) * 
            #(CON_RD * tav.sel({LEV_ERA:hl_ref_star-1},method='nearest')) * 
            (np.log(p_ref) - np.log(p_ref_star))
    )

    # remove multi-dimensional coordinates
    if HLEV_ERA in phi_ref.coords:
        del phi_ref[HLEV_ERA]
    if LEV_ERA in phi_ref.coords:
        del phi_ref[LEV_ERA]
    if PLEV_GCM in phi_ref.coords:
        del phi_ref[LEV_ERA]

    return(phi_ref)


##############################################################################
##### CLIMATE DELTA COMPUTATION AND INTERPOLATION
##############################################################################
def load_delta(delta_input_dir, var_name, era5_date_time, 
               target_date_time=None,
               name_base=file_name_bases['SCEN-HIST']):
    """
    Load a climate delta and if target_date_time is given,
    interpolate it to that date and time of the year.
    """
    ## full climate delta (either daily or monthly)
    full_delta = xr.open_dataset(os.path.join(delta_input_dir,
                            name_base.format(var_name)))
    ## convert time values to standard Pandas Datetimes
    ## this is to catch error arising from cftime.DatetimeNoLeap time format
    ## (https://stackoverflow.com/questions/54462798/cftime-datetimenoleap-object-fails-to-convert-with-pandas-to-datetime)
    ## that is used by some models (e.g. GFDL, but not MPI-ESM1-2-HR)
    ## but only do this if we do not already have a datetime index
    if not isinstance(
        full_delta.indexes['time'], 
        pd.core.indexes.datetimes.DatetimeIndex
    ):
        # ignore warning with calendar conversion
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # convert time to datetime index to convert calendar to standard pandas format
            # to convert non-leap-year calendar to standard calendar
            full_delta = full_delta.assign_coords(
                time=full_delta.indexes['time'].to_datetimeindex()
            )

    ## remove leap year february 29th if in delta
    leap_day = None
    for tind,dt64 in enumerate(full_delta.time.values):
        dt = dt64_to_dt(dt64)
        if (dt.month == 2) and (dt.day == 29):
            leap_day = dt64
    if leap_day is not None:
        full_delta = full_delta.drop_sel(time=leap_day)

    ## if climate delta should be interpolated to a specific time
    if target_date_time is not None:
        # replace delta year values with year of current target_date_time
        for i in range(len(full_delta.time)):
            full_delta.time.values[i] = dt64_to_dt(
                        full_delta.time.values[i]).replace(    #? FIX (added .values[])
                                year=target_date_time.year)

        # find time index of climate delta before target time
        # (and implement periodicity if necessary)
        is_before = (full_delta.time.values <= 
                    np.datetime64(target_date_time))

        # target time is within year 
        # --> periodicity not necessary
        if np.sum(is_before) > 0:
            ind_before = np.argwhere(is_before)[-1].squeeze()
            before = full_delta.isel(time=ind_before)

        # target time is before the first delta time step 
        # --> periodicity necessary
        else:
            ind_before = -1
            before = full_delta.isel(time=ind_before)
            before.time.values = dt64_to_dt(
                        before.time.values).replace(    #? FIX (added .values[])
                                year=target_date_time.year-1)

        # find time index of climate delta after target time
        # (and implement periodicity if necessary)
        is_after =(full_delta.time.values >= 
                    np.datetime64(target_date_time)) 

        # target time is within year 
        # --> periodicity not necessary
        if np.sum(is_after) > 0:
            ind_after = np.argwhere(is_after)[0].squeeze()
            after = full_delta.isel(time=ind_after)

        # target time is after the last delta time step 
        # --> periodicity necessary
        else:
            ind_after = 0
            after = full_delta.isel(time=ind_after)
            after.time.values = dt64_to_dt(
                        after.time.values).replace(    #? FIX (added .values[])
                                year=target_date_time.year+1)

        # if target time is exactly contained in climate delta
        # just take that value (from "before" which is arbitrary)
        if ind_before == ind_after:
            delta = before[var_name].expand_dims(dim='time', axis=0)

        # if interpolation is necessary, concate "before" and 
        # "after" and interpolate 
        else:
            full_delta = xr.concat([before, after], dim='time')
            # interpolate in time and select variable
            delta = full_delta[var_name].interp(time=target_date_time, 
                                        method='linear', 
                                    ).expand_dims(dim='time', axis=0)

        # make sure time is in the same format as in ERA5 file
        # ERA5 has "seconds since xyz" while delta has np.datetime64
        delta['time'] = era5_date_time

    ## if full climate delta should be returned without 
    ## time interpolation
    else:
        delta = full_delta[var_name]

    return(delta)


def load_delta_interp(delta_input_dir, var_name, target_P,
                    era5_date_time, target_date_time,
                    ignore_top_pressure_error=False):
    """
    Does the following:
        - load a climate delta
        - for specific variables (ta and hur) also load surface
          climate delta,
          as well as HIST surface pressure. This is to extend
          the 3D climate deltas with surface values which makes
          the interpolation to the ERA5 model levels more precise.
        - vertically interpolate climate deltas to ERA5 model levels
    """
    delta = load_delta(delta_input_dir, var_name, 
                        era5_date_time, target_date_time)

    ## for ta and hur variables also load climate delta for near surface
    ## values and the historical surface pressure to improve vertical
    ## interpolation
    if var_name in ['ta','hur']:
    #if var_name in []:
        sfc_var_name = var_name + 's'
        delta_sfc = load_delta(delta_input_dir, sfc_var_name, 
                            era5_date_time, target_date_time)
        ps_hist = load_delta(delta_input_dir, 'ps', 
                            era5_date_time, target_date_time,
                            name_base=file_name_bases['HIST'])
    else:
        delta_sfc = None
        ps_hist = None

    # interpolate climate delta onto ERA5 model levels
    delta = vert_interp_delta(delta, target_P, delta_sfc, ps_hist,
                            ignore_top_pressure_error)
    return(delta)


def replace_delta_sfc(source_P, ps_hist, delta, delta_sfc):
    """
    In the 3D climate deltas, replace the value just below
    the surface by the surface climate delta value and insert
    it at HIST surface pressure. This improves the precision
    of the climate deltas during interpolation to the ERA5 model levels.
    All 3D climate delta values below the historical surface pressure
    are set to the surface value (constant extrapolation). This is
    because within the orography the GCM climate delta is assumed
    to be incorrect.
    """
    out_source_P = source_P.copy()
    out_delta = delta.copy()
    if ps_hist > np.max(source_P):
        sfc_ind = len(source_P) - 1
        out_source_P[sfc_ind] = ps_hist
        out_delta[sfc_ind] = delta_sfc
    elif ps_hist < np.min(source_P):
        raise ValueError()
    else:
        sfc_ind = np.max(np.argwhere(ps_hist > source_P))
        out_delta[sfc_ind:] = delta_sfc
        out_source_P[sfc_ind] = ps_hist
    return(out_source_P, out_delta)


def vert_interp_delta(delta, target_P, delta_sfc=None, ps_hist=None,
                       ignore_top_pressure_error=False):
    """
    Vertically interpolate climate delta onto ERA5 model levels.
    If delta_sfc and ps_hist are given, surface values will
    be inserted into the 3D climate delta at the height of
    the surface pressure. This gives a more precise interpolation.
    Climate delta values below the surface are set to the surface
    climate delta because below the surface, the GCM climate delta
    is considered unreliable and thus constant extrapolation
    seems more reasonable.
    """

    # sort delta dataset from top to bottom (pressure ascending)
    delta = delta.reindex(
                {PLEV_GCM:list(reversed(delta[PLEV_GCM]))})

    # create 4D source pressure with GCM pressure levels
    source_P = delta[PLEV_GCM].expand_dims(
                    dim={LON_GCM:delta[LON_GCM],
                         LAT_GCM:delta[LAT_GCM],
                         TIME_GCM:delta[TIME_GCM]}).transpose(
                                TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM)

    ## if surface values are given, replace them at the
    ## level of the surface pressure
    if delta_sfc is not None:
        source_P, delta = xr.apply_ufunc(
                replace_delta_sfc, source_P, 
                ps_hist, 
                delta, delta_sfc,
                input_core_dims=[[PLEV_GCM],[],[PLEV_GCM],[]],
                output_core_dims=[[PLEV_GCM],[PLEV_GCM]],
                vectorize=True)
        source_P = source_P.transpose(TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM)
        delta = delta.transpose(TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM)

    # make sure all arrays contain the required dimensions
    if source_P.dims != (TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM):
        raise ValueError()
    if delta.dims != (TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM):
        raise ValueError()
    if target_P.dims != (TIME_ERA, LEV_ERA, LAT_ERA, LON_ERA):
        raise ValueError()

    # make sure there is no extrapolation at the model top
    # unless these levels are anyways not important for the user
    # and they manually set ignore_top_pressure_error=True
    if np.min(target_P) < np.min(source_P):
        if not ignore_top_pressure_error:
            raise ValueError('ERA5 top pressure is lower than '+
                             'climate delta top pressure. If you are ' +
                             'certain that you do not need the data ' +
                             'beyond to upper-most pressure level of the ' +
                             'climate delta, you can set the flag ' +
                             '--ignore_top_pressure_error and re-run the ' +
                             'script.')
                             

    # run interpolation
    delta_interp = interp_logp_4d(delta, source_P, target_P,
                        extrapolate='constant')
    return(delta_interp)


def interp_logp_4d(var, source_P, targ_P, extrapolate='off',
                   time_key=None, lat_key=None, lon_key=None):
    """
    Interpolate 3D array in vertical (pressure) dimension using the
    logarithm of pressure.
    extrapolate:
        - off: no extrapolation
        - linear: linear extrapolation
        - constant: constant extrapolation
        - nan: set to nan
    """
    if extrapolate not in ['off', 'linear', 'constant', 'nan']:
        raise ValueError('Invalid input value for "extrapolate"')

    # set default values (ERA5 values) for dimension keys
    if time_key is None: time_key = TIME_ERA
    if lat_key is None: lat_key = LAT_ERA
    if lon_key is None: lon_key = LON_ERA

    #print(var.shape)
    #print(source_P.shape)
    #print(targ_P.shape)

    if ( (var.shape[0] != source_P.shape[0]) or
         (var.shape[0] != targ_P.shape[0])   ):
         raise ValueError('Time dimension of input files is inconsistent!')
    if ( (var.shape[2] != source_P.shape[2]) or
         (var.shape[2] != targ_P.shape[2])   ):
         raise ValueError('Lat dimension of input files is inconsistent!')
    if ( (var.shape[3] != source_P.shape[3]) or
         (var.shape[3] != targ_P.shape[3])   ):
         raise ValueError('Lon dimension of input files is inconsistent!')

    targ = xr.zeros_like(targ_P)
    tmp = np.zeros_like(targ.values)
    interp_1d_for_timelatlon(var.values,
                np.log(source_P).values,
                np.log(targ_P).values, 
                tmp,
                len(targ_P[time_key]), len(targ_P[lat_key]),    
                len(targ_P[lon_key]),
                extrapolate)
    targ.values = tmp
    return(targ)

@njit()
def interp_1d_for_timelatlon(orig_array, src_p, targ_p, interp_array,
                        ntime, nlat, nlon, extrapolate):
    """
    Vertical interpolation helper function with numba njit for 
    fast performance.
    Loop over time lat and lon dimensions and interpolate each column
    individually
    extrapolate:
        - off: no extrapolation
        - linear: linear extrapolation
        - constant: constant extrapolation
        - nan: set to nan
    """
    for time_ind in range(ntime):
        for lat_ind in range(nlat):
            for lon_ind in range(nlon):
                src_val_col = orig_array[time_ind, :, lat_ind, lon_ind]
                src_p_col = src_p[time_ind, :, lat_ind, lon_ind]
                targ_p_col = targ_p[time_ind, :, lat_ind, lon_ind]
                
                if src_p_col[-1] < src_p_col[0]:
                    raise ValueError('Source pressure values must be ascending!')
                if targ_p_col[-1] < targ_p_col[0]:
                    raise ValueError('Target pressure values must be ascending!') 

                # call 1D interpolation function for current column
                interp_col = interp_extrap_1d(src_p_col, src_val_col, 
                                            targ_p_col, extrapolate)
                interp_array[time_ind, :, lat_ind, lon_ind] = interp_col


@njit()
def interp_extrap_1d(src_x, src_y, targ_x, extrapolate):
    """
    Numba helper function for interpolation of 1d vertical column.
    Does constant extrapolation which is used for the climate deltas.
    extrapolate:
        - off: no extrapolation
        - linear: linear extrapolation
        - constant: constant extrapolation
        - nan: set to nan
    """
    targ_y = np.zeros(len(targ_x))
    for ti in range(len(targ_x)):
        i1 = -1
        i2 = -1
        require_extrap = False
        for si in range(len(src_x)):
            ty = np.nan
            # extrapolate lower end
            if (si == 0) and src_x[si] > targ_x[ti]:
                if extrapolate == 'linear':
                    i1 = si
                    i2 = si + 1
                elif extrapolate == 'constant':
                    i1 = si
                    i2 = si
                require_extrap = True
                break
            # exact match
            elif src_x[si] == targ_x[ti]:
                i1 = si
                i2 = si
                break
            # upper src_x found (interpolation)
            elif src_x[si] > targ_x[ti]:
                i1 = si - 1
                i2 = si
                break
            # we are still smaller than targ_x[ti] 
            else:
                pass

        # extrapolate upper end
        if i1 == -1:
            if extrapolate == 'linear':
                i1 = len(src_x) - 2
                i2 = len(src_x) - 1
            elif extrapolate == 'constant':
                i1 = len(src_x) - 1 
                i2 = len(src_x) - 1
            require_extrap = True

        # raise value if extrapolation is required but not enabled.
        if require_extrap and extrapolate == 'off':
            raise ValueError('Extrapolation deactivated but data '+
                             'out of bounds.')

        # interpolate/extrapolate values
        if require_extrap and extrapolate == 'nan':
            targ_y[ti] = np.nan
        else:
            if i1 == i2:
                targ_y[ti] = src_y[i1]
            else:
                targ_y[ti] = (
                    src_y[i1] + (targ_x[ti] - src_x[i1]) * 
                    (src_y[i2] - src_y[i1]) / (src_x[i2] - src_x[i1])
                )

    return(targ_y)


def determine_p_ref(p_min_era, p_min_pgw, p_ref_opts, p_ref_last=None):
    """
    Find lowest GCM pressure level among p_ref_opts that lies 
    above the minimum pressure (e.g. currently set to 90% surface pressure)
    in both ERA and PGW climate.
    Also ensure that during the iterations, no reference pressure level 
    at lower altitude than during last iterations is used. This is to
    prevent the iteration algorithm to oscillate between two reference
    pressure levels and not converge.
    """
    for p in p_ref_opts:
        if (p_min_era > p) & (p_min_pgw > p):
            if p_ref_last is None:
                return(p)
            else:
                return(min(p, p_ref_last))




##############################################################################
##### SMOOTHING OF ANNUAL CYCLE FOR DAILY CLIMATE DELTAS
##############################################################################
def filter_data(annualcycleraw, variablename_to_smooth, outputpath):

	"""
	This function performs a temporal smoothing of an annual timeseries 
    (typically daily resolution) using a spectral filter 
    (Bosshard et al. 2011).

	Input:
		Input 1: Path to a netcdf file of the annual cycle to be smoothed. 
        Normally this is the change in a specific variable between 
        two simulations (e.g. warming). 
        Can be 4 or 3 dimensional, where the time is one dimension 
        and the others are space dimensions.
		Input 2: The name of the variable within the given netcdf file
		Input 3: Path to the output file
		
	Output:
		A netcdf file containing the smoothed annual cycle. 
	"""	

	Diff = xr.open_dataset(annualcycleraw
                )[variablename_to_smooth].squeeze()
	coords = Diff.coords

	print('Dimension that is assumed to be time dimension is called: ', 
            Diff.dims[0])
	print('shape of data: ', Diff.shape)

	Diff = Diff.data

	#create an array to store the smoothed timeseries
	#Diff_smooth=np.zeros_like(Diff, dtype=np.float32) 

	if len(Diff.shape) == 4:
		times = Diff.shape[0] 
		levels = Diff.shape[1]
		ygrids = Diff.shape[2]
		xgrids = Diff.shape[3]
	elif len(Diff.shape) == 3:
		times = Diff.shape[0]
		ygrids = Diff.shape[1]
		xgrids = Diff.shape[2]
		levels = 0
	else:
		sys.exit('Wrong dimensions of input file should be 3 or 4-D')


	if len(Diff.shape) == 4:
        #loop over levels to smooth the timeseries on every level
		for i in range(levels):
			for yy in range(ygrids):
				for xx in range(xgrids):	
                    #reconstruct the smoothed timeseries using function below
					Diff[:,i,yy,xx] = harmonic_ac_analysis(Diff[:,i,yy,xx]) 



	if len(Diff.shape) == 3:		
		for yy in range(ygrids):
			for xx in range(xgrids):	
            #dump the smoothed timeseries in the array on the original level
				Diff[:,yy,xx] = harmonic_ac_analysis(Diff[:,yy,xx]) 
			

	print('Done with smoothing')

	#del Diff

	Diff = xr.DataArray(Diff, coords=coords, name=variablename_to_smooth)
	Diff.to_netcdf(outputpath, mode='w')


def harmonic_ac_analysis(ts):
    """
    Estimation of the harmonics according to formula 12.19 -
    12.23 on p. 264 in Storch & Zwiers

    Is incomplete since it is only for use in surrogate smoothing 
    --> only the part of the formulas that is needed there

    Parameters
    ----------
    ts :  1D numpy array
        Passes a time series

    Returns
    -------
        hcts : a reconstructed smoothed timeseries 
                (the more modes are summed the less smoothing)
        mean : the mean of the timeseries (needed for reconstruction)
    """

    if np.any(np.isnan(ts) == True): #if there are nans, return nans
        smooths = np.full_like(ts, np.nan) #sys.exit('There are nan values')
        return smooths
    else:
        #calculate the mean of the timeseries (used for reconstruction)
        mean = ts.mean() 

        lt = len(ts) #how long is the timeseries?
        P = lt

        #initialize the output array. 
        #we will use at max 4 modes for reconstruction 
        #(for peformance reasons, it can be increased)
        hcts = np.zeros((4,lt))

        timevector=np.arange(1,lt+1,1)	#timesteps used in calculation	

        #a measure that is to check that the performed calculation 
        # is justified.
        q = math.floor(P/2.) 

        #create the reconstruction timeseries, mode by mode 
        #(starting at 1 until 5, if one wants more smoothing 
        #this number can be increased.)
        for i in range(1,4): 
            if i < q: #only if this is true the calculation is valid
            
                #these are the formulas from Storch & Zwiers
                bracket = 2.*math.pi*i/P*timevector
                a = 2./lt*(ts.dot(np.cos(bracket))) 
                #dot product (Skalarprodukt) for scalar number output!
                b = 2./lt*(ts.dot(np.sin(bracket))) 
                
                #calculate the reconstruction time series
                hcts[i-1,:] = a * np.cos(bracket) + b * np.sin(bracket) 
            
            else: #abort if the above condition is not fulfilled. In this case more programming is needed.
                sys.exit('Whooops that should not be the case for a yearly '+
                'timeseries! i (reconstruction grade) is larger than '+
                'the number of timeseries elements / 2.')

        smooths = sum(hcts[0:3,:]) + mean
        return smooths




##############################################################################
##### BILINEAR REGRIDDING
##############################################################################
def regrid_lat_lon(ds_gcm, ds_era5, var_name,
                    method='bilinear', i_use_xesmf=0):
    """
    Method to do lat/lon bilinear interpolation for periodic or non-periodic
    grid either with xesmf (i_use_xesmf), or with an xarray-only 
    implementation if the xesmf package is not installed.
    
    Parameters
    ----------
    
    ds_gcm: xr.array
        Passes the original grid and data from the gcm
    ds_era5 : xr.array
        Passes the target grid data to interpolate to
    var_name : string
        Passes the name of the variable to interpolate
    
    Returns
    -------
    ds_regridded : np.ndarray
        ds_gcm interpolated to the ds_era5 grid
    """

    if method != 'bilinear':
        NotImplementedError()

    targ_lon = ds_era5[LON_ERA]
    targ_lat = ds_era5[LAT_ERA]

    ## determine if GCM data set is periodic
    dlon_gcm = np.median(np.diff(ds_gcm[LON_GCM].values))
    dlat_gcm = np.median(np.diff(ds_gcm[LAT_GCM].values))
    if (dlon_gcm + np.max(ds_gcm[LON_GCM].values) - 
                  np.min(ds_gcm[LON_GCM].values)) >= 359.9:
        periodic_lon = True
        if i_debug >= 1:
            print('Regridding: Use periodic boundary conditions for GCM ' +
                    'input data as it ' +
                    'appears to be periodic in longitudinal ' +
                    'direction.')
    else:
        periodic_lon = False


    #### IMPLEMENTATION WITH XESMF
    ##########################################################################
    ## XESMF sometimes alters the exact values of the latitude coordinate
    ## a little bit which was found to be problematic. Therefore, there is an
    ## xarray-only implmenetation below.
    if i_use_xesmf:
        ds_in = ds_gcm
        regridder = xe.Regridder(ds_in, ds_era5, "bilinear", 
                                 periodic=periodic_lon)
        print(regridder)
        ds_out = regridder(ds_in[var_name])
        ds_gcm = ds_out.to_dataset(name=var_name)
        # keep attributes of variables and coordinates
        for field in [var_name, TIME_GCM, PLEV_GCM, LAT_GCM, 
                      LON_GCM, 'height']:
            if field in ds_in:
                ds_gcm[field].attrs = ds_in[field].attrs
        # keep global attributes
        ds_gcm.attrs = ds_in.attrs

    #### IMPLEMENTATION WITH XARRAY
    ##########################################################################
    ## The results should be identical to XESMF
    ## except for tiny differences that appear to originate from
    ## numerical precision.
    else:
        #### LATITUDE INTERPOLATION
        ######################################################################
        ## make sure latitude is increasing with index
        if (ds_gcm[LAT_GCM].isel({LAT_GCM:0}).values > 
            ds_gcm[LAT_GCM].isel({LAT_GCM:-1}).values):
            if i_debug >= 1:
                print('Regridding: GCM data has opposite ' +
                      'order of latitude. Apply reindexing.')
            # flip latitude dimension
            ds_gcm = ds_gcm.reindex(
                    {LAT_GCM:list(reversed(ds_gcm[LAT_GCM]))})

        ## If GCM dataset almost reaches poles, add a pole grid point
        ## with the zonal average of the values closest to the pole
        if (np.max(targ_lat.values) + dlat_gcm > 89.9) and ds_gcm[LAT_GCM].max()<90.0 :
            north = ds_gcm.isel({LAT_GCM:-1})
            north[LAT_GCM].values = 90
            north[var_name] = north[var_name].mean(dim=[LON_GCM])
            ds_gcm = xr.concat([ds_gcm,north], dim=LAT_GCM)
        if (np.min(targ_lat.values) - dlat_gcm < -89.9) and ds_gcm[LAT_GCM].min()>-90.0 :
            south = ds_gcm.isel({LAT_GCM:0})
            south[LAT_GCM].values = -90
            south[var_name] = south[var_name].mean(dim=[LON_GCM])
            ds_gcm = xr.concat([south,ds_gcm], dim=LAT_GCM)

        ## make sure there is no extrapolation to the North and South
        if ( (np.max(targ_lat.values) > np.max(ds_gcm[LAT_GCM].values)) |
             (np.min(targ_lat.values) < np.min(ds_gcm[LAT_GCM].values))):
            print('GCM lat: min {} max {}'.format(
                            np.min(ds_gcm[LAT_GCM].values),
                            np.max(ds_gcm[LAT_GCM].values))) 
            print('ERA5 lat: min {} max {}'.format(
                            np.min(targ_lat.values),
                            np.max(targ_lat.values))) 
            raise ValueError('ERA5 dataset extends further North or South ' +
                              'than GCM dataset!. Perhaps consider using ' +
                              'ERA5 on a subdomain only if global coverage ' +
                              'is not required?') 
        ## run interpolation
        ds_gcm = ds_gcm.interp({LAT_GCM:targ_lat})

        #### LONGITUDE INTERPOLATION
        ######################################################################
        ### Implement periodic boundary conditions
        ### This is also a check and fix for a different longitude format 
        ### e.g. GCM -180:180 while ERA5 0:360 (or vice versa)
        if periodic_lon:
        
            lon_above = ds_gcm.assign_coords({LON_GCM:ds_gcm[LON_GCM] + 360}).isel({LON_GCM:0})
            lon_below = ds_gcm.assign_coords({LON_GCM:ds_gcm[LON_GCM] - 360}).isel({LON_GCM:-1})

            if np.max(targ_lon.values) > np.max(ds_gcm[LON_GCM].values):
                ds_gcm = xr.concat([ds_gcm, lon_above], dim=LON_GCM)
            if np.min(targ_lon.values) < np.min(ds_gcm[LON_GCM].values):
                ds_gcm = xr.concat([lon_below, ds_gcm], dim=LON_GCM)

        ## make sure there is no extrapolation to the East and West
        if ( (np.max(targ_lon.values) > np.max(ds_gcm[LON_GCM].values)) |
             (np.min(targ_lon.values) < np.min(ds_gcm[LON_GCM].values))):
            print('GCM lon: min {} max {}'.format(
                            np.min(ds_gcm[LON_GCM].values),
                            np.max(ds_gcm[LON_GCM].values))) 
            print('ERA5 lon: min {} max {}'.format(
                            np.min(targ_lon.values),
                            np.max(targ_lon.values))) 
            raise ValueError('ERA5 dataset extends further East or West ' +
                              'than GCM dataset!. Perhaps consider using ' +
                              'ERA5 on a subdomain only if global coverage ' +
                              'is not required?') 

        ## run interpolation
        ds_gcm = ds_gcm.interp({LON_GCM:targ_lon})

    ## test for NaN
    #if np.sum(np.isnan(ds_gcm[var_name])).values > 0:
    #    raise ValueError('NaN in GCM dataset after interpolation.')

    return(ds_gcm)

def nan_ignoring_interp(da_era5_land_fr, da_delta, kernel_radius, sharpness):
    """Point cloud based NaN-ignoring interpolation from unstructured to regular lon/lat grid

    This function takes an xarray 2D array with individual lon/lat coordinates and 
    uses pyvista.interpolate() to interpolate that grid to a regular lon/lat grid. This function
    will ignore all NaN values present in the unstructured grid.

    Parameters
    ----------

    da_era5_land_fr : xr.array
        Passes the ERA5 land fraction and grid data
    da_delta : xr.array
        Passes the unstructured array from which to interpolate
    radius : scaler
        Passes a maximal distance for interpolation 
    
    Returns
    -------
    new_era5_grid : np.ndarray
        da_delta interpolated to the ERA5 grid
    """

    #-------------------
    #PREPROCESSING GCM DATA
    #-------------------

    #Construct lon and lat arrays for point cloud by flattening them into an array
    gcm_lat_raw = da_delta.coords[LAT_GCM_OCEAN].values.reshape(-1)
    gcm_lon_raw = da_delta.coords[LON_GCM_OCEAN].values.reshape(-1)
    gcm_val_raw = da_delta.values.reshape(-1)

    # Change to -180,180 range. This simplifies the wrap-around from 0-360
    for i in range(len(gcm_lon_raw)):
        if gcm_lon_raw[i] > 180:
            gcm_lon_raw[i] -= 360

    #Remove NaN values and save NaN mask
    nan_mask_curv = ~np.isnan(gcm_val_raw)
    gcm_val = gcm_val_raw[nan_mask_curv]
    gcm_lon = gcm_lon_raw[nan_mask_curv]
    gcm_lat = gcm_lat_raw[nan_mask_curv]


    # Convert lon/lat degrees to lon/lat kilometers

    #First save the sign of the lon and lat array
    sign_map_lat = np.sign(gcm_lat)
    sign_map_lon = np.sign(gcm_lon)
    #Specify which Earth model to use
    geod = Geod(ellps="WGS84")
    #Convert to km based scale. First two return arguments can be ignored as they are not relevant for this task
    #and this function is still faster then alternatives due to numpy vectorizing
    del1, del2, gcm_lat_meter = geod.inv(gcm_lon, np.zeros(len(gcm_lat)), gcm_lon, gcm_lat)
    del1, del2, gcm_lon_meter = geod.inv(np.zeros(len(gcm_lat)), gcm_lat, gcm_lon, gcm_lat)
    del1, del2, lon_offset    = geod.inv(np.zeros(len(gcm_lat)), gcm_lat, np.ones(len(gcm_lat))*180, gcm_lat)
    
    #Add the signs back to the distances, so that we can have minus coords
    gcm_lat_meter = np.multiply(gcm_lat_meter, sign_map_lat)
    gcm_lon_meter = np.multiply(gcm_lon_meter, sign_map_lon)
    
    #Implement Boundary points
    gcm_val_bd = np.empty(len(gcm_val)*3)
    gcm_lon_bd = np.empty(len(gcm_lon_meter)*3)
    gcm_lat_bd = np.empty(len(gcm_lat_meter)*3)
    
    gcm_val_bd[:] = np.tile(gcm_val, 3)
    gcm_lat_bd[:] = np.tile(gcm_lat_meter, 3)
    gcm_lon_bd[:] = np.tile(gcm_lon_meter, 3)

    gcm_lon_bd[:len(gcm_lon_meter)] -= lon_offset * 2
    gcm_lon_bd[2*len(gcm_lon_meter):] += lon_offset * 2


    #Create grid matrix for further use in the interpolation scheme
    #3rd dim is needed for the interpolation library but can be ignored
    curv_points = np.zeros((gcm_lat_bd.shape[0], 3))
    curv_points[:,0] = gcm_lat_bd
    curv_points[:,1] = gcm_lon_bd

    #-------------------
    #PREPROCESSING ERA5 DATA
    #-------------------

    #Extract Era5 values from xarray
    era5_val_raw = da_era5_land_fr.values.reshape(-1)
    era5_lat_raw = da_era5_land_fr.coords[LAT_ERA].values
    era5_lon_raw = da_era5_land_fr.coords[LON_ERA].values

    # Change to -180,180 range if needed
    for i in range(len(era5_lon_raw)):
        if era5_lon_raw[i] > 180:
            era5_lon_raw[i] -= 360

    #Create flatten lon and lat arrays
    era5_lat_flatten = np.repeat(era5_lat_raw, len(era5_lon_raw))
    era5_lon_flatten = np.tile(era5_lon_raw, len(era5_lat_raw))

    #Repeat the conversion to km
    sign_map_lat = np.sign(era5_lat_flatten)
    sign_map_lon = np.sign(era5_lon_flatten)
    del1, del2, era5_lat_meter = geod.inv(era5_lon_flatten, np.zeros(len(era5_lat_flatten)), 
            era5_lon_flatten, era5_lat_flatten)
    del1, del2, era5_lon_meter = geod.inv(np.zeros(len(era5_lon_flatten)), 
            era5_lat_flatten, era5_lon_flatten, era5_lat_flatten) 
    era5_lat_meter = np.multiply(era5_lat_meter,sign_map_lat)
    era5_lon_meter = np.multiply(era5_lon_meter,sign_map_lon)

    #Create grid matrix for further use in the interpolation scheme
    #3rd dim is needed for the interpolation library but can be ignored
    era5_points = np.zeros((era5_lat_meter.shape[0], 3))
    era5_points[:,0] = era5_lat_meter
    era5_points[:,1] = era5_lon_meter

    # Create masks pure_ocean and frac_ocean for use of the masking land values later on
    land_map = np.where(era5_val_raw > 0.7)[0]

    #-------------------
    #Computing the interpolation
    #-------------------
    #Set up pyvista unstructured grids for both gcm and era5 data
    grid = PolyData(era5_points)
    points = PolyData(curv_points)
    #Fill gcm sst into the grid
    points['values'] = gcm_val_bd
    #Core interpolation: mask land values with empty and only consider points in a certain radius
    grid = grid.interpolate(
        points, 
        null_value=np.nan,
        radius=kernel_radius, 
        sharpness=sharpness
    )
    #-------------------
    #POSTPROCESSING ERA5 DATA
    #-------------------
    
    #Save result
    result = grid['values']
    #Overlay land mask to negate sst propagation onto land
    result[land_map] = np.nan
    
    # Reshape back to matrix
    return result.reshape((len(era5_lat_raw),len(era5_lon_raw)))

def interp_wrapper(
    origin_grid, 
    target_grid, 
    var_name, 
    i_use_xesmf=0, 
    nan_interp_kernel_radius=300000,
    nan_interp_sharpness=3,
    ):
    """Interpolation wrapper that allows for each variable to be assigned a custom scheme

    This function implements for different variables different kinds of interpoaltion. Default is bi-linear
    interpolation from regular grid to regular grid.

    Parameters
    ----------

    origin_grid : xr.array
        Passes the GCM original grid structure and values
    target_grid : xr.array
        Passes target grid structure to interpolate
    i_use_xesmf = boolean
        Passes booloean to adjust which bi-linear scheme is used, only applicable for athmospheric variables
    radius : scaler
        Passes a maximal distance for interpolation, only applicable for ocean variables 
    
    Returns
    -------
    new_era5_grid : xr.DataSet
        da_delta interpolated to the ERA5 grid
    """
    #Custom interpolation for TOS
    if var_name == 'tos':
        land_fraction = target_grid["FR_LAND"][0,:,:]
        tos_values = origin_grid['tos']
        
        #Interpolate all 12 months indivdually
        result = np.empty((12,len(target_grid[LAT_ERA]), len(target_grid[LON_ERA])))
        for i in range(12):
            result[i,:,:] = nan_ignoring_interp(
                land_fraction, 
                tos_values[i,:,:], 
                kernel_radius=nan_interp_kernel_radius,
                sharpness=nan_interp_sharpness
            )
        
        #Save into xr.Dataset
        ds = xr.Dataset(
            data_vars=dict(
                tos=(["time","lat", "lon"], result),
            ),
            coords=dict(
                lat=(["lat"], target_grid[LAT_ERA].data),
                lon=(["lon"], target_grid[LON_ERA].data),
                time=origin_grid[TIME_GCM]
            ),
            attrs=dict(description="SST on ERA5 grid", units="K", long_name="Sea Surface Temperature"),
        )
    #Default interpolation
    else:
        ds = regrid_lat_lon(origin_grid, target_grid, var_name,
                        method='bilinear',
                        i_use_xesmf=i_use_xesmf)
    return ds

def integrate_tos(tos_field, ts_field, land_frac, ice_frac):
    """Combines TOS and TS temperature as a weighted according to land and ocean contribution
    
    This functions assumes that all fields are on the same grid!

    Parameters
    ----------
    tos_field : 2d.nparray
        Passes TOS field on ERA5 grid
    ts_field : 2d.nparray
        Passes TS field on ERA5 grid
    land_frac : 2d.nparray
        Passes land fraction on ERA5 grid
    ice_field : 2d.nparray
        Passes sea ice fraction on ERA5 grid

    
    Returns
    -------
    combined_ts : 2d.nparray
        Returns combined temperature fields
    """
    #Save dims to reshape the original grid
    dims = tos_field.shape
    #Ocean mask
    ice_frac = ice_frac.reshape(-1)
    tos_field = tos_field.reshape(-1)
    #TODO add cases for true mask

    mask = ~np.isnan(ice_frac) & ~np.isnan(tos_field)
    #Turn all fields into arrays for easier accessing
    ice_frac_masked  = ice_frac[mask]
    land_frac_masked = land_frac.reshape(-1)[mask]
    tos_field_masked = tos_field.reshape(-1)[mask]
    ts_field_masked = ts_field.reshape(-1)[mask]
    #Arrays to save output in
    output = np.ones(len(tos_field.reshape(-1)))
    output[:] = ts_field.reshape(-1)
    
    ts_frac = np.clip(ice_frac_masked + land_frac_masked, 0, 1)
    output[mask] = np.add(np.multiply(ts_frac, ts_field_masked), np.multiply(1-ts_frac, tos_field_masked))
    
    return output.reshape(dims)
