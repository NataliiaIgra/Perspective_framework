import netCDF4 as nc
import numpy as np
import xarray as xr

input = xr.open_dataset('met_em.d02.2018-11-08_06_00_00.nc')
input["UU"].values = np.full(shape=(1, 38, 270, 361), fill_value=-5, dtype=float)
input["VV"].values = np.full(shape=(1, 38, 271, 360), fill_value=-5, dtype=float)
#input.to_netcdf('copy_met_em.d02.2018-11-08_06_00_00.nc')