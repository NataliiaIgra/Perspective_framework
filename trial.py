import netCDF4 as nc
import numpy as np
import xarray as xr
import csv

input = xr.open_dataset('met_em.d02.2018-11-08_06_00_00.nc')
input["UU"].values = np.full(shape=(1, 38, 270, 361), fill_value=-5, dtype=float)
input["VV"].values = np.full(shape=(1, 38, 271, 360), fill_value=-5, dtype=float)
input["TT"].values = np.full(shape=(1, 38, 270, 360), fill_value=300, dtype=float)
input["RH"].values = np.full(shape=(1, 38, 270, 360), fill_value=25, dtype=float)
input["GHT"].values = np.full(shape=(1, 38, 270, 360), fill_value=0, dtype=float)
input.to_netcdf('copy_met_em.d02.2018-11-08_06_00_00.nc')


# input = xr.open_dataset('wrfinput_d01')
# input["W"].values = np.full(shape=(1, 46, 270, 360), fill_value=5, dtype=float)
# input.to_netcdf('copy_wrfinput_d01')

# rootgrp = nc.Dataset("test.nc", "w", format="NETCDF4")
# time = rootgrp.createDimension("Time", 1)
# num_metgrid_levels = rootgrp.createDimension("num_metgrid_levels", 38)
# south_north = rootgrp.createDimension("south_north", 270)
# west_east = rootgrp.createDimension("west_east", 360)
# new_var = rootgrp.createVariable("WW", "f4", ("Time", "num_metgrid_levels", "south_north", "west_east"))

# ncfile = 'met_em.d02.2018-11-08_06_00_00.nc'
# ds = nc.Dataset(ncfile, mode="r")
# print(ds["RH"][0][0])

# newfile = f'C:\PythonProjects\Perspective_Framework\wrfinput_d01'
# ds2 = nc.Dataset(newfile, mode="r")
# lat = ds2.variables['XLAT'][0]
# lon = ds2.variables['XLONG'][0]
# print(len(lat[0]))
# print(ds2["XLONG"][0][269])
# print(ds2.variables)


# with open('employee_file.csv', mode='w', newline="") as employee_file:
#     employee_writer = csv.writer(employee_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#
#     employee_writer.writerow(['LATITUDE', 'LONGITUDE'])
#
#     for j in range(0, 270, 10):
#         for i in range(0, 360, 10):
#             employee_writer.writerow([lat[j][i], lon[j][i]])
#
#     for q in range(0, 360, 10):
#         employee_writer.writerow([lat[269][q], lon[269][q]])
#
#     for p in range(0, 270, 10):
#         employee_writer.writerow([lat[p][359], lon[p][359]])
#
#     employee_writer.writerow([lat[269][359], lon[269][359]])





