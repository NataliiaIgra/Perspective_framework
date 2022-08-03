import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import math
import wrf


def wrf_out_read (dir, time):
    mins= time % 60
    hrs = math.floor(time / 60)
    wrfout = Dataset('{}wrfout_d02_2018-11-08_{:02d}:{:02d}:00'.format(dir, hrs, mins))
    return wrfout

def relax_zone_remover (input, sr):
    output = input
    for _ in range(sr):
        output = np.delete(output, -1, 0)
        output = np.delete(output, -1, 1)
    return output

def wind_plot (data, xcoords, ycoords, field_1, field_2, height_value, color):
    u = wrf.getvar(data, field_1, timeidx=wrf.ALL_TIMES, method='cat', meta=False)
    v = wrf.getvar(data, field_2, timeidx=wrf.ALL_TIMES, method='cat', meta=False)
    u = u[height_value, :, :]
    v = v[height_value, :, :]

    wf = plt.quiver(xcoords[::8,::8], ycoords[::8,::8], u[::8,::8], v[::8,::8], color = color, scale = 8, scale_units = 'xy', pivot = 'tail', width = 0.002) 
    return wf
   

out_time = 45 #Output time to plot in minutes
sr = 4

plt.style.use('wrf')

#wrfout files with reinit
outs_folder_reinit = 'C:/PythonProjects/Perspective_Framework/'
wrfouts_reinit = wrf_out_read(outs_folder_reinit, out_time)


#reading coordiantes
x = wrf.getvar(wrfouts_reinit, 'XLONG', timeidx=wrf.ALL_TIMES, method='cat', meta=False) / 1000   #converting coordinates to km
y = wrf.getvar(wrfouts_reinit, 'XLAT', timeidx=wrf.ALL_TIMES, method='cat', meta=False) / 1000
xf = wrf.getvar(wrfouts_reinit, 'FXLONG', timeidx=wrf.ALL_TIMES, method='cat', meta=False) / 1000   #converting coordinates to km
yf = wrf.getvar(wrfouts_reinit, 'FXLAT', timeidx=wrf.ALL_TIMES, method='cat', meta=False) / 1000
hgt = wrf.getvar(wrfouts_reinit, 'HGT', timeidx=wrf.ALL_TIMES, method='cat', meta=False) / 1000

#reading data to single array
lfn_reinit = wrf.getvar(wrfouts_reinit, 'LFN', timeidx=wrf.ALL_TIMES, method='cat', meta=False)

#removing the relaxation zones of the level-set function
lfn_reinit = relax_zone_remover(lfn_reinit, sr)
xf = relax_zone_remover(xf, sr)
yf = relax_zone_remover(yf, sr)

fig = plt.figure()
ax = plt.subplot2grid((1,1), (0,0))

ax.contour(xf, yf, lfn_reinit, 0, colors='r')
ax.plot([], [], color='r', label='Fire Line ($\phi>0)$')

CS = ax.contourf(x, y, hgt)

#plotting the wind arrows
wf = wind_plot (wrfouts_reinit, x, y, 'ua', 'va', 0, 'w')
plt.quiverkey(wf, 0.7, 0.9, U=5, label=r'$5 \frac{m}{s}$', labelpos='E', coordinates='figure', color = 'k')

cbar = plt.colorbar()
cbar.set_label('Terrain Height (m)')

ax.tick_params(direction='in')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

plt.ylabel('Y (km)')
plt.xlabel('X (km)')
plt.legend()
plt.xlim(0, 5)
plt.xticks(np.arange(0, 5.5, 0.5))
plt.ylim(0, 5)
plt.yticks(np.arange(0, 5.5, 0.5))


plt.show()


