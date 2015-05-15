#! /usr/bin/env python
from matplotlib import use
use('agg')
import pyart
from matplotlib import pyplot as plt
import matplotlib
import sys
import numpy as np
from time import time
from netCDF4 import num2date, date2num
import datetime
import os


def get_closest(file_list, target_dateobj, fmt, skip):
    dateobjs = []
    valid_dates = []
    for i in range(len(files)):
        try:
            dateobj = datetime.datetime.strptime(file_list[i][skip::], fmt)
            valid_dates.append(file_list[i])
            dateobjs.append(dateobj)
        except:
            pass

    b_d = target_dateobj

    def func(x):
        d =  x
        delta =  d - b_d if d > b_d else datetime.timedelta.max
        return delta

    good_file = valid_dates[np.where(np.array(dateobjs) == min(dateobjs, key = func))[0]]
    return good_file


def do_grid_map_gates_to_grid(radar):
    grids = pyart.map.grid_from_radars(
         radar, grid_shape=(36, 1041, 1041),
        grid_limits=((0, 17000.0),(-900000, 900000), (-1000000, 1000000)),
        gridding_algo="map_gates_to_grid",
        weighting_function='BARNES',
        grid_origin = [35, -98.5])
    return grids


if __name__ == "__main__":
    for t in np.arange(0,55, 5):
        print(t)
        desired_date = datetime.datetime(2015,5, 13, 19, t)
        outdir = sys.argv[2]
        #files = os.listdir('/data/ok_rain/temp/level2/raw/KVNX')
        filenames = []
        top_dir = sys.argv[1]
        sub_dirs = os.listdir(top_dir)
        for ddir in sub_dirs:
            try:
                files = os.listdir(top_dir + '/' + ddir)
                filenames.append(top_dir + \
                    ddir + '/' + get_closest(files, desired_date, '%Y%m%d_%H%M', 5))
            except:
                pass
        print(filenames)
        print(len(filenames))
        radars = []
        for fl in filenames:
            print('reading ', fl)
            try:
                radars.append( pyart.io.read(fl) )
            except:
                print('fail')
        print('gridding')
        t1 = time()
        grid = do_grid_map_gates_to_grid(radars)
        print time() - t1
        display = pyart.graph.GridMapDisplay(grid)

        # create the figurei
        font = {'size': 16}
        matplotlib.rc('font', **font)
        fig = plt.figure(figsize=[15, 8])

        # panel sizes
        map_panel_axes = [0.05, 0.05, .4, .80]
        x_cut_panel_axes = [0.55, 0.10, .4, .30]
        y_cut_panel_axes = [0.55, 0.50, .4, .30]
        colorbar_panel_axes = [0.05, 0.90, .4, .03]

        # parameters
        level = 2
        vmin = -8
        vmax = 64
        lat = 35
        lon = -98.5

        # panel 1, basemap, radar reflectivity and NARR overlay
        ax1 = fig.add_axes(map_panel_axes)
        display.plot_basemap(lon_lines = np.arange(-104, -93, 2) )
        display.plot_grid('reflectivity', level=level, vmin=vmin, vmax=vmax)
        display.plot_crosshairs(lon=lon, lat=lat)

        # plot the reanalysis on the basemap
        # colorbar
        cbax = fig.add_axes(colorbar_panel_axes)
        display.plot_colorbar(cax=cbax)

        # panel 2, longitude slice.
        ax2 = fig.add_axes(x_cut_panel_axes)
        display.plot_longitude_slice('reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax)
        ax2.set_xlabel('Distance from SGP CF (km)')

        # panel 3, latitude slice
        ax3 = fig.add_axes(y_cut_panel_axes)
        display.plot_latitude_slice('reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax)

        # add a title
        slc_height = grid.axes['z_disp']['data'][level]
        dts = num2date(grid.axes['time']['data'], grid.axes['time']['units'])
        datestr = dts[0].strftime('%H:%M Z on %Y-%m-%d')
        title = 'Sliced at ' + str(slc_height) + ' meters at ' + datestr
        fig.text(0.5, 0.9, title)
        sstr = desired_date.strftime('%Y%m%d_%H%M%S')
        print(sstr)
        plt.savefig(outdir + '/diag_'+sstr+'.png')


