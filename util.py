import numpy as np 
import matplotlib.pyplot as plt
from datetime import datetime
from glob import glob 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import matplotlib as mpl
import os


def __get_data2(fn): 
    """read the netcdf4 data
    
    Args:
        fn: Path to OCO2 data file

    returns: 
        xco2: XCO2 data
        latitude, longitude: location 
        time: date & time when OCO2 observation is recorded (which is not used when visualizing the data)
    """
    # debug
    ds = nc.Dataset(fn)
    xco2 = ds['xco2'][:]
    latitude = ds['latitude'][:]
    longitude = ds['longitude'][:]
    time = ds['time'][:]
    ds.close()
    return xco2, latitude, longitude, time

def __great_circle(lon1, lat1, lon2, lat2):
    """calculate the great circle distance in km
    
    Args: 
        lon1, lat1: location at point1 
        lon2, lat2: location at point2 
    
    Returns: 
        great distance between point1 and point2 in km
    """
    
    lon1, lat1, lon2, lat2 = map(np.deg2rad, [lon1, lat1, lon2, lat2])
    return 6371 * (
        np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2))
    )

def __gen_filelist(path_dir): 
    """ fetch files under a certain folder
    Args: 
        path_dir: path to the folder where data is placed (e.g., path_dir = '../../data/data_xco2/oco2_LtCO2_*.nc4')

    Returns: 
        flist_all: all the filepaths 
        date_infor_all: date stamp of all the files (datetime.datetime(2014, 9, 6, 0, 0)) 
    
    """
    if path_dir[-1] != '/': 
        path_dir = path_dir + '/'
    flist2 = np.array(glob(path_dir+'*.nc4'))
    date_infor = np.array([  datetime.strptime('20'+ os.path.split(x)[1][-37:-31], '%Y%m%d')  for x in flist2])
    flist_all = flist2[np.argsort(date_infor)]
    date_infor_all = np.sort(date_infor)

    return flist_all, date_infor_all

def __get_useful_file(path_dir, time_yy, years=1, \
    fov_lon=[139.3, 140.4], fov_lat=[35.15, 36.05]): 

    """Find the file that data points are within in the specified FOV. 
    Args: 
        path_dir: path where the data are placed
        time_y: year 
        years: duration (in unit of year)
        fov_lon, fov_lat: the specified FOV ([139.3, 140.4] and [35.15, 36.05] by default)

    output: 
        useful_file: 
    """
    flist_all, date_infor_all = __gen_filelist(path_dir)


    start_time = datetime.strptime('20'+'{}0101'.format(time_yy), '%Y%m%d')
    end_time   = datetime.strptime('20'+'{}1231'.format(time_yy+years-1) , '%Y%m%d')

    need_ind   = np.where((date_infor_all <= end_time) * (date_infor_all >=start_time ))[0]

    useful_file = []        
    for i, fn in enumerate(flist_all[need_ind]): 
        xco2, latitude, longitude, time, = __get_data2(fn) 
        fov_ind = (longitude.data>=fov_lon[0]) * (longitude.data<=fov_lon[1]) * (latitude.data>=fov_lat[0]) * (latitude.data<=fov_lat[1])
        if np.sum(fov_ind) > 0: 
            useful_file.append(fn)

    return useful_file


def __make_meshgrid(useful_file, grid_num, xco2_min, xco2_max, fov_lon, fov_lat): 
    """Mapping the Xco2 observation into a mesh grid array. Values are scaled into (0, 155).  
    Grid with multiple data points in it are filled with the averaged values over the data.

    Args: 
        useful_file: Files with data points within in the specified FOV. 
        grid_num: the number that slicing latitude and longitude.
        [xco2_min, xco2_max]: value range in y-axis. 
        fov_lon, fov_lat: FOV of the map. By default FOV is within a [139.3, 140.4,35.15, 36.05] window

    returns: 
        x, y: x- and y-axis of the mesh grid map.
        z_map: resulting map    
    
    """

    lat_slice = np.linspace(fov_lat[0], fov_lat[1], grid_num)
    lon_slice = np.linspace(fov_lon[0], fov_lon[1], grid_num)
    x, y = np.meshgrid(lon_slice, lat_slice)
    z  = np.empty(x.flatten().shape) 

    z[:] = np.nan

    # fill in the data 
    for fn in useful_file: 
        xco2, latitude, longitude, time, = __get_data2(fn) 
        fov_ind = (longitude.data>=fov_lon[0]) * (longitude.data<=fov_lon[1]) * (latitude.data>=fov_lat[0]) * (latitude.data<=fov_lat[1])
        if np.sum(fov_ind) > 0: 
            color_value = (xco2-xco2_min)*(254.0/(xco2_max-xco2_min))
            for i in np.where(fov_ind)[0]:  
                gcd = __great_circle(longitude[i], latitude[i], x.flatten(), y.flatten())
                gcd_arg = gcd.argmin()
                if np.isnan(z[gcd_arg]) :
                    z[gcd_arg] = color_value[i]
                else: 
                    z[gcd_arg] = np.mean( [color_value[i] ,z[gcd_arg]] )                
    z_map = z.reshape(x.shape)

    return x, y, z_map


def __cbar(cax,cs, xco2_min, xco2_max): 
    # # color bar is only located next to the second axis 
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="5%", pad=1.3)
    cbar = plt.colorbar(cs, cax=cax, orientation='horizontal', ticks=np.linspace(0, 254, 5)) 
    cbar.set_alpha(1)
    cbar.draw_all()     # Avoid stange strips appearing in the colorbar
    # cax_tick_label = cbar.ax.get_xticks()     # # adjust the cax tick labels 
    # cbar.ax.set_xticks(np.linspace(0, 254, 5))
    cbar.ax.set_xticklabels(['{0:.0f}'.format(x)  for x in np.linspace(xco2_min, xco2_max, 5)])
    cax.set_xlabel('XCO2 (ppm)', fontsize=20)

    return


def __fig_map_sigle(ax, fov_lon, fov_lat): 

    """plot a base world map on the spacified axes object.
    Args: 
        ax: 
        fov_lon, fov_lat: FOV 
    
    Returns:
        m: basemap object 

    """

    m = Basemap(llcrnrlon=fov_lon[0], llcrnrlat=fov_lat[0], urcrnrlon=fov_lon[1], urcrnrlat=fov_lat[1],\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc', ax=ax)

    m.drawcoastlines(color='gray',zorder=5)
    m.drawparallels(np.arange(10,90,0.5),labels=[1,1,0,1], zorder=5)
    m.drawmeridians(np.arange(120,160,0.5),labels=[1,1,0,1], zorder=5)
    m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')

    # x_tok, y_tok = m(Tok_lon,Tok_lat )
    # m.plot(x_tok, y_tok, '^r', markersize=10)
    return m


def fig_check_value_range(path_dir, time_yy, years=1, \
    xco2_min=380, xco2_max=435, fov_lon=[139.3, 140.4], fov_lat=[35.15, 36.05] ): 
    
    """ generating a figure for confirming the data range during the year.
    Args: 
        path_dir: path of the stored data file.
        time_yy: year of the file would be fetched. 
        [xco2_min, xco2_max]: value range in y-axis. (in the [380, 435] window by default)
        fov_lon, fov_lat: the specified FOV of the plotting ([139.3, 140.4] and [35.15, 36.05] by default)
    """    
    
    plt.clf()
    plt.close()

    flist_all, date_infor_all = __gen_filelist(path_dir)
    start_time = datetime.strptime('20'+'{}0101'.format(time_yy), '%Y%m%d')
    end_time = datetime.strptime('20'+'{}1231'.format(time_yy+years-1) , '%Y%m%d')
    need_ind = np.where((date_infor_all <= end_time) * (date_infor_all >=start_time ))[0]

    ax = plt.subplot()

    useful_file = []        
    for i, fn in enumerate(flist_all[need_ind]): 
        xco2, latitude, longitude, time, = __get_data2(fn) 
        fov_ind = (longitude.data>=fov_lon[0]) * (longitude.data<=fov_lon[1]) * (latitude.data>=fov_lat[0]) * (latitude.data<=fov_lat[1])
        if np.sum(fov_ind) > 0: 
            useful_file.append(fn)
            legend_timestamp = datetime.strftime(date_infor_all[need_ind[i]],'%Y-%m-%d')
            ax.plot(time[fov_ind], xco2[fov_ind], marker='o', markersize=4 , linestyle='none', label=legend_timestamp)

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.axhline(y=xco2_min, color='k', lw=0.5, ls=':' )
    ax.axhline(y=xco2_max, color='k', lw=0.5, ls=':')

    ax.set_xlabel('time')
    ax.set_ylabel('xco2')

    time_yy = os.path.split(useful_file[0])[1][-37:-31]

    ax.set_title('Data collected during 20{}'.format(time_yy))

    plt.tight_layout()

    figname = 'fig/xcoc2_20{}.png'.format(time_yy)
    plt.savefig(figname, dpi=250)

    print('{} is generated.'.format(figname))

    return 


def fig_point_plot(path_dir, time_yy, years=1 ,\
     xco2_min=380, xco2_max=435, jump=100, fov_lon=[139.3, 140.4], fov_lat=[35.15, 36.05]): 
    
    """ Generating a figure for visualizing 2D distribution of OCO2 observation by simply making a dot plot.
    Args: 
        path_dir: path of the stored data file.
        time_yy: year of the file would be fetched. 
        years: data duration (unit: year)
        jump:only pick up data at a 100-datapoint span (i.e., the 0th, 100th, 200th ... data point) so the the plotting process would not be too costy. 
        [xco2_min, xco2_max]: data value range for visualizing. 
        fov_lon, fov_lat: FOV of the map. By default FOV is within a [139.3, 140.4,35.15, 36.05] window.
    """    
    
    useful_file = __get_useful_file(path_dir, time_yy, years=years)

    plt.close()
    plt.clf()
    
    ax = plt.subplot()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)

    # color map information
    N = 255 ; cmap = plt.get_cmap('jet', N)

    for fn in useful_file: 
        xco2, latitude, longitude, time= __get_data2(fn) 
        xco2.data[xco2.data == 1e+20] = np.nan
        testing_ind = (longitude.data>=125) * (longitude.data<=150) * (latitude.data>=20) * (latitude.data<=60)
        if testing_ind.sum() != 0: 
            color_value = (xco2-xco2_min)*(254.0/(xco2_max-xco2_min))
            testing_ind2 = np.where(testing_ind != False)
            for ind in np.arange(0, len(testing_ind2[0]), jump) :
                i = testing_ind2[0][ind]  
                v = color_value[i]
                lon = longitude.data[i]
                lat = latitude.data[i]
                ax.plot(lon, lat, c=cmap(int(v)), marker='o')
                
    ax.set_xlim(fov_lon)
    ax.set_ylim(fov_lat)
    ax.set_xlabel('Longitude ($^o$)')
    ax.set_ylabel('Latitude ($^o$)')
    ax.set_title('Data collected during 20{}'.format(time_yy))

    norm = mpl.colors.Normalize(vmin=xco2_min, vmax=xco2_max)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, cax=cax)
    cbar.ax.set_ylabel('XCO2 (ppm)')

    plt.tight_layout()

    figname = 'fig/xco2_sat_point_20{}.png'.format(time_yy)
    plt.savefig(figname, dpi=250)

    print('{} is generated.'.format(figname))

    return


def fig_map_meshgrid(path_dir, time_yy, grid_num, years=1, \
     xco2_min=380, xco2_max=435,fov_lon=[139.3, 140.4], fov_lat=[35.15, 36.05]): 
    
    """basically share the same function with `fig_point_plot` but in a mesh grid manner. 
    Grid with multiple data points in it are filled with an averaged value.
    
    Args: 
        path_dir: path to the folder where data is stored (e.g., path_dir = '../../data/data_xco2/oco2_LtCO2_*.nc4').
        time_yy: target year (start).
        grid_num: a number denoting how many slices that latitude and longitude should be divided into.
        years: duration of the data.
        [xco2_min, xco2_max]: value range in y-axis. 
        fov_lon, fov_lat: FOV of the map. By default FOV is within a [139.3, 140.4,35.15, 36.05] window.


    """

    useful_file = __get_useful_file(path_dir, time_yy, years=years)

    x, y, z_map = __make_meshgrid(useful_file, grid_num, xco2_min=xco2_min, xco2_max=xco2_max, fov_lon=fov_lon, fov_lat=fov_lat)
    
    plt.clf()
    plt.close()

    fig, ax = plt.subplots(figsize=(12,12))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=1.3)
    m = Basemap(llcrnrlon=x.min(), llcrnrlat=y.min(), urcrnrlon=x.max(), urcrnrlat=y.max(), \
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc', ax=ax)

    m.drawcoastlines(color='gray',zorder=5)
    m.drawparallels(np.arange(10,90,0.5),labels=[1,1,0,1], zorder=5)
    m.drawmeridians(np.arange(120,160,0.5),labels=[1,1,0,1], zorder=5)
    m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')

    
    # mark Tokyo's location
    # Tok_lat = 35.652832;  Tok_lon = 139.839478 
    # x_tok, y_tok = m(Tok_lon,Tok_lat )
    # m.plot(x_tok, y_tok, '^r', markersize=10)

    x2 = np.linspace(0, m.urcrnrx, z_map.shape[1])
    y2 = np.linspace(0, m.urcrnry, z_map.shape[0])

    xx, yy = np.meshgrid(x2, y2)

    cs = m.pcolormesh(xx, yy, z_map, cmap='jet', zorder=10, alpha=0.5,vmin=0, vmax=254)

    cbar = plt.colorbar(cs, cax=cax)

    # Avoid strips appearing in the colorbar
    cbar.set_alpha(1)
    cbar.draw_all()
    

    # # adjust the cax tick labels 
    cax_tick_label = cax.get_yticks()
    # cax.set_yticklabels(np.linspace(xco2_min, xco2_max, len(cax_tick_label)))
    cax.set_yticklabels(['{0:.0f}'.format(x)   for x in np.linspace(xco2_min, xco2_max, len(cax_tick_label))])
    cax.set_ylabel('XCO2 (ppm)', fontsize=20)

    # time_yy = os.path.split(useful_file[0])[1][-37:-31]
    grid_size = round(100/z_map.shape[0])
    ax.set_title('Data collected during 20{}'.format(time_yy))
    figname = 'fig/xcoc2_map_meshgrid_20{}_{}.png'.format(time_yy, grid_size)


    plt.tight_layout()
    plt.savefig(figname, dpi=250)
    print('{} is generated.'.format(figname))


    return 


def fig_map_meshgrid_2panels(path_dir, time_yy1, time_yy2, grid_num, years=1, \
                             xco2_min=380, xco2_max=435, fov_lon=[139.3, 140.4], fov_lat=[35.15, 36.05]): 


    plt.clf()
    plt.close()

    useful_file1= __get_useful_file(path_dir, time_yy1, years=years, fov_lon=fov_lon, fov_lat=fov_lat)
    useful_file2= __get_useful_file(path_dir, time_yy2, years=years, fov_lon=fov_lon, fov_lat=fov_lat)
    
    x, y, z_map1 = __make_meshgrid(useful_file1, grid_num, xco2_min=xco2_min, xco2_max=xco2_max, fov_lon=fov_lon, fov_lat=fov_lat )
    x, y, z_map2 = __make_meshgrid(useful_file2, grid_num, xco2_min=xco2_min, xco2_max=xco2_max, fov_lon=fov_lon, fov_lat=fov_lat)

    # split the subplots 
    fig, axs = plt.subplots(1,2, figsize=(19,12))

    # plot each subplot
    for ax, z_map in  zip(axs, [z_map1, z_map2] ): 
        m = __fig_map_sigle(ax, fov_lon=fov_lon, fov_lat=fov_lat)

        x2 = np.linspace(0, m.urcrnrx, z_map.shape[1])
        y2 = np.linspace(0, m.urcrnry, z_map.shape[0])
        xx, yy = np.meshgrid(x2, y2)
        cs = m.pcolormesh(xx, yy, z_map, cmap='jet', zorder=10, alpha=0.5,vmin=0, vmax=254)

        grid_size = round(200/z_map.shape[0]) # fov is 200X200 km
    
    # add color bar
    cax = fig.add_axes([0.2, 0.08, 0.5, 0.03]) 
    __cbar(cax,cs, xco2_min, xco2_max)


    # add titles
    if years == 1:
        axs[0].set_title('20{}'.format(time_yy1), fontsize=20)
        axs[1].set_title('20{}'.format(time_yy2), fontsize=20 )
    else: 
        axs[0].set_title('20{}-20{}'.format(time_yy1, time_yy1+years-1), fontsize=20    )
        axs[1].set_title('20{}-20{}'.format(time_yy2, time_yy2+years-1), fontsize=20    )
        

    # final 
    plt.tight_layout()

    figname = 'fig/xcoc2_map_meshgrid_20{}_20{}_{}_{}years.png'.format(time_yy1, time_yy2, grid_size, years)
    plt.savefig(figname, dpi=250)
    print('{} is generated.'.format(figname))

    return







        