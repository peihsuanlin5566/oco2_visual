import numpy as np 
import matplotlib as mpl
import util 


if __name__ == "__main__": 

    params = {'mathtext.default': 'regular',
          'xtick.direction': 'in',
          'ytick.direction': 'in',
          'xtick.top': True,
          'ytick.right': True,
          'axes.labelsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'axes.titlesize':12 }
    mpl.rcParams.update(params)

    font = {'size'   : 15}
    mpl.rc('font', **font)

    # place where the data are placed
    path_dir = '../../data/data_xco2/oco2_LtCO2_*.nc4'
    
    # grid size (10: ~10 km per mesh grid)
    grid_size = 5
    grid_num = int(200/5)

    # year
    time_yys = 15

    # plotting 
    for time_yy in time_yys: 
        util.fig_map_meshgrid(path_dir, time_yy, grid_num,)

