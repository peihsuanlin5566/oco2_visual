import numpy as np 
import matplotlib as mpl
import util import fig_check_value_range, fig_point_plot, fig_map_meshgrid_2panels


if __name__ == "__main__": 

    # plotting settings 
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
    path_dir = './data/'
    
    # grid size (10: ~10 km per mesh grid)
    # 5 km/ per mesh grid by default
    grid_size = 5
    grid_num = int(200/5)

    # denoting the time 
    time1_start = '2016-01-01'
    time1_last = '2016-12-01'

    time2_start = '2017-01-01'
    time2_last = '2017-12-31'

    # plotting 
    # step 1: confirming that the specified range of values can present the data points correctly
    # here we use min=380, max=435   
    xco2_min = 380
    xco2_max = 435
    fig_check_value_range(path_dir, time1_start, time1_last,
                            xco2_min=xco2_min, xco2_max=xco2_max,)
    fig_check_value_range(path_dir, time1_start, time1_last,
                            xco2_min=xco2_min, xco2_max=xco2_max,) 

    # step 2: take a quick look of the plot: 
    fig_point_plot(path_dir, time1_start, time1_last,
                    xco2_min=xco2_min, xco2_max=xco2_max,)
    fig_point_plot(path_dir, time2_start, time2_last,
                    xco2_min=xco2_min, xco2_max=xco2_max,) 

    
    # step 3: plot the 2 panel mesh map overlaid w map
    fig_map_meshgrid_2panels(path_dir, time1_start, time1_last, time2_start, time2_last, grid_num=grid_num)

    # then you are done! 

