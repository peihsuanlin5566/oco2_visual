import numpy as np 
import os
import argparse
import time
from util import *
from datetime import datetime

def get_args_parser(known=False):
    
    parser = argparse.ArgumentParser('OCO2 visualization tool')

    # indicating the time duration for plotting the CO2 observaton 
    parser.add_argument('--date_init', default="2015-01-01", type=str, 
            help='the stat data of the data (format: yyyy-mm-dd)')
    parser.add_argument('--date_lat', default="2015-02-01", type=str, 
            help='the last date of the data (format: yyyy-mm-dd)')

    # optional: mesh grid arguments 
    parser.add_argument('--data_path', help='path for placing the OCO2 data')
    parser.add_argument('--grid_size', default=5, type=float, help='grid size in km')
    parser.add_argument('--single_mesh_map', default=True, type=bool, help='Plot a meshgrid co2 map')
    
    # optional: FOV srguments
    parser.add_argument('--fov_lon_left', default=False, type=float, help='Indicating the FOV')
    parser.add_argument('--fov_lon_right', default=False, type=float, help='Indicating the FOV')
    parser.add_argument('--fov_lat_up', default=False, type=float, help='Indicating the FOV')
    parser.add_argument('--fov_lat_bottom', default=False, type=float, help='Indicating the FOV')
   
    # parser.add_argument('--load_checkpoint', 
    #     default=0, type=int, help='load the checkpoint file under ./model/expXX/checkpoint/. Specify the XX here')

    # others
    opt = parser.parse_known_args()[0] if known else parser.parse_args()
    # opt = parser.parse_known_args()[0]
    return opt

if __name__ == "__main__": 

    args = get_args_parser() 

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
    # path_dir = '../../data/data_xco2/'
    path_dir = args.data_path
    
    # grid size (10: ~10 km per mesh grid)
    grid_size = args.grid_size
    grid_num = int(200/grid_size)

    # time duration specification 
    date_str1 = args.time_init
    date_str2 = args.time_last

    date_obj1 = datetime.strptime(date_str1, '%Y-%m-%d')
    date_obj2 = datetime.strptime(date_str2, '%Y-%m-%d')

    # duration
    duration = args.duration
    

    # plotting 
    if args.single_mesh_map: 
        fig_map_meshgrid(path_dir, time_yy, grid_num,)
    else : 
        print('double mesh co2 map')
        # util.fig_map_meshgrid(path_dir, time_yy, grid_num,)

    