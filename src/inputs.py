import os
from stl import mesh
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib import cm
import pandas as pd
from os.path import join
import open3d as o3d

work_dir = os.getcwd()
base_dir = join(work_dir, '../')

input_dir = join(base_dir, 'data/')
input_files = [f'S0{x}.ply' for x in range(1, 6)]

out_dir = join(base_dir, 'out/')
    
names = [x[:-4] for x in input_files]

colors = ['g', 'r', 'b', 'k'] 
rec_max_n = [1, 5, 10, 25]
max_n = 25    
