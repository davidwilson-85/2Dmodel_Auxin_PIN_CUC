#!/usr/bin/env python

from os import system
import params as pr

pr.is_series = False

param_space = np.linspace(pr.series_param_a['min'], pr.series_param_a['max'], pr.series_param_a['num_points'])
# Calculate total number of simulations in the series and set parameter
series_num_total = len(param_space)

for key,val in enumerate(param_space):
    system(f"python3 run_model.py {key} {val}.")
    #read text file with auxin vals



system('python3 run_model.py 1 1')

#python3 func_graph.py datafile