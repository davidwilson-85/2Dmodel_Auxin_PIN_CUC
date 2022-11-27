#!/usr/bin/env python

# Wrapper to create a series of simulatoins to automate parameter exploration
# Calls run_model.py several times sequentially, varying parameters each time
# Only one (or two) paramaters are variable; the remaining are fixed and are read from params.py
# Keeps the output (or part of it) of each run
# Makes plot(s) comparing all outputs
# 

import os, shutil

import numpy as np

import run_model
import params as pr

# ###
pr.is_series = True

# Cleanup destination folder (remove and create)
shutil.rmtree(pr.img_dest_folder)
os.mkdir(pr.img_dest_folder)

# Define parameter that varies and range of values that it takes
param_a_name = 'k_pin1_effi_basal'
param_a_values = {
	'min': 0,
	'max': .05,
	'num_points': 6
}
param_a_space = np.linspace(param_a_values['min'], param_a_values['max'], param_a_values['num_points'])
# Calculate total number of simulations in the series and set parameter
pr.series_num_total = len(param_a_space)

# Run series
for key, val in enumerate(param_a_space):
	exec('pr.' + param_a_name + '=' + str(val)) # exec() converts str to code
	print('series_num: ' + str(key))
	run_model.run(series_num = key)

