#!/usr/bin/env python

# Wrapper to automate parameter exploration
# Calls run_model.py several times sequentially, varying parameters each time
# Only one (or two) paramaters are variable; the remaining are fixed and are read from params.py
# Keeps the output (or part of it) of each run
# Makes plot(s) comparing all outputs
# 

import numpy as np

import run_model
import params as pr

# Override to suppress creation of imags/vids
pr.cell_plot_frequency = False
pr.create_video = False 

# Define parameter that varies and values that it takes
param_a_name = 'k_cuc_auxin_synth'
param_a_values = {
	'min': 0,
	'max': 1,
	'num_points': 11
}
param_a_space = np.linspace(param_a_values['min'], param_a_values['max'], param_a_values['num_points'])

for i in param_a_space:
	exec('pr.' + param_a_name + '=' + str(i)) # exec() converts str to code
	run_model.run()

