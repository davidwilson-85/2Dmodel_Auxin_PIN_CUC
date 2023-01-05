#!/usr/bin/env python

import time, os, random, shutil, datetime
import numpy as np
from scipy.integrate import odeint

import params as pr
#from sim_logs import params_2022_12_01_21_23_19 as pr # To re-run logged simulation
import inputs as ip
import regulatory_network as rn
import func_graph
import func_auxin
import func_pin

import auxiliary as aux
import tests.check as check

# Do checks
check.check_dirs()

def run(series_num = False):

	"""
	params:
		series_num: When model is run through run_model_series.py wrapper, this indicates the simulation number in the series.
	"""

	# Calculate number of iterations based on simulation time and step size
	nbr_iterations = int(pr.simulation_time / pr.euler_h)

	current_datetime = str(datetime.datetime.now())[:19].replace(':','-').replace(' ','_')

	# Write initial state and simulation parameters to log
	aux.write_to_log(current_datetime)

	# Time execution of simulation
	start_time = time.time()

	# Perform simulation cycles
	for iteration in range(nbr_iterations + 1):
		
		sim_time = iteration * pr.euler_h
		
		# Print iteration to terminal
		if iteration < nbr_iterations:
			print(str(iteration) + ' / ' + str(nbr_iterations), end='\r')
		else:
			print(str(iteration) + ' / ' + str(nbr_iterations), end='\n')
		
		# OPTIONAL: DRAW CELL PLOT
		if pr.is_series == False:
			if (iteration * pr.euler_h) % pr.cell_plot_frequency == 0:
				func_graph.create_cell_plot(current_datetime, iteration)
		if pr.is_series == True:
			if iteration  == nbr_iterations:
				func_graph.create_cell_plot(current_datetime, iteration, series_num = series_num)
		
		# SOLVE MODEL REGULATORY NETWORK
		rn.solve_model()

		# SOLVE REMAINING PROCESSES BY FORWARD EULER METHOD
		func_auxin.auxin_custom_manipulation(iteration, sim_time)
		if pr.k_auxin_diffusion > 0:
			func_auxin.auxin_diffusion()
		func_pin.pin_polarity()
		func_auxin.pin_on_auxin()
		
		# Limit / correct values out of bound (e.g. auxin < 0)?
		for y in range(ip.tissue_rows):
			for x in range(ip.tissue_columns):
				if ip.auxin[y,x] > 250: ip.auxin[y,x] = 250

		# Track simulation
		aux.track_simulation(iteration, nbr_iterations)

		# FOR TEMPORARY/TESTING FUNCTIONALY
		'''
		ip.auxin[1,5] = 250
		
		if sim_time >= 40:
			ip.cuc[4:7,4:7] = 8
			ip.auxin[4:7,4:7] += (1 * pr.euler_h)
		if sim_time >= 60:
			ip.auxin[:,5] += (4 * pr.euler_h)
		'''
		
	print("%s seconds" % (time.time() - start_time))
		
	# Track series
	if pr.is_series == True:
		aux.track_series(series_num, series_num_total)

	# Create video/gif files
	if pr.is_series == False:
		if pr.create_video == True:
			func_graph.create_video(current_datetime)
		if pr.create_gif == True:
			func_graph.create_gif(current_datetime)


if __name__ == '__main__':

	# Cleanup destination folder (remove and create)
	shutil.rmtree(pr.img_dest_folder) 
	os.mkdir(pr.img_dest_folder)
	
	if pr.is_series == False:

		# Run simulation
		run()

	if pr.is_series == True:

		# Retrieve parameter that varies and range of values that it takes
		param_a_space = np.linspace(pr.series_param_a['min'], pr.series_param_a['max'], pr.series_param_a['num_points'])
		# Calculate total number of simulations in the series and set parameter
		series_num_total = len(param_a_space)

		for i in param_a_space: print(i)

		# Run series
		for key, val in enumerate(param_a_space):
			exec('pr.' + pr.series_param_a['name'] + '=' + str(val)) # exec() converts str to code
			print('series_num: ' + str(key) + '; ' + pr.series_param_a['name'] + ' = ' + str(val))
			run(series_num = key)
		

'''
TO DO:
* Try staggered cells
* Try PD growth?
* Try to integrate custom synth and degradation in ODEint model
* Check if the order in which the functions are called has an effect on output of simulations

'''
