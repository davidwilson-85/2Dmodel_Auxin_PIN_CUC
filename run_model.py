#!/usr/bin/env python

import time, os, random, shutil, datetime
import numpy as np
from scipy.integrate import odeint

import importlib

import params as pr
#from sim_logs import params_2022_12_20_22_53_59 as pr # To re-run logged simulation
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
		series_num: When model is run as part of a series, this indicates the simulation number in the series.
	"""

	# When running as series, force reload the module inputs to reset all values to those in the parameters file. Then reassign the parameter being tested in the series.
	'''if pr.is_series == True:
		importlib.reload(pr)
		importlib.reload(ip)
		exec('pr.' + pr.series_param_a['name'] + '=' + str(param_val)) # exec() converts str to code'''

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
		#for y in range(ip.tissue_rows):
		#	for x in range(ip.tissue_columns):
		#		if ip.auxin[y,x] < 0: ip.auxin[y,x] = 0

		# Track simulation
		aux.track_simulation(iteration, nbr_iterations)

		# FOR TEMPORARY/TESTING FUNCTIONALY
		ip.cuc[4:7,4:7] = 9
		
	print("%s seconds" % (time.time() - start_time))
		
	# Graph auxin profile in central column of cells
	if pr.is_series == True:
		aux.create_line_plot_2(series_num, series_num_total)
	if pr.is_series == False:
		aux.create_line_plot(0, 1)

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

		# Run single simulation
		run()

	if pr.is_series == True:

		# Retrieve parameter that varies and range of values that it takes
		param_a_space = np.linspace(pr.series_param_a['min'], pr.series_param_a['max'], pr.series_param_a['num_points'])
		series_num_total = len(param_a_space)

		# Run simulation series
		for key, val in enumerate(param_a_space):
			#val = float(val)

			# Force reload the module inputs to reset all values to those in the parameters file (this is critical). Then reassign the parameter being tested in the series.
			importlib.reload(pr)
			importlib.reload(ip)
			exec('pr.' + pr.series_param_a['name'] + '=' + str(val)) # exec() converts str to code
			
			# Run
			print('series_num: ' + str(key) + '; ' + pr.series_param_a['name'] + ' = ' + str(val))
			run(series_num = key)

'''
TO DO:
* Try staggered cells
* Try PD growth?
* Try to integrate custom synth and degradation in ODEint model
* Check if the order in which the functions are called has an effect on output of simulations

'''
