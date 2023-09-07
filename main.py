#!/usr/bin/env python

import time, datetime, importlib, sys
import numpy as np
from scipy.integrate import odeint

try:
	pr = importlib.import_module(sys.argv[1].split('.')[0], package=None)
except IndexError:
	print('Argument missing. Please specify a parameters file as the first argument.')
	quit()

#import params as pr
import inputs as ip
import regulatory_network as rn
import func_graph
import func_auxin
import func_pin
import func_cuc

import auxiliary as aux


def run(series_num_total = False, series_num = False):

	"""
	Simulation loop function. Called by run_model.py 
	
	Params:
		series_num_total: When model is run as part of a series, this indicates the total number of simulations in the series.
		series_num: When model is run as part of a series, this indicates the simulation number in the series.
	
	Input: 
		All input is stored in file params.py, imported here
	
	Output:
		A series of files (graphs, images, videos, text files)
	"""

	# Calculate number of iterations based on simulation time and step size
	nbr_iterations = int(pr.simulation_time_total / pr.euler_h)

	timestamp = str(datetime.datetime.now())[:19].replace(':','-').replace(' ','_')

	# Write params.py (initial state and simulation parameters) to log
	aux.write_to_log(timestamp)

	# Time execution of simulation
	start_time = time.time()

	# Perform simulation cycles
	for iteration in range(nbr_iterations + 1):
		
		simulation_time_current = iteration * pr.euler_h
		
		# Print iteration to terminal
		if iteration < nbr_iterations:
			print(str(iteration) + ' / ' + str(nbr_iterations), end='\r')
		else:
			print(str(iteration) + ' / ' + str(nbr_iterations), end='\n')
		
		# OPTIONAL: DRAW CELL PLOT
		if pr.is_series == False:
			if (iteration * pr.euler_h) % pr.cell_plot_frequency == 0:
				func_graph.create_cell_plot(timestamp, iteration)
		if pr.is_series == True:
			if iteration == nbr_iterations:
				func_graph.create_cell_plot(timestamp, iteration, series_num = series_num)
		
		# SOLVE MODEL REGULATORY NETWORK
		rn.solve_rn_model()

		# SOLVE REMAINING PROCESSES
		func_auxin.auxin_custom_manipulation(iteration, simulation_time_current)
		func_auxin.auxin_diffusion()
		func_pin.pin_polarity()
		func_auxin.pin_on_auxin()

		if pr.cuc_noise['limit'] > 0: func_cuc.cuc_custom_manipulation(iteration)
		
		# Limit / correct values out of bound (e.g. auxin < 0)?
		#for y in range(ip.tissue_rows):
		#	for x in range(ip.tissue_columns):
		#		if ip.auxin[y,x] < 0: ip.auxin[y,x] = 0

		# Track simulation
		aux.track_simulation(iteration, nbr_iterations)
		#aux.make_kymograph(iteration, nbr_iterations)

		# FOR TEMPORARY/TESTING FUNCTIONALY
		#if simulation_time_current < 30: ip.cuc[:] = 0
		
	print("%s seconds" % (time.time() - start_time))
		
	# Graph auxin profile in central column of cells
	if pr.is_series == False:
		aux.create_line_plot_single(timestamp, 0, 1)
	if pr.is_series == True:
		aux.create_line_plot_multi(timestamp, series_num, series_num_total)

	# Create video/gif files
	if pr.is_series == False:
		if pr.create_video == True:
			func_graph.create_video(timestamp)
		if pr.create_gif == True:
			func_graph.create_gif(timestamp)
