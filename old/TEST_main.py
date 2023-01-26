#!/usr/bin/env python

import time, os, random, shutil, datetime
import numpy as np
from scipy.integrate import odeint

from sys import argv

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

def run(series_num = False, series_num_total = False):

	"""
	params:
		series_num: When model is run as part of a series, this indicates the simulation number in the series.
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

		'''if iteration  == nbr_iterations:
			func_graph.create_cell_plot(current_datetime, iteration, series_num = series_num)'''
		
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
		#ip.auxin[1,5] = 250
		
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


run()