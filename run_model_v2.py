#!/usr/bin/env python

import time, os, random, shutil, datetime

import numpy as np

import params_v2 as pr
import inputs_v2 as ip

import func_graph_v2
import func_auxin_v2
import func_cuc
import func_pin

import auxiliary as aux
import tests.check as check

import integrator as itg
from scipy.integrate import odeint

# Setup checks
check.check_dirs()

# Calculate number of interations based on simulation time and step size
nbr_iterations = int(pr.simulation_time / pr.euler_h)

#print('shape', ip.auxin.shape)

current_datetime = str(datetime.datetime.now())[:19].replace(':','-').replace(' ','_')

# Write initial state and simulation parameters to log
aux.write_to_log(current_datetime)

# Time execution of simulation
start_time = time.time()

time_points = np.linspace(0, 1, 2)

# =====================================================================================

# Cleanup destination folder (remove and create)
shutil.rmtree(pr.img_dest_folder) 
os.mkdir(pr.img_dest_folder)

# Perform simulation cycles
for iteration in range(nbr_iterations + 1):
	sim_time = iteration * pr.euler_h
	
	# Print iteration to terminal
	if iteration < nbr_iterations:
		print(str(iteration) + ' / ' + str(nbr_iterations), end='\r')
	else:
		print(str(iteration) + ' / ' + str(nbr_iterations), end='\n')
	
	# DRAW CELL PLOT
	if iteration % pr.cell_plot_frequency == 0:
		func_graph_v2.create_cell_plot(current_datetime, iteration)
	
	# SOLVE STEP

	# Update ip.auxin_neighbours
	func_auxin_v2.update_auxin_neighbours()

	# Compute cell by cell
	for y in range(ip.tissue_rows):
		for x in range(ip.tissue_columns):
			
			model_init_values = [
				ip.auxin[y,x],
				ip.auxin_neighbours[0,y,x],
				ip.auxin_neighbours[1,y,x],
				ip.auxin_neighbours[2,y,x],
				ip.auxin_neighbours[3,y,x],
				ip.cuc[x,y]
			]
			cell_solution = odeint(itg.model, model_init_values, time_points)
			#print('cell_solution: ')
			#print(cell_solution[-1,0])
			ip.auxin_tmp[y,x] = cell_solution[-1,0]
	
	# Consolidate auxin changes
	ip.auxin = np.copy(ip.auxin_tmp)

	print('--')
	print(ip.auxin)
	print('--')

print("%s seconds" % (time.time() - start_time))
# =====================================================================================

# Create video/gif files

if pr.create_video == True:
	func_graph_v2.create_video(current_datetime)

if pr.create_gif == True:
	func_graph_v2.create_gif(current_datetime)

print("%s seconds" % (time.time() - start_time))
