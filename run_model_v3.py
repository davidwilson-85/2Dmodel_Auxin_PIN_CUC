#!/usr/bin/env python

import time, os, random, shutil, datetime
import numpy as np
from scipy.integrate import odeint

import params_v3 as pr
import inputs_v3 as ip
import regulatory_network as rn
import func_graph_v3
import func_auxin_v3
import func_pin_v3

import auxiliary as aux
import tests.check as check

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

# =====================================================================================

# Cleanup destination folder (remove and create)
shutil.rmtree(pr.img_dest_folder) 
os.mkdir(pr.img_dest_folder)

# Perform simulation cycles
for iteration in range(nbr_iterations + 1):
	sim_time = iteration * pr.euler_h
	#print(sim_time)
	
	# Print iteration to terminal
	if iteration < nbr_iterations:
		print(str(iteration) + ' / ' + str(nbr_iterations), end='\r')
	else:
		print(str(iteration) + ' / ' + str(nbr_iterations), end='\n')
	
	# DRAW CELL PLOT
	if iteration % pr.cell_plot_frequency == 0:
		func_graph_v3.create_cell_plot(current_datetime, iteration)
	
	# SOLVE MODEL REGULATORY NETWORK
	rn.solve_model()
	
	# AUXIN CUSTOM MANIPULATION
	func_auxin_v3.auxin_custom_manipulation(iteration, sim_time)

	# Correct values out of bound (e.g. auxin < 0)
	######

	#*************************************************************************************
	# SOLVE REMAINING PROCESSES BY FORWARD EULER METHOD
	# AUXIN DIFFUSION
	if pr.k_auxin_diffusion > 0:
		func_auxin_v3.auxin_diffusion()
	# PIN1 POLARIZATION
	func_pin_v3.pin_polarity()
	# PIN1-MEDIATED AUXIN EFFLUX
	if pr.k_pin1_transp > 0:
		func_auxin_v3.pin_on_auxin()
	#*************************************************************************************

	# Track simulation
	aux.track_simulation(iteration, nbr_iterations)

	# FOR TEMPORARY/TESTING FUNCTIONALY

	#func_graph.create_heatmap(ip.auxin, iteration)

	#if iteration > 5:
		#ip.cuc[5:8,5:8] = 8
		#ip.auxin[8,6] += 5

	#if iteration == 1900:
		#print(ip.auxin)
		#print(np.array2string(ip.auxin, separator=','))
		#with open('templates/2D/template_auxin_1', 'w') as file:
			#aux.save_ndarray()
			#pass
	
print("%s seconds" % (time.time() - start_time))
# =====================================================================================

# Create video/gif files

if pr.create_video == True:
	func_graph_v3.create_video(current_datetime)

if pr.create_gif == True:
	func_graph_v3.create_gif(current_datetime)

print("%s seconds" % (time.time() - start_time))

'''
TO DO:
* Try to implement ODEint solving for each iteration. At least for the equations that describe expression of each component. 
* Check if the order in which the functions are called has an effect on output of simulations
* Create quick way of switching ON/OFF features like UTG, WTF, CUC2
* Set diffusion to 0 at faces that are outer boundaries
* Add pin1_UTGresponsiveness to function! It is not used at all at the moment!

'''



'''
FOR REFERENCE:

matrix_shape[0] -> number of rows = y
matrix_shape[1] -> number of columns = x

auxin[6,11] = 5 refers to 7th row, 12th column
pin1[0,0,5] = 5 refers to top wall of cell in 1st row, 6th column

'''
