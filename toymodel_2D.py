#!/usr/bin/env python

import time, os, random, shutil, datetime

import numpy as np

import params as pr
import inputs as ip

import func_graph
import func_auxin
import func_cuc
import func_pin
import func_aux

# Calculate number of interations based on simulation time and step size
nbr_iterations = int(pr.simulation_time / pr.euler_h)

print('shape', ip.auxin.shape)

current_datetime = str(datetime.datetime.now())[:19].replace(':','-').replace(' ','_')

# Write initial state and simulation parameters to log
func_aux.write_to_log(current_datetime)

# Time execution of simulation
start_time = time.time()

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
	#*************************************************************************************
	# DRAW CELL PLOT
	if iteration % pr.cell_plot_frequency == 0:
		func_graph.create_cell_plot(current_datetime, iteration)
	#*************************************************************************************
	# PIN1 EXPRESSION (AUXIN AND CUC EFFECT)
	if pr.k_auxin_pin1 > 0 or pr.k_cuc_pin1 > 0 or pr.k_pin1_decay > 0:
		func_pin.pin_expression()
	#*************************************************************************************
	# CUC EXPRESSION
	if pr.k_cuc > 0:
		func_cuc.cuc_expression()
	#*************************************************************************************
	# AUXIN HOMEOSTASIS
	func_auxin.auxin_homeostasis(iteration, sim_time)
	#*************************************************************************************
	# AUXIN DIFFUSION
	if pr.k_auxin_diffusion > 0:
		func_auxin.auxin_diffusion()
	#*************************************************************************************	
	# PIN1 POLARIZATION
	#func_pin.pin_polarity(pr.pin1_polarity)
	#*************************************************************************************
	# PIN1-MEDIATED AUXIN EFFLUX
	#if pr.k_pin1_transp > 0:
	#	func_auxin.pin_on_auxin(pr.k_pin1_transp)
	#*************************************************************************************

	# FOR TEMPORARY/TESTING FUNCTIONALY

	#func_graph.create_heatmap(ip.auxin, iteration)

	#if iteration > 1000:
	#	ip.cuc[5:8,5:8] = 8
	#	#ip.auxin[8,6] += 5

	#if iteration == -1:
		#print(np.array2string(ip.auxin, separator=','))
		#with open('templates/2D/template_auxin_1', 'w') as file:
			#func_aux.save_ndarray()
			#pass
	
print("%s seconds" % (time.time() - start_time))
# =====================================================================================

# Create video/gif files

if pr.create_video == True:
	func_graph.create_video(current_datetime)

if pr.create_gif == True:
	func_graph.create_gif(current_datetime)

print("%s seconds" % (time.time() - start_time))

'''
TO DO:

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
