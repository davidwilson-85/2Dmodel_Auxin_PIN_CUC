#!/usr/bin/env python

import time
import os, shutil, random
import numpy as np
#import matplotlib.pyplot as plt
#from PIL import Image, ImageDraw, ImageFont

import params as pr
import inputs as ip

import func_graph
import func_auxin
import func_cuc
import func_pin


print('shape', ip.auxin.shape)
print('cols: ', ip.tissue_columns)
print('rows: ', ip.tissue_rows)

start_time = time.time()

# =====================================================================================

# NEW SIMULATION

# Cleanup destination folder (remove and create)
shutil.rmtree(pr.img_dest_folder) 
os.mkdir(pr.img_dest_folder)


for iteration in range(pr.nbr_iterations):

	# Print iteration to terminal
	if iteration < pr.nbr_iterations - 1:
		print(iteration + 1, end='\r')
	else:
		print(iteration + 1, end='\n')
	
	# Draw cell plot
	if iteration % pr.cell_plot_frequency == 0:
		func_graph.create_cell_plot(iteration)

	# Apply noise to auxin concentrations
	if pr.auxin_noise_factor > 0 and iteration == 0:
		func_auxin.auxin_noise()

	# PIN1 polarity
	func_pin.pin_polarity(pr.pin1_polarity)

	# PIN1 expression
	if pr.k_auxin_pin1 > 0 or pr.k_cuc_pin1 > 0 or pr.k_pin1_decay > 0:
		func_pin.pin_expression()

	# Auxin diffusion
	if pr.k_auxin_diffusion > 0:
		func_auxin.auxin_diffusion()

	# PIN1-mediated auxin transport
	if pr.k_pin1_transp > 0:
		func_auxin.pin_on_auxin(pr.k_pin1_transp)

	#if pr.k_auxin_synth > 0 or pr.k_cuc_yuc > 0 or pr.k_auxin_decay > 0:
		#func_auxin.auxin_homeostasis(pr.euler_h, auxin, cuc, pr.k_auxin_synth, pr.k_cuc_yuc, pr.k_auxin_decay)

	# Custom auxin modification
	#ip.auxin[5,6] += 1.5
	#ip.auxin[11,6] -= 1.5


print("%s seconds" % (time.time() - start_time))

quit()

# =====================================================================================

# Cleanup destination folder (remove and create)
shutil.rmtree(pr.img_dest_folder) 
os.mkdir(pr.img_dest_folder)

# Perform simulation cycles
for iteration in range(pr.nbr_iterations):
	# Print iteration to terminal
	if iteration < pr.nbr_iterations - 1:
		print(iteration + 1, end='\r')
	else:
		print(iteration + 1, end='\n')
	#*************************************************************************************
	# DRAW CELL PLOT
	if iteration % pr.cell_plot_frequency == 0:
		func_graph.create_cell_plot(iteration)
	#*************************************************************************************
	# Apply noise to auxin [THIS IS BEING MIGRATED TO AUXIN HOMEOSTASIS]
	if pr.auxin_noise_factor > 0:
		func_auxin.auxin_noise()
	#*************************************************************************************
	# AUXIN AND CUC EFFECT ON PIN1 EXPRESSION
	if pr.k_auxin_pin1 > 0 or pr.k_cuc_pin1 > 0 or pr.k_pin1_decay > 0:
		func_pin.pin_expression()
	#*************************************************************************************
	# CUC EXPRESSION
	if pr.k_md_cuc > 0 or pr.k_auxin_cuc > 0 or pr.k_cuc_decay > 0:
		func_cuc.cuc_expression()
	#*************************************************************************************
	# AUXIN HOMEOSTASIS
	#if pr.k_auxin_synth > 0 or pr.k_cuc_yuc > 0 or pr.k_auxin_degr > 0 or pr.auxin_noise_factor > 0:
	#	func_auxin.auxin_homeostasis(iteration)
	#*************************************************************************************
	# AUXIN DIFFUSION
	if pr.k_auxin_diffusion > 0:
		func_auxin.auxin_diffusion()
	#*************************************************************************************	
	# PIN1 EXPRESSION
	if pr.k_auxin_pin1 > 0 or pr.k_cuc_pin1 > 0 or pr.k_pin1_decay > 0:
		func_pin.pin_expression()
	#*************************************************************************************	
	# PIN1 POLARIZATION
	func_pin.pin_polarity(pr.pin1_polarity)
	#*************************************************************************************
	# PIN1-MEDIATED AUXIN EFFLUX
	if pr.k_pin1_transp > 0:
		func_auxin.pin_on_auxin(pr.k_pin1_transp
	#*************************************************************************************
	# Custom auxin modification [THIS IS GOING SOON TO AUXIN HOMEOSTASIS]
	#ip.auxin[5,6] += 1.5
	#ip.auxin[11,6] -= 1.5

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
