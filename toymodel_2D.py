#!/usr/bin/env python

import os, shutil, random
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

import params

import func_graph
import func_auxin
import func_cuc
import func_pin
import func_pinB


# Local synthesis or degradation (absolute or relative)
	# Here define a list of elements, each specifying the cell coordinates, synth/degr, abs/rel, cycles, etc

	# Example: marginSynthesis = new Local(type='synth', mode='abs', coords=(2,4,0), cycles=(0-4)) 


# === LOAD TEMPLATE DATA

# Templates
auxin_template = 'templates/auxin_template'
pin1_template = 'templates/pin1_template'
cuc_template = 'templates/cuc_template'
middle_domain_template = 'templates/middle_domain_template'

auxin = np.loadtxt(auxin_template, delimiter=',', unpack=False)
auxin = auxin * 10
auxin_matrix_shape = auxin.shape
tissue_rows, tissue_columns = auxin.shape[0], auxin.shape[1]
pin1 = np.loadtxt(pin1_template, delimiter=',', unpack=False).reshape((4,tissue_rows,tissue_columns)) # Format is [z,y,x]
pin1_matrix_shape = pin1.shape
cuc = np.loadtxt(cuc_template, delimiter=',', unpack=False)
middle_domain = np.loadtxt(middle_domain_template, delimiter=',', unpack=False)

auxin_fluxes = np.zeros(shape=(10,tissue_rows,tissue_columns)) # Z: T_out, T_in, R_out, R_in...
#array_auxin_net_fluxes = np.zeros(shape=(2,tissue_rows,tissue_columns)) # where z[0] => dx and z[1] => dy

# LUTs
lut_auxin = np.loadtxt('luts/lut_red_sat.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_pin1 = np.loadtxt('luts/lut_green.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_cuc = np.loadtxt('luts/lut_green_sat.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)

print('shape', auxin.shape)
print('cols: ', tissue_columns)
print('rows: ', tissue_rows)


# =====================================================================================

# Tests

for iteration in range(params.nbr_iterations):

	# Print iteration to terminal
	if iteration < params.nbr_iterations - 1:
		print(iteration + 1, end='\r')
	else:
		print(iteration + 1, end='\n')

	# PIN1 polarity
	func_pinB.pin_polarity(auxin, pin1, params.k_UTG, cuc, params.cuc_threshold_pin1)

	if params.k_auxin_diffusion > 0:

		func_auxin.auxin_diffusion(params.euler_h, params.k_auxin_diffusion, auxin_matrix_shape, tissue_columns, tissue_rows, auxin, auxin_fluxes, iteration)

	# Cleanup destination folder (remove and create)
	if iteration == 0:
		shutil.rmtree(params.img_dest_folder) 
		os.mkdir(params.img_dest_folder)
	
	# Draw cell plot 
	if iteration % params.cell_plot_frequency == 0:
		func_graph.create_cell_plot(
			auxin_matrix_shape,
			auxin,
			params.auxin_range,
			lut_auxin,
			pin1,
			params.pin1_range,
			lut_pin1,
			cuc,
			params.cuc_range,
			lut_cuc,
			iteration,
			auxin_fluxes,
			params.img_dest_folder
		)

	if params.k_pin1_transp > 0:

		func_auxin.pin_on_auxin(params.euler_h, auxin, pin1, params.k_pin1_transp, tissue_rows, tissue_columns, pin1_matrix_shape)


	if params.k_auxin_synth > 0 or params.k_cuc_yuc > 0 or params.k_auxin_decay > 0:
		
		func_auxin.auxin_homeostasis(params.euler_h, auxin, cuc, params.k_auxin_synth, params.k_cuc_yuc, params.k_auxin_decay)

	auxin[0,4] = auxin[0,4] + 2.5

quit()


# =====================================================================================


# === PROCESS DATA


# Perform simulation cycles

for iteration in range(params.nbr_iterations):
	
	# Print iteration to terminal
	if iteration < params.nbr_iterations - 1:
		print(iteration + 1, end='\r')
	else:
		print(iteration + 1, end='\n')
	
	#*************************************************************************************

	# Plot data in heatmap
	#create_heatmap(data=data)
	
	#*************************************************************************************

	# Plot data in cell tilling
	
	# Cleanup destination folder (remove and create)
	if iteration == 0:
		shutil.rmtree(params.img_dest_folder) 
		os.mkdir(params.img_dest_folder)
	
	# Draw cell plot 
	if iteration % params.cell_plot_frequency == 0:
		func_graph.create_cell_plot(
			auxin_matrix_shape,
			auxin,
			params.auxin_range,
			lut_auxin,
			pin1,
			params.pin1_range,
			lut_pin1,
			cuc,
			params.cuc_range,
			lut_cuc,
			iteration,
			auxin_fluxes,
			params.img_dest_folder
		)
	
	#*************************************************************************************
	
	# Apply noise to auxin
	
	if params.auxin_noise_factor > 0:
	
		for y in range(tissue_rows):
			for x in range(tissue_columns):
			
				auxin[y,x] = auxin[y,x] * random.uniform(-params.auxin_noise_factor, params.auxin_noise_factor)
			
				if auxin[y,x] < 0:
					auxin[y,x] = float(0.0000001)
	
	#*************************************************************************************
	
	# AUXIN AND CUC EFFECT ON PIN1 EXPRESSION
	
	if params.k_auxin_pin1 > 0 or params.k_cuc_pin1 > 0 or params.k_pin1_decay > 0:
		
		func_pin.pin_expression(pin1, auxin, cuc, params.k_auxin_pin1, params.k_cuc_pin1, params.k_pin1_decay)
	
	#*************************************************************************************
	
	# CUC EXPRESSION
	
	if params.k_md_cuc > 0 or params.k_auxin_cuc > 0 or params.k_cuc_decay > 0:
		
		func_cuc.cuc_expression(params.euler_h, middle_domain, auxin, cuc, params.k_md_cuc, params.k_auxin_cuc, params.k_cuc_decay)
	
	#*************************************************************************************

	# AUXIN HOMEOSTASIS
	
	if params.k_auxin_synth > 0 or params.k_cuc_yuc > 0 or params.k_auxin_decay > 0:
		
		func_auxin.auxin_homeostasis(params.euler_h, auxin, cuc, params.k_auxin_synth, params.k_cuc_yuc, params.k_auxin_decay)
	
	# Integrate local synthesis etc in the function...
	#auxin[0,0] = 9
	#auxin[4,4] = 0
	
	#*************************************************************************************

	# AUXIN DIFFUSION

	if params.k_auxin_diffusion > 0:

		func_auxin.auxin_diffusion(params.euler_h, params.k_auxin_diffusion, auxin_matrix_shape, tissue_columns, tissue_rows, auxin, auxin_fluxes, iteration)

	#*************************************************************************************	
	
	# PIN1 POLARIZATION

	if params.pin1_pol_mode == 'smith2006':
		
		# Smith 2006
		if params.k_UTG > 1:
			func_pin.pin_utg_smith2006(auxin, pin1, params.k_UTG, cuc, params.cuc_threshold_pin1)

	if params.pin1_pol_mode == 'ratio':
		
		# My method:
		if params.k_UTG > 0:
			func_pin.pin_utg_ratio(auxin, pin1, params.k_UTG, cuc, params.cuc_threshold_pin1)

	if params.pin1_pol_mode == 'wtf':

		func_pin.pin_wtf_p(auxin_fluxes, pin1, params.k_WTF)

	#*************************************************************************************

	# PIN1-MEDIATED AUXIN EFFLUX

	if params.k_pin1_transp > 0:

		func_auxin.pin_on_auxin(params.euler_h, auxin, pin1, params.k_pin1_transp, tissue_rows, tissue_columns, pin1_matrix_shape)
	
	auxin[3,3] = auxin[3,3] + 5
	auxin[10,3] = auxin[10,3] + 5
	auxin[0,] = 0
	auxin[13,] = 0
	auxin[:,0] = 0
	auxin[:,6] = 0

	




#===========================================
# TO DO:
# -Implement transport of auxin by PIN and see if instabilities in auxin distribution trigger patterning
# -Create quick way of switching ON/OFF features like UTG, WTF, CUC2
# 
# 
# matrix_shape[0] -> number of rows = y
# matrix_shape[1] -> number of columns = x
#
# auxin[6,11] = 5 refers to 7th row, 12th column
# pin1[0,0,5] = 5 refers to top wall of cell in 1st row, 6th column
# 
# 
# Set diffusion to 0 at faces that are outer boundaries
# 
# Check inner working of func_pin.auxin_on_pin_polarity()
# 
# Add pin1_UTGresponsiveness to function!!!!!!!!!!!!!!! It is not used at all at the moment!!!!!!!!!!!!!!!!
# 

