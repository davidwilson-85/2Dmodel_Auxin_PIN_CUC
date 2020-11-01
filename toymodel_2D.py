#!/usr/bin/env python

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


# =====================================================================================

# Tests

for iteration in range(pr.nbr_iterations):

	# Print iteration to terminal
	if iteration < pr.nbr_iterations - 1:
		print(iteration + 1, end='\r')
	else:
		print(iteration + 1, end='\n')

	# PIN1 polarity
	func_pin.pin_polarity()


	if pr.k_auxin_pin1 > 0 or pr.k_cuc_pin1 > 0 or pr.k_pin1_decay > 0:
		func_pin.pin_expression()


	if pr.k_auxin_diffusion > 0:

		func_auxin.auxin_diffusion(pr.euler_h, pr.k_auxin_diffusion, ip.auxin_matrix_shape, ip.tissue_columns, ip.tissue_rows, ip.auxin, ip.auxin_fluxes, iteration)

	# Cleanup destination folder (remove and create)
	if iteration == 0:
		shutil.rmtree(pr.img_dest_folder) 
		os.mkdir(pr.img_dest_folder)
	
	# Draw cell plot 
	if iteration % pr.cell_plot_frequency == 0:
		func_graph.create_cell_plot(
			ip.auxin_matrix_shape,
			ip.auxin,
			pr.auxin_range,
			ip.lut_auxin,
			ip.pin1,
			pr.pin1_range,
			ip.lut_pin1,
			ip.cuc,
			pr.cuc_range,
			ip.lut_cuc,
			iteration,
			ip.auxin_fluxes,
			pr.img_dest_folder
		)

	#if pr.k_pin1_transp > 0:

		#func_auxin.pin_on_auxin(pr.euler_h, auxin, pin1, pr.k_pin1_transp, tissue_rows, tissue_columns, pin1_matrix_shape)


	#if pr.k_auxin_synth > 0 or pr.k_cuc_yuc > 0 or pr.k_auxin_decay > 0:
		
		#func_auxin.auxin_homeostasis(pr.euler_h, auxin, cuc, pr.k_auxin_synth, pr.k_cuc_yuc, pr.k_auxin_decay)

	ip.auxin[0,4] = ip.auxin[0,4] + 2.5

quit()


# =====================================================================================


# === PROCESS DATA


# Perform simulation cycles

for iteration in range(pr.nbr_iterations):
	
	# Print iteration to terminal
	if iteration < pr.nbr_iterations - 1:
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
		shutil.rmtree(pr.img_dest_folder) 
		os.mkdir(pr.img_dest_folder)
	
	# Draw cell plot 
	if iteration % pr.cell_plot_frequency == 0:
		func_graph.create_cell_plot(
			auxin_matrix_shape,
			auxin,
			pr.auxin_range,
			lut_auxin,
			pin1,
			pr.pin1_range,
			lut_pin1,
			cuc,
			pr.cuc_range,
			lut_cuc,
			iteration,
			auxin_fluxes,
			pr.img_dest_folder
		)
	
	#*************************************************************************************
	
	# Apply noise to auxin
	
	if pr.auxin_noise_factor > 0:
	
		for y in range(tissue_rows):
			for x in range(tissue_columns):
			
				auxin[y,x] = auxin[y,x] * random.uniform(-pr.auxin_noise_factor, pr.auxin_noise_factor)
			
				if auxin[y,x] < 0:
					auxin[y,x] = float(0.0000001)
	
	#*************************************************************************************
	
	# AUXIN AND CUC EFFECT ON PIN1 EXPRESSION
	
	if pr.k_auxin_pin1 > 0 or pr.k_cuc_pin1 > 0 or pr.k_pin1_decay > 0:
		
		#func_pin.pin_expression(pin1, auxin, cuc, pr.k_auxin_pin1, pr.k_cuc_pin1, pr.k_pin1_decay)
		func_pin.pin_expression()

	#*************************************************************************************
	
	# CUC EXPRESSION
	
	if pr.k_md_cuc > 0 or pr.k_auxin_cuc > 0 or pr.k_cuc_decay > 0:
		
		func_cuc.cuc_expression(pr.euler_h, middle_domain, auxin, cuc, pr.k_md_cuc, pr.k_auxin_cuc, pr.k_cuc_decay)
	
	#*************************************************************************************

	# AUXIN HOMEOSTASIS
	
	if pr.k_auxin_synth > 0 or pr.k_cuc_yuc > 0 or pr.k_auxin_decay > 0:
		
		func_auxin.auxin_homeostasis(pr.euler_h, auxin, cuc, pr.k_auxin_synth, pr.k_cuc_yuc, pr.k_auxin_decay)
	
	# Integrate local synthesis etc in the function...
	#auxin[0,0] = 9
	#auxin[4,4] = 0
	
	#*************************************************************************************

	# AUXIN DIFFUSION

	if pr.k_auxin_diffusion > 0:

		func_auxin.auxin_diffusion(pr.euler_h, pr.k_auxin_diffusion, auxin_matrix_shape, tissue_columns, tissue_rows, auxin, auxin_fluxes, iteration)

	#*************************************************************************************	
	
	# PIN1 POLARIZATION

	if pr.pin1_pol_mode == 'smith2006':
		
		# Smith 2006
		if pr.k_UTG > 1:
			func_pin.pin_utg_smith2006(auxin, pin1, pr.k_UTG, cuc, pr.cuc_threshold_pin1)

	if pr.pin1_pol_mode == 'ratio':
		
		# My method:
		if pr.k_UTG > 0:
			func_pin.pin_utg_ratio(auxin, pin1, pr.k_UTG, cuc, pr.cuc_threshold_pin1)

	if pr.pin1_pol_mode == 'wtf':

		func_pin.pin_wtf_p(auxin_fluxes, pin1, pr.k_WTF)

	#*************************************************************************************

	# PIN1-MEDIATED AUXIN EFFLUX

	if pr.k_pin1_transp > 0:

		func_auxin.pin_on_auxin(pr.euler_h, auxin, pin1, pr.k_pin1_transp, tissue_rows, tissue_columns, pin1_matrix_shape)
	
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

