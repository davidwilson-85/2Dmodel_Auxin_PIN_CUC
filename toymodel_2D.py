#!/usr/bin/env python

#import tilings
import os, shutil, random
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

import func_graph
import func_auxin
import func_cuc
import func_pin


# === SET PARAMETERS

# Templates
auxin_template = 'templates/auxin_template'
pin1_template = 'templates/pin1_template'
cuc_template = 'templates/cuc_template'

# Switches
PIN1_UTG = True
PIN1_WTF = 'linear' # linear, cuadratic...
AUX_LAX_transport = False
CUC = False

# Parameters
auxin_range = (0, 9)            # This is only to map variable values to heatmap values
k_auxin_diffusion = 0.2			# Relative amount of molecules that cross between two adjacent cells per cycle
auxin_synthesis = 0.1     		# Absolute amount of molecules synthesized per cycle
auxin_destruction = 0.1     	# Absolute amount of molecules destroyed per cycle
auxin_noise_factor = 0

pin1_range = (0, 9)
k_pin1_UTGresponsiveness = 0.01	# Relative amount of molecules that can change cell face per cycle
###################################
#
# Add pin1_UTGresponsiveness to function!!!!!!!!!!!!!!! It is not used at all at the moment!!!!!!!!!!!!!!!!
#
###################################
k_pin1_transp = 0.01			# = Nbr auxin molecules transported / ( PIN1 molecule * cycle )
cuc_range = (0, 9)				# This is only to map variable values to heatmap values
auxin_on_cuc = 0
cuc_on_pin1Pol = 0

nbr_iterations = 500
img_dest_folder = 'images/test'
cell_plot_frequency = 499

# Local synthesis or degradation (absolute or relative)
	# Here define a list of elements, each specifying the cell coordinates, synth/degr, abs/rel, cycles, etc

	# Example: marginSynthesis = new Local(type='synth', mode='abs', coords=(2,4,0), cycles=(0-4)) 


# === LOAD DATA

auxin = np.loadtxt(auxin_template, delimiter=',', unpack=False)
auxin_matrix_shape = auxin.shape
tissue_rows, tissue_columns = auxin.shape[0], auxin.shape[1]
pin1 = np.loadtxt(pin1_template, delimiter=',', unpack=False).reshape((4,tissue_rows,tissue_columns)) # Format is [z,y,x]
pin1_matrix_shape = pin1.shape

cuc = np.loadtxt(cuc_template, delimiter=',', unpack=False)

array_auxin_fluxes = np.zeros(shape=(10,tissue_rows,tissue_columns)) # Z: T_out, T_in, R_out, R_in...
#array_auxin_net_fluxes = np.zeros(shape=(2,tissue_rows,tissue_columns)) # where z[0] => dx and z[1] => dy

# LUTs
lut_auxin = np.loadtxt('luts/lut_red.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_pin1 = np.loadtxt('luts/lut_green.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_cuc = np.loadtxt('luts/lut_fire.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)

print 'shape', auxin.shape
print 'cols: ', tissue_columns
print 'rows: ', tissue_rows


# =====================================================================================
# =====================================================================================


# === PROCESS DATA


# Perform simulation cycles

for iteration in range(nbr_iterations):

	print iteration

	# Plot data in heatmap
	#create_heatmap(data=data)
	
	#*************************************************************************************

	# Plot data in cell tilling
	
	# Cleanup destination folder (remove and create)
	if iteration == 0:
		shutil.rmtree(img_dest_folder) 
		os.mkdir(img_dest_folder)
	
	# Draw cell plot 
	if iteration % cell_plot_frequency == 0:
		func_graph.create_cell_plot(auxin_matrix_shape, auxin, auxin_range, lut_auxin, pin1, pin1_range, lut_pin1, cuc, cuc_range, lut_cuc, iteration, array_auxin_fluxes, img_dest_folder)
	
	#*************************************************************************************
	
	# Apply noise to auxin
	
	if auxin_noise_factor > 0:
	
		for y in range(tissue_rows):
			for x in range(tissue_columns):
			
				auxin[x,y] = auxin[x,y] + random.uniform(-auxin_noise_factor, auxin_noise_factor)
			
				if auxin[x,y] < 0:
					auxin[x,y] = float(0.0000001)
	
	#*************************************************************************************

	# AUXIN DIFFUSION

	if k_auxin_diffusion > 0:

		func_auxin.auxin_diffusion(k_auxin_diffusion, auxin_matrix_shape, tissue_columns, tissue_rows, auxin, array_auxin_fluxes, iteration)

	#*************************************************************************************

	# AUXIN SYNTHESIS / DEGRADATION

	source = 7
	sink = 3

	#auxin[5,5] = auxin[5,5] + 0.5 

	#*************************************************************************************	
	
	# PIN1 UTG

	if k_pin1_UTGresponsiveness > 0:

		func_pin.auxin_on_pin_polarity(auxin, pin1, k_pin1_UTGresponsiveness, tissue_rows, tissue_columns)

	#*************************************************************************************

	# PIN1-MEDIATED AUXIN EFFLUX

	if k_pin1_transp > 0:

		func_auxin.pin_on_auxin(auxin, pin1, k_pin1_transp, tissue_rows, tissue_columns, pin1_matrix_shape)
	
	

	




#===========================================
# TO DO:
# -Implement transport of auxin by PIN and see if instabilities in auxin distribution trigger patterning
# -Create quick way of swutching ON/OFF features like UTG, WTF, CUC2
# -Add auxin effect on PIN1 expression
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

