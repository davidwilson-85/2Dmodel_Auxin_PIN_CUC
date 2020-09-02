#!/usr/bin/env python

#import tilings
import os, shutil, math, random
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

import functions
import func_pin
import func_cuc


# === SET PARAMETERS

# Templates
auxin_template = 'templates/auxin_template'
pin1_template = 'templates/pin1_template'
cuc_template = 'templates/cuc_template'

# Switches
Diffusion = True
PIN1_UTG = True
PIN1_WTF = 'linear' # linear, cuadratic...
PIN1_transport = True
AUX_LAX_transport = False
CUC = False
Auxin_noise = False

# Parameters
auxin_range = (0, 9)
auxin_diffusionFactor = 0.2		# Relative amount of molecules that cross between two adjacent cells per cycle
auxin_synthesis = 0.1     		# Absolute amount of molecules synthesized per cycle
auxin_destruction = 0.1     	# Absolute amount of molecules destroyed per cycle
auxin_noise_factor = 0.1

pin1_range = (0, 9)
pin1_UTGresponsiveness = 0.001	# Relative amount of molecules that can change cell face per cycle
###################################
#
# Add pin1_UTGresponsiveness to function!!!!!!!!!!!!!!! It is not used at all at the moment!!!!!!!!!!!!!!!!
#
###################################
pin1_transp = 0.01				# = Nbr auxin molecules transported / ( PIN1 molecule * cycle )

cuc_range = (0, 9)
auxin_on_cuc = 0
cuc_on_pin1Pol = 0

nbr_iterations = 500
img_dest_folder = 'images/test'
cell_plot_frequency = 1

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


# === FUNCTIONS (To move to functions.py)

# Create a heatmap using matplotlib's imshow()
def create_heatmap(data):
	
	fig = plt.imshow(data, cmap="plasma", vmin=0, vmax=1, interpolation='none')
	plt.savefig('images/test/image' + str(iteration) + '.png')
	plt.close()


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
		functions.create_cell_plot(auxin_matrix_shape, auxin, auxin_range, lut_auxin, pin1, pin1_range, lut_pin1, cuc, cuc_range, lut_cuc, iteration, array_auxin_fluxes, img_dest_folder)
	
	#*************************************************************************************
	
	# Apply noise to auxin
	
	if Auxin_noise == True:
	
		for y in range(tissue_rows):
			for x in range(tissue_columns):
			
				auxin[x,y] = auxin[x,y] + random.uniform(-auxin_noise_factor, auxin_noise_factor)
			
				if auxin[x,y] < 0:
					auxin[x,y] = float(0.0000001)
	
	#*************************************************************************************

	# AUXIN DIFFUSION

	if Diffusion == True:

		functions.auxin_diffusion(auxin_diffusionFactor, auxin_matrix_shape, tissue_columns, tissue_rows, auxin, array_auxin_fluxes, iteration)

	#*************************************************************************************

	# AUXIN SYNTHESIS / DEGRADATION

	source = 7
	sink = 3

	#auxin[5,5] = auxin[5,5] + 0.5 


	#*************************************************************************************
	
	# PIN1 UTG

	# Calculate updated PIN1 allocation
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			# Top	
			if y > 0:
				utg_auxinRatioT = math.sqrt(auxin[y-1,x] / auxin[y,x])
			else:
				utg_auxinRatioT = 1.0
			
			# Right
			if x < tissue_columns - 1:
				utg_auxinRatioR = math.sqrt(auxin[y,x+1] / auxin[y,x])
			else:
				utg_auxinRatioR = 1.0
				
			# Bottom
			if y < tissue_rows - 1:
				utg_auxinRatioB = math.sqrt(auxin[y+1,x] / auxin[y,x])
			else:
				utg_auxinRatioB = 1.0
				
			# Left
			if x > 0:
				utg_auxinRatioL = math.sqrt(auxin[y,x-1] / auxin[y,x])
			else:
				utg_auxinRatioL = 1.0
			
			# Calculate updated normalized PIN1 concentration values per cell face
			total_pin1 = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]

			updated_pin1_T = pin1[0,y,x] * utg_auxinRatioT
			updated_pin1_R = pin1[1,y,x] * utg_auxinRatioR
			updated_pin1_B = pin1[2,y,x] * utg_auxinRatioB
			updated_pin1_L = pin1[3,y,x] * utg_auxinRatioL

			total_pin1_updated = updated_pin1_T + updated_pin1_R + updated_pin1_B + updated_pin1_L

			normalization_ratio = total_pin1 / total_pin1_updated

			# Normalize (force the sum of PIN1 in all faces in the cell to remain unchanged)
			pin1[0,y,x] = updated_pin1_T * normalization_ratio
			pin1[1,y,x] = updated_pin1_R * normalization_ratio
			pin1[2,y,x] = updated_pin1_B * normalization_ratio
			pin1[3,y,x] = updated_pin1_L * normalization_ratio

			#pin1[0,0,5] = 9

			#print pin1[0,y,x], pin1[1,y,x], pin1[2,y,x], pin1[3,y,x]

			#print utg_auxinRatioT, utg_auxinRatioR, utg_auxinRatioB, utg_auxinRatioL
	
	
	
	# PIN1-MEDIATED AUXIN TRANSPORT
	# Apply PIN1 transport to auxin concentration values
	#
	# PIN1_Tr = Nbr auxin molecules / ( PIN1 molecule * cycle )
	#
	#
	#

	# Create absolute efflux transport values for each cell
	transpVectors = np.zeros(pin1_matrix_shape, dtype=(float,1)) # 3D array = (cell_face, column, row)
	
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			auxin_molecules = auxin[y,x]
			total_pin1 = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
			transported_molecules_total = auxin_molecules * total_pin1 * pin1_transp

			# Manually limit transport if it exceeds the amount of auxin molecules
			# Later on: calculate excess ratio and then use to reduce the the value of the vector below
			if transported_molecules_total > auxin_molecules:
				transported_molecules_total = auxin_molecules
				print 'warning: transported molecules had to be manually adjusted in cell ' + str(y), str(x)

			# To top (y,x -> y-1,x)
			if y > 0:
				transpVectors[0,y,x] = auxin_molecules * pin1[0,y,x] * pin1_transp
			else:
				transpVectors[0,y,x] = 0

			# To right (y,x -> y,x+1)
			if x < tissue_columns - 1:
				transpVectors[1,y,x] = auxin_molecules * pin1[1,y,x] * pin1_transp
			else:
				transpVectors[1,y,x] = 0

			# To bottom (y,x -> y+1,x)
			if y < tissue_rows - 1:
				transpVectors[2,y,x] = auxin_molecules * pin1[2,y,x] * pin1_transp
			else:
				transpVectors[2,y,x] = 0

			# To left (y,x -> y,x-1)
			if x > 0:
				transpVectors[3,y,x] = auxin_molecules * pin1[3,y,x] * pin1_transp
			else:
				transpVectors[3,y,x] = 0

	for y in range(tissue_rows):
		for x in range(tissue_columns):

			# From top (y,x <- y-1,x)
			if y > 0:
				transpFromTop = transpVectors[2,y-1,x]
			else:
				transpFromTop = 0

			# From right (y,x <- y,x+1)
			if x < tissue_columns - 1:
				transpFromRight = transpVectors[3,y,x+1]
			else:
				transpFromRight = 0 

			# From bottom (y,x <- y+1,x)
			if y < tissue_rows - 1:
				transpFromBottom = transpVectors[0,y+1,x]
			else:
				transpFromBottom = 0

			# From left (y,x <- y,x-1)
			if x > 0:
				transpFromLeft = transpVectors[1,y,x-1]
			else:
				transpFromLeft = 0

			# Update the auxin concentration in each cell
			total_efflux = transpVectors[0,y,x] + transpVectors[1,y,x] + transpVectors[2,y,x] + transpVectors[3,y,x]
			auxin[y,x] = auxin[y,x] - total_efflux + transpFromTop + transpFromRight + transpFromBottom + transpFromLeft

	




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
# 

