#!/usr/bin/env python

#import tilings
import math
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw


# === SET PARAMETERS

auxin_template = 'templates/auxin_template'
pin1_template = 'templates/pin1_template'
cuc_template = 'templates/...'

# Parameters
auxin_range = (0, 9)
auxin_diffusionFactor = 0.1		# Relative amount of molecules that cross between two adjacent cells per cycle
auxin_synthesis = 0.1     		# Absolute amount of molecules synthesized per cycle
auxin_destruction = 0.1     	# Absolute amount of molecules destroyed per cycle
auxin_lossRate = 0.75

pin1_range = (0, 9)
pin1_UTGresponsiveness = 0.01	# Relative amount of molecules that can change cell face per cycle
pin1_transp = 0.1				# = Nbr auxin molecules transported / ( PIN1 molecule * cycle )

cuc_range = (0, 1)

nbr_iterations = 500


# === LOAD DATA

auxin = np.loadtxt(auxin_template, delimiter=',', unpack=False)
auxin_matrix_shape = auxin.shape
pin1 = np.loadtxt(pin1_template, delimiter=',', unpack=False).reshape((4,11,11)) # Format is [z,y,x]
pin1_matrix_shape = pin1.shape

# LUTs
lut_auxin = np.loadtxt('luts/lut_red.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_pin1 = np.loadtxt('luts/lut_green.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_cuc = np.loadtxt('luts/lut_fire.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)


matrix_num_columns, matrix_num_rows = auxin.shape[0], auxin.shape[1]


# === FUNCTIONS

# Create a heatmap using matplotlib's imshow()
def create_heatmap(data):
	
	fig = plt.imshow(data, cmap="plasma", vmin=0, vmax=1, interpolation='none')
	plt.savefig('images/test/image' + str(iteration) + '.png')
	plt.close()


# Maps an integer representing the amount of a magnitude (e.g. [auxin]) and translates it to the corresponding RGB triplet of the selected LUT
def index_to_rgb(lut, level, range):

	# Clip values out of range
	if level < range[0]: level = range[0]
	if level > range[1]: level = range[1]

	# Rescale range to 0-255 (this is typical lut range)
	rescaled_level = int( ( level / range[1] ) * 255 )
	
	# Create RGB triplet
	rgb_triplet = (lut[1,rescaled_level],lut[2,rescaled_level],lut[3,rescaled_level])
	return rgb_triplet


# Create cell plot using pil
def create_cell_plot(matrix_shape):
	
	im = Image.new('RGBA', size=(700,700))
	x_origin = 0
	y_origin = 0
	cellSide = 50
	y = y_origin

	for i in range(matrix_shape[0]):
		x = x_origin
		
		for j in range(matrix_shape[1]):
			
			# Draw cell outline and fill
			ImageDraw.Draw(im).polygon([(x,y),(x+cellSide,y),(x+cellSide,y+cellSide),(x,y+cellSide)], outline=(50,50,50,255), fill=index_to_rgb(lut_auxin, auxin[i,j], auxin_range))
			
			# Draw PIN1 (top, right, bottom, left)
			#ImageDraw.Draw(im).line([(x+1,y+2),(x+cellSide-1,y+2)], fill=index_to_rgb(lut_pin1, pin1[0,i,j], pin1_range), width=3)
			#ImageDraw.Draw(im).line([(x+cellSide-2,y+2),(x+cellSide-2,y+cellSide-1)], fill=index_to_rgb(lut_pin1, pin1[1,i,j], pin1_range), width=3)
			#ImageDraw.Draw(im).line([(x+1,y+cellSide-2),(x+cellSide-1,y+cellSide-2)], fill=index_to_rgb(lut_pin1, pin1[2,i,j], pin1_range), width=3)
			#ImageDraw.Draw(im).line([(x+2,y+2),(x+2,y+cellSide-2)], fill=index_to_rgb(lut_pin1, pin1[3,i,j], pin1_range), width=3)

			ImageDraw.Draw(im).line([(x+4,y+3),(x+cellSide-4,y+3)], fill=index_to_rgb(lut_pin1, pin1[0,i,j], pin1_range), width=3)
			ImageDraw.Draw(im).line([(x+cellSide-3,y+4),(x+cellSide-3,y+cellSide-4)], fill=index_to_rgb(lut_pin1, pin1[1,i,j], pin1_range), width=3)
			ImageDraw.Draw(im).line([(x+4,y+cellSide-3),(x+cellSide-4,y+cellSide-3)], fill=index_to_rgb(lut_pin1, pin1[2,i,j], pin1_range), width=3)
			ImageDraw.Draw(im).line([(x+3,y+4),(x+3,y+cellSide-4)], fill=index_to_rgb(lut_pin1, pin1[3,i,j], pin1_range), width=3)
			
			# Draw CUC		
			#ImageDraw.Draw(im).ellipse([(x+15,y+15),(x+cellSide-15,y+cellSide-15)], fill='green')
			
			x = x + cellSide
		
		y = y + cellSide

	x = x_origin
	y = y_origin

	#ImageDraw.Draw(im).line([(x+1,y+2),(x+cellSide-1,y+2)], fill='red', width=3)

	im.save('images/test/image' + str(iteration) +'.png')

def auxin_diffusion(auxin_diffusionFactor):
	
	#
	# [auxin](i,j) = [auxin](i,j) - out_diff + in_diff_T + in_diff_R + in_diff_B + in_diff_L
	#
	# Rate of diffusion from cell i to j: dD(i->j)/dt = [auxin(i)] * k
	#
	# k = diffusion factor constant
	#

	# Create absolute efflux diffusion values for each cell
	diffusionVectors = np.zeros(auxin_matrix_shape, dtype=(float,1))	
	for y in range(matrix_num_columns):
		for x in range(matrix_num_rows):
			#print auxin[i,j]
			diffusionVectors[y,x] = auxin[y,x] * auxin_diffusionFactor
			
	#print diffusionVectors
	
	# Apply diffusion to auxin concentration values
	for y in range(matrix_num_columns):
		for x in range(matrix_num_rows):

			# Count how many neighbours the cell has, to calculate the amount of efflux diffusion
			nbr_cell_neighbours = 0
			
			# numpy arrays accept negative indexes, so testing if an index exists at the left and top of the grid, and at the right and bottom, must be done in a different way: 
			
			# Left
			if y-1 >= 0:
				diffusionFromLeft = diffusionVectors[y-1,x]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromLeft = 0
			
			# Right
			try: 
				diffusionFromRight = diffusionVectors[y+1,x]
				nbr_cell_neighbours	+=1
			except IndexError:
				diffusionFromRight = 0
			
			# Top	
			if x-1 >= 0:
				diffusionFromTop = diffusionVectors[y,x-1]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromTop = 0
			
			# Bottom
			try: 
				diffusionFromBottom = diffusionVectors[y,x+1]
				nbr_cell_neighbours	+=1
			except IndexError:
				diffusionFromBottom = 0
			
			# Update the concentration in each cell
			auxin[y,x] = auxin[y,x] - ( diffusionVectors[y,x] * nbr_cell_neighbours ) + diffusionFromLeft + diffusionFromRight + diffusionFromTop + diffusionFromBottom


	
# =====================================================================================
# =====================================================================================


# === PROCESS DATA

for iteration in range(nbr_iterations):

	# Plot data in heatmap
	#create_heatmap(data=data)

	# Plot data in cell tilling
	create_cell_plot(auxin_matrix_shape)

	
	# AUXIN DIFUSSION

	auxin_diffusion(auxin_diffusionFactor)


	# Apply PIN1 transport to auxin concentration values
	#
	# PIN1_Tr = Nbr auxin molecules / ( PIN1 molecule * cycle )
	#
	#
	#

	# Create absolute efflux transport values for each cell
	transpVectors = np.zeros(pin1_matrix_shape, dtype=(float,1)) # 3D array = (cell_face, column, row)
	
	for y in range(matrix_num_columns):
		for x in range(matrix_num_rows):
			#print auxin[i,j]

			total_molecules = auxin[y,x]
			total_pin1 = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
			transported_molecules = total_pin1 * pin1_transp

			ratio = total_molecules / transported_molecules
			if ratio > 1: ratio = 1


			#transpVectors[0,j,i] = auxin[y,x] * auxin_diffusionFactor
	
	for y in range(pin1_matrix_shape[0]):
		for x in range(pin1_matrix_shape[1]):

			# Top
			transpFromTop = 0


			# Update the auxin concentration in each cell
			#auxin[i,j] = auxin[i,j] - ( transpVectors[i,j] * nbr_cell_neighbours ) + transpFromLeft + transpFromRight + transpFromTop + transpFromBottom

	



	# Calculate updated PIN1 allocation
	for y in range(auxin_matrix_shape[0]):
		for x in range(auxin_matrix_shape[1]):

			# Left
			if x-1 >= 0:
				utg_auxinRatioL = math.sqrt(auxin[y,x-1] / auxin[y,x])
			else:
				utg_auxinRatioL = 1.0
			
			# Right
			try:
				utg_auxinRatioR = math.sqrt(auxin[y,x+1] / auxin[y,x])
			except IndexError:
				utg_auxinRatioR = 1.0
			
			# Top	
			if y-1 >= 0:
				utg_auxinRatioT = math.sqrt(auxin[y-1,x] / auxin[y,x])
			else:
				utg_auxinRatioT = 1.0
			
			# Bottom
			try:
				utg_auxinRatioB = math.sqrt(auxin[y+1,x] / auxin[y,x])
			except IndexError:
				utg_auxinRatioB = 1.0

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



	
	auxin[3,3] = auxin[3,3] + auxin_synthesis
	auxin[10,5] = 0.5
	#auxin[9,11] = auxin[9,11] + auxin_synthesis



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
# 
# 
# 


