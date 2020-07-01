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
auxin_synthesis = 1     		# Absolute amount of molecules synthesized per cycle
auxin_lossRate = 0.75

pin1_range = (0, 9)
pin1_UTGresponsiveness = 0.01	# Relative amount of molecules that can change cell face per cycle

cuc_range = (0, 1)

nbr_iterations = 500


# === LOAD DATA

auxin = np.loadtxt(auxin_template, delimiter=',', unpack=False)
matrix_shape = auxin.shape
pin1 = np.loadtxt(pin1_template, delimiter=',', unpack=False).reshape((4,12,12)) # Format is [z,y,x]

# LUTs
lut_auxin = np.loadtxt('luts/lut_red.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_pin1 = np.loadtxt('luts/lut_green.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_cuc = np.loadtxt('luts/lut_fire.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)


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



# === PROCESS DATA

for iteration in range(nbr_iterations):

	# Plot data in heatmap
	#create_heatmap(data=data)

	# Plot data in cell tilling
	create_cell_plot(matrix_shape)
	
	# Create absolute diffusion efflux values for each cell
	diffusionVectors = np.zeros(matrix_shape, dtype=(float,1))	
	for i in range(matrix_shape[0]):
		for j in range(matrix_shape[1]):
			#print auxin[i,j]
			diffusionVectors[i,j] = auxin[i,j] * auxin_diffusionFactor
			
	#print diffusionVectors
	
	# Calculate updated auxin concentration values
	for i in range(matrix_shape[0]):
		for j in range(matrix_shape[1]):

			# Count how many neighbours the cell has, to calculate the amount of efflux diffusion
			nbr_cell_neighbours = 0
			
			# numpy arrays accept negative indexes, so testing if an index exists at the left and top of the grid, and at the right and bottom, must be done in a different way: 
			
			# Left
			if i-1 >= 0:
				diffusionFromLeft = diffusionVectors[i-1,j]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromLeft = 0
			
			# Right
			try: 
				diffusionFromRight = diffusionVectors[i+1,j]
				nbr_cell_neighbours	+=1
			except IndexError:
				diffusionFromRight = 0
			
			# Top	
			if j-1 >= 0:
				diffusionFromTop = diffusionVectors[i,j-1]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromTop = 0
			
			# Bottom
			try: 
				diffusionFromBottom = diffusionVectors[i,j+1]
				nbr_cell_neighbours	+=1
			except IndexError:
				diffusionFromBottom = 0
			
			# Update the concentration in each cell
			auxin[i,j] = auxin[i,j] - ( diffusionVectors[i,j] * nbr_cell_neighbours ) + diffusionFromLeft + diffusionFromRight + diffusionFromTop + diffusionFromBottom


	# Calculate updated PIN1 allocation
	for i in range(matrix_shape[0]):
		for j in range(matrix_shape[1]):

			# Left
			if i-1 >= 0:
				utg_auxinRatioL = math.sqrt(auxin[i-1,j] / auxin[i,j])
			else:
				utg_auxinRatioL = 1.0
			
			# Right
			try:
				utg_auxinRatioR = math.sqrt(auxin[i+1,j] / auxin[i,j])
			except IndexError:
				utg_auxinRatioR = 1.0
			
			# Top	
			if j-1 >= 0:
				utg_auxinRatioT = math.sqrt(auxin[i,j-1] / auxin[i,j])
			else:
				utg_auxinRatioT = 1.0
			
			# Bottom
			try:
				utg_auxinRatioB = math.sqrt(auxin[i,j+1] / auxin[i,j])
			except IndexError:
				utg_auxinRatioB = 1.0

			# Calculate updated PIN1 concentration values per cell face
			total_pin1 = pin1[0,i,j] + pin1[1,i,j] + pin1[2,i,j] + pin1[3,i,j]

			updated_pin1_T = pin1[0,j,i] * utg_auxinRatioT
			updated_pin1_R = pin1[1,j,i] * utg_auxinRatioR
			updated_pin1_B = pin1[2,j,i] * utg_auxinRatioB
			updated_pin1_L = pin1[3,j,i] * utg_auxinRatioL

			total_pin1_updated = updated_pin1_T + updated_pin1_R + updated_pin1_B + updated_pin1_L

			normalization_ratio = total_pin1 / total_pin1_updated

			pin1[0,j,i] = updated_pin1_T * normalization_ratio
			pin1[1,j,i] = updated_pin1_R * normalization_ratio
			pin1[2,j,i] = updated_pin1_B * normalization_ratio
			pin1[3,j,i] = updated_pin1_L * normalization_ratio


	
	auxin[3,3] = auxin[3,3] + auxin_synthesis
	auxin[9,11] = auxin[9,11] + auxin_synthesis



