#!/usr/bin/env python

#import tilings
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

#import PIL
#import glob

#tilings.draw_tiling(tilings.generate_hexagons,  filename='hexagons.png')
#data = np.full((5, 5), 0)
#data[2,2] = 100


# === SET PARAMETERS

template = 'templates/2Dcell_template'
luts_to_load = ('luts/lut_gem.csv','luts/lut_fire.csv')

# Parameters for auxin
auxin_range = (0, 1)
diffusionFactor = 0.05		# Proportion of molecules that cross between two adjacent cells per cycle
synthesis = 0.0005	     	# Absolute amount of molecules synthesised per cycle
lossRate = 0.75
nbr_iterations = 20


# === LOAD DATA

data = np.loadtxt(template, delimiter=',', unpack=True)
matrix_shape = data.shape

lut1 = np.loadtxt(luts_to_load[1], delimiter=',', unpack=True, dtype=('int'))


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
	
	im = Image.new('RGB', size=(501,501))
	x_origin = 0
	y_origin = 0
	cellSide = 50
	y = y_origin

	for i in range(matrix_shape[0]):
		x = x_origin
		
		for j in range(matrix_shape[1]):
			
			# Draw cell outline
			ImageDraw.Draw(im).polygon([(x,y),(x+cellSide,y),(x+cellSide,y+cellSide),(x,y+cellSide)], outline='white', fill=index_to_rgb(lut1, data[i,j], auxin_range))
			# Draw PIN1
			#ImageDraw.Draw(im).line([(x+1,y+2),(x+cellSide-1,y+2)], fill='red', width=3)
			# Draw CUC		
			#ImageDraw.Draw(im).ellipse([(x+15,y+15),(x+cellSide-15,y+cellSide-15)], fill='green')
			
			x = x + cellSide
		
		y = y + cellSide

	x = x_origin
	y = y_origin

	#ImageDraw.Draw(im).line([(x+1,y+2),(x+cellSide-1,y+2)], fill='red', width=3)

	im.save('images/test2/image' + str(iteration) +'.png')



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
			#print data[i,j]
			diffusionVectors[i,j] = data[i,j] * diffusionFactor
			
	#print diffusionVectors
	
	
	# Calculate updated concentration values
	for i in range(matrix_shape[0]):
		for j in range(matrix_shape[1]):

			nbr_cell_neighbours = 0
			
			try: 
				diffusionFromLeft = diffusionVectors[i-1,j]
				nbr_cell_neighbours	
			except IndexError:
				diffusionLeft = 0
	
			try: 
				diffusionFromRight = diffusionVectors[i+1,j]
			except IndexError:
				diffusionRight = 0
				
			try: 
				diffusionFromTop = diffusionVectors[i,j-1]
			except IndexError:
				diffusionTop = 0
				
			try: 
				diffusionFromBottom = diffusionVectors[i,j+1]
			except IndexError:
				diffusionBottom = 0
	
			data[i,j] = data[i,j] - diffusionVectors[i,j] * 4 + diffusionFromLeft + diffusionFromRight + diffusionFromTop + diffusionFromBottom
			data[4,4] = data[4,4] + synthesis



