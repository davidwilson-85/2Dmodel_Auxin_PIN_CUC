#!/usr/bin/env python

import tilings
import numpy as np
import matplotlib.pyplot as plt

#import PIL
#import glob

#tilings.draw_tiling(tilings.generate_hexagons,  filename='hexagons.png')
#data = np.full((5, 5), 0)
#data[2,2] = 100

data = np.loadtxt('2-D_template', delimiter=',', unpack=True)
matrix_shape = data.shape

diffusionFactor = 0.035    # Proportion of molecules that cross between two adjacent cells per cycle
synthesisRate = 0.0001     # Absolute amount of molecules synthesised per cycle
lossRate = 0.75

nbr_iterations = 500

# =======================================================

for j in range(nbr_iterations):

	# Plot data
	fig = plt.imshow(data, cmap="plasma", vmin=0, vmax=1)
	plt.savefig('graphs/image' + str(j) + '.png')
	plt.close()
	
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
			
			try: 
				diffusionLeft = diffusionVectors[i-1,j]
			except IndexError:
				diffusionLeft = 0
	
			try: 
				diffusionRight = diffusionVectors[i+1,j]
			except IndexError:
				diffusionRight = 0
				
			try: 
				diffusionTop = diffusionVectors[i,j-1]
			except IndexError:
				diffusionTop = 0
				
			try: 
				diffusionBottom = diffusionVectors[i,j+1]
			except IndexError:
				diffusionBottom = 0
	
			data[i,j] = data[i,j] - diffusionVectors[i,j] * 4 + diffusionLeft + diffusionRight + diffusionTop + diffusionBottom
			data[60,60] = data[60,60] + synthesisRate



quit()

fig = plt.imshow(data, cmap="plasma")
plt.savefig('image.png')
