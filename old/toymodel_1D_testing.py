#!/usr/bin/env python

import tilings
import numpy as np
import matplotlib.pyplot as plt

#import PIL
#import glob

#tilings.draw_tiling(tilings.generate_hexagons,  filename='hexagons.png')
#data = np.full((5, 5), 0)
#data[2,2] = 100

data = np.loadtxt('1-D_template', delimiter=',', unpack=True)
nbr_cells = len(data)

diffusionFactor = 0.025

synthesisRate = 0.05

lossRate = 0.75

nbr_iterations = 1000

# =======================================================

for j in range(nbr_iterations):

	# Plot data
	fig = plt.figure()
	ax = fig.add_axes([0,0,1,1])
	x = range(nbr_cells)
	ax.set_ylim((0,5))
	ax.bar(x, data)
	plt.savefig('graphs/image' + str(j) + '.png')
	plt.close()
	
	# Create absolute diffusion efflux values for each cell
	diffusionVectors = np.zeros(nbr_cells, dtype=(float,2))
	#vectors = np.zeros((5,5), dtype=(float,4))
	for i in range(nbr_cells):
		diffusionVectors[i][0] = data[i] * diffusionFactor
		diffusionVectors[i][1] = data[i] * diffusionFactor

	#print data
	#print vectors
	
	# Calculate updated concentration values
	for i in range(nbr_cells):
	
		try: 
			diffusionLeft = diffusionVectors[i-1][1]
		except IndexError:
			diffusionLeft= 0
	
		try: 
			diffusionRight = diffusionVectors[i+1][0]
		except IndexError:
			diffusionRight = 0
	
		data[i] = data[i] - diffusionVectors[i][0] - diffusionVectors[i][1] + diffusionLeft + diffusionRight
		
	data[0] = data[0] + synthesisRate
	data[15] = data[15] - data[15] * lossRate




quit()

fig = plt.imshow(data, cmap="plasma")
plt.savefig('image.png')
