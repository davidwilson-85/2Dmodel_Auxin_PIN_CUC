#!/usr/bin/env python

import numpy as np

def auxin_diffusion(auxin_diffusionFactor, gridShape, columns, rows, auxin):
	
	#
	# [auxin](i,j) = [auxin](i,j) - out_diff + in_diff_T + in_diff_R + in_diff_B + in_diff_L
	#
	# Rate of diffusion from cell i to j: dD(i->j)/dt = [auxin(i)] * k
	#
	# k = diffusion factor constant
	#

	# Create absolute efflux diffusion values for each cell
	diffusionVectors = np.zeros(gridShape, dtype=(float,1))	
	for y in range(rows):
		for x in range(columns):
			# Amount of auxin that is lost per cell face
			diffusionVectors[y,x] = auxin[y,x] * auxin_diffusionFactor
			
	#print diffusionVectors
	
	# Apply diffusion to auxin concentration values
	for y in range(rows):
		for x in range(columns):

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
			
			# Update the auxin concentration in each cell
			auxin[y,x] = auxin[y,x] - ( diffusionVectors[y,x] * nbr_cell_neighbours ) + diffusionFromLeft + diffusionFromRight + diffusionFromTop + diffusionFromBottom