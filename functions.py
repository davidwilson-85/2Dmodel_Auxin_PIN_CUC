#!/usr/bin/env python

import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os, shutil



def create_cell_plot(matrix_shape, auxin, auxin_range, lut_auxin, pin1, pin1_range, lut_pin1, iteration, diff_vectors):

	#
	# Create cell plot using PIL
	#
	
	im = Image.new('RGBA', size=(1100,1100))
	x_origin = 0
	y_origin = 0
	cellSide = 50
	y = y_origin

	for i in range(matrix_shape[0]):
		x = x_origin
		
		for j in range(matrix_shape[1]):
			
			# Draw cell outline and fill
			ImageDraw.Draw(im).polygon([(x,y),(x+cellSide,y),(x+cellSide,y+cellSide),(x,y+cellSide)], outline=(50,50,50,255), fill=index_to_rgb(lut_auxin, auxin[i,j], auxin_range))

			# Draw PIN1
			ImageDraw.Draw(im).line([(x+4,y+3),(x+cellSide-4,y+3)], fill=index_to_rgb(lut_pin1, pin1[0,i,j], pin1_range), width=3)
			ImageDraw.Draw(im).line([(x+cellSide-3,y+4),(x+cellSide-3,y+cellSide-4)], fill=index_to_rgb(lut_pin1, pin1[1,i,j], pin1_range), width=3)
			ImageDraw.Draw(im).line([(x+4,y+cellSide-3),(x+cellSide-4,y+cellSide-3)], fill=index_to_rgb(lut_pin1, pin1[2,i,j], pin1_range), width=3)
			ImageDraw.Draw(im).line([(x+3,y+4),(x+3,y+cellSide-4)], fill=index_to_rgb(lut_pin1, pin1[3,i,j], pin1_range), width=3)
			
			# Draw CUC		
			#ImageDraw.Draw(im).ellipse([(x+15,y+15),(x+cellSide-15,y+cellSide-15)], fill='green')

			# Draw auxin concentration
			ImageDraw.Draw(im).text((x+20,y+20), str(round(auxin[i,j],1)), fill=(255, 255, 0))

			# Draw auxin diffussion vector
			ImageDraw.Draw(im).line([(x+20,y+20),(x+20+5*diff_vectors[0,i,j],y+20+5*diff_vectors[1,i,j])], fill='white', width=2)
			
			x = x + cellSide
		
		y = y + cellSide

	x = x_origin
	y = y_origin

	#ImageDraw.Draw(im).line([(x+1,y+2),(x+cellSide-1,y+2)], fill='red', width=3)

	# Set destination folder
	dest_folder = 'images/test'

	# Cleanup destination folder
	if iteration == 0:
		shutil.rmtree(dest_folder) 
		os.mkdir(dest_folder)


	im.save(dest_folder + '/image' + str(iteration) +'.png')




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




def auxin_diffusion(auxin_diffusionFactor, gridShape, columns, rows, auxin, diff_vectors):
	
	#
	# [auxin](i,j) = [auxin](i,j) - out_diff + in_diff_T + in_diff_R + in_diff_B + in_diff_L
	#
	# Rate of diffusion from cell i to j: dD(i->j)/dt = [auxin(i)] * k
	#
	# k = diffusion factor constant
	#

	# Calculate absolute efflux diffusion values for each cell
	diffusionVectors = np.zeros(gridShape, dtype=(float,1))	
	for y in range(rows):
		for x in range(columns):
			# Amount of auxin that is lost per cell face
			diffusionVectors[y,x] = auxin[y,x] * auxin_diffusionFactor
	
	# Calculate absolute incoming auxin from each neighbour. Update the auxin concentration in each cell
	for y in range(rows):
		for x in range(columns):

			# Count how many neighbours the cell has, to calculate the amount of efflux diffusion
			nbr_cell_neighbours = 0
			
			# numpy arrays accept negative indexes, so testing if an index exists at the left and top of the grid, and at the right and bottom, must be done in a different way: 
			
			# From Left
			if y-1 >= 0:
				diffusionFromLeft = diffusionVectors[y-1,x]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromLeft = 0
			
			# From Right
			try: 
				diffusionFromRight = diffusionVectors[y+1,x]
				nbr_cell_neighbours	+=1
			except IndexError:
				diffusionFromRight = 0
			
			# From Top	
			if x-1 >= 0:
				diffusionFromTop = diffusionVectors[y,x-1]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromTop = 0
			
			# From Bottom
			try: 
				diffusionFromBottom = diffusionVectors[y,x+1]
				nbr_cell_neighbours	+=1
			except IndexError:
				diffusionFromBottom = 0

			# Update the auxin concentration in each cell
			auxin[y,x] = auxin[y,x] - ( diffusionVectors[y,x] * nbr_cell_neighbours ) + diffusionFromLeft + diffusionFromRight + diffusionFromTop + diffusionFromBottom

			# Calculate vectors
			dx = 0 + diffusionFromRight + 0 - diffusionFromLeft 
			dy = diffusionFromTop + 0 - diffusionFromBottom

			diff_vectors[0,y,x] = dx
			diff_vectors[1,y,x] = dy