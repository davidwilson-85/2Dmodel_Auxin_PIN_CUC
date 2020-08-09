#!/usr/bin/env python

import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os, shutil



def create_cell_plot(matrix_shape, auxin, auxin_range, lut_auxin, pin1, pin1_range, lut_pin1, iteration, diff_vectors, img_dest_folder):

	#
	# Create cell plot using PIL
	#
	
	im = Image.new('RGBA', size=(300,300))
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
			ImageDraw.Draw(im).line([(x+25,y+25),(x+25+100*diff_vectors[0,i,j],y+25-100*diff_vectors[1,i,j])], fill='white', width=2)
			#ImageDraw.Draw(im).line([(x+25,y+25),(x+25+100,y+25+100)], fill='white', width=2)
			
			x = x + cellSide
		
		y = y + cellSide

	x = x_origin
	y = y_origin

	#ImageDraw.Draw(im).line([(x+1,y+2),(x+cellSide-1,y+2)], fill='red', width=3)

	# Save image
	im.save(img_dest_folder + '/image' + str(iteration+1000) +'.png')




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




def auxin_diffusion(auxin_diffusionFactor, gridShape, columns, rows, auxin, diff_vectors, iteration):
	
	#
	# [auxin](i,j) = [auxin](i,j) - out_diff + in_diff_T + in_diff_R + in_diff_B + in_diff_L
	#
	# Rate of diffusion from cell i to j: dD(i->j)/dt = [auxin(i)] * k
	#
	# k = diffusion factor constant
	#


	auxin_fluxes = np.zeros(shape=(8,rows,columns)) # Z: T_out, T_in...

	
	for y in range(rows):
		for x in range(columns):
			
			# Top face
			if y > 0:
				T_out = auxin[x,y] * auxin_diffusionFactor
				T_in = auxin[x,y-1] * auxin_diffusionFactor
			else:
				T_out, T_in = 0, 0

			# Right face
			if x < columns-1:
				R_out = auxin[x,y] * auxin_diffusionFactor
				R_in = auxin[x+1,y] * auxin_diffusionFactor
			else:
				R_out, R_in = 0, 0

			# Bottom face
			if y < rows-1:
				B_out = auxin[x,y] * auxin_diffusionFactor
				B_in = auxin[x,y+1] * auxin_diffusionFactor
			else:
				B_out, B_in = 0, 0

			# Left face
			if x > 0:
				L_out = auxin[x,y] * auxin_diffusionFactor
				L_in = auxin[x-1,y] * auxin_diffusionFactor
			else:
				L_out, L_in = 0, 0

			auxin_fluxes[0,x,y] = T_out
			auxin_fluxes[1,x,y] = T_in
			auxin_fluxes[2,x,y] = R_out
			auxin_fluxes[3,x,y] = R_in
			auxin_fluxes[4,x,y] = B_out
			auxin_fluxes[5,x,y] = B_in
			auxin_fluxes[6,x,y] = L_out
			auxin_fluxes[7,x,y] = L_in


			# Calculate net fluxes (outbound) to draw vectors
			T_net_flux = T_out - T_in
			R_net_flux = R_out - R_in
			B_net_flux = B_out - B_in
			L_net_flux = L_out - L_in
			
			vector_x_component = L_net_flux + R_net_flux
			vector_y_component = T_net_flux + B_net_flux

			diff_vectors[0,y,x] = vector_x_component
			diff_vectors[1,y,x] = vector_y_component


	# Update the auxin concentrations after calculating all the fluxes to avoid polarity effect of looping through numpy array
	for y in range(rows):
		for x in range(columns):

			#auxin[x,y] = auxin[x,y] - (T_out + R_out + B_out + L_out) + (T_in + R_in + B_in + L_in)

			auxin[x,y] = auxin[x,y] - (auxin_fluxes[0,x,y] + auxin_fluxes[2,x,y] + auxin_fluxes[4,x,y] + auxin_fluxes[6,x,y]) + (auxin_fluxes[1,x,y] + auxin_fluxes[3,x,y] + auxin_fluxes[5,x,y] + auxin_fluxes[7,x,y])




			
	
	


def auxin_diffusion_old(auxin_diffusionFactor, gridShape, columns, rows, auxin, diff_vectors, iteration):
	
	#
	# [auxin](i,j) = [auxin](i,j) - out_diff + in_diff_T + in_diff_R + in_diff_B + in_diff_L
	#
	# Rate of diffusion from cell i to j: dD(i->j)/dt = [auxin(i)] * k
	#
	# k = diffusion factor constant
	#


	# Calculate absolute outgoing auxin to each neighbour (= per cell face)
	
	array_auxinLossPerCellFace = np.zeros(gridShape, dtype=(float,1))
	
	for y in range(rows):
		for x in range(columns):
			# Amount of auxin that is lost per cell face
			array_auxinLossPerCellFace[y,x] = auxin[y,x] * auxin_diffusionFactor
	
	# Calculate absolute incoming auxin from each neighbour (= per cell face). Then, update the auxin concentration in each cell
	
	for y in range(rows):
		print 'y', y, 'rows', rows
		for x in range(columns):
			
			# REWRITE THIS TO MAKE IT SIMPLER (BOUNDARY ALSO AFFECTS OUTGOING AUXIN)

			# Count how many neighbours the cell has, to calculate the amount of efflux diffusion
			nbr_cell_neighbours = 0
			
			# numpy arrays accept negative indexes, so testing if an index exists at the left and top of the grid, and at the right and bottom, must be done in a different way: 
			
			# From top	
			if y > 0:
				diffusionFromTop = array_auxinLossPerCellFace[y,x-1]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromTop = 0
			
			# From right
			if x < columns - 1:
				diffusionFromRight = array_auxinLossPerCellFace[y+1,x]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromRight = 0
			
			# From bottom
			if y < rows - 1:
				diffusionFromBottom = array_auxinLossPerCellFace[y,x+1]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromBottom = 0
			
			# From left
			if x > 0:
				diffusionFromLeft = array_auxinLossPerCellFace[y-1,x]
				nbr_cell_neighbours	+=1
			else:
				diffusionFromLeft = 0

			# Update the auxin concentration in each cell
			auxin[y,x] = auxin[y,x] - ( array_auxinLossPerCellFace[y,x] * nbr_cell_neighbours ) + diffusionFromLeft + diffusionFromRight + diffusionFromTop + diffusionFromBottom

			
			# Calculate direction vectors
			
			# BUG HERE: I NEED TO CHECK IF THERE IS NEIGHBOR BEFORE CALCULATING AUXIN LOSS.......
			
			# Net flux top face
			netFlux_T_face = array_auxinLossPerCellFace[y,x] - diffusionFromTop
			
			# Net flux right face
			netFlux_R_face = array_auxinLossPerCellFace[y,x] - diffusionFromRight
			
			# Net flux bottom face
			netFlux_B_face = - array_auxinLossPerCellFace[y,x] + diffusionFromBottom
			
			# Net flux left face
			netFlux_L_face = - array_auxinLossPerCellFace[y,x] + diffusionFromLeft
			
			x_component = netFlux_L_face + netFlux_R_face
			y_component = netFlux_T_face + netFlux_B_face

			diff_vectors[0,y,x] = x_component
			diff_vectors[1,y,x] = y_component
			
			if y == 1 and x == 0 and iteration == 0:
				print 'diffusionFromTop', diffusionFromTop
				print netFlux_T_face, netFlux_R_face, netFlux_B_face, netFlux_L_face			
