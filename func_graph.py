#!/usr/bin/env python

import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os, shutil



def create_cell_plot(matrix_shape, auxin, auxin_range, lut_auxin, pin1, pin1_range, lut_pin1, cuc, cuc_range, lut_cuc, iteration, array_af, img_dest_folder):

	#
	# Create cell plot using PIL
	#

	# Vector magnification factor (only changes visualization)
	vector_mag = 10

	#img = Image.new('RGB', (100, 100))
	#draw = ImageDraw.Draw(img, 'RGBA')
	#drw.polygon([(50, 0), (100, 100), (0, 100)], (255, 0, 0, 125))
	
	im = Image.new('RGB', size=(400,700))
	draw = ImageDraw.Draw(im, 'RGBA')
	x_origin = 0
	y_origin = 0
	cellSide = 50
	y = y_origin

	for i in range(matrix_shape[0]):
		
		x = x_origin
		for j in range(matrix_shape[1]):
			
			# Draw cell outline and auxin fill
			draw.polygon([(x,y),(x+cellSide,y),(x+cellSide,y+cellSide),(x,y+cellSide)], outline=(50,50,50,255), fill=index_to_rgb(lut_auxin, auxin[i,j], auxin_range))

			# Draw PIN1
			draw.line([(x+4,y+3),(x+cellSide-4,y+3)], fill=index_to_rgb(lut_pin1, pin1[0,i,j], pin1_range), width=3)
			draw.line([(x+cellSide-3,y+4),(x+cellSide-3,y+cellSide-4)], fill=index_to_rgb(lut_pin1, pin1[1,i,j], pin1_range), width=3)
			draw.line([(x+4,y+cellSide-3),(x+cellSide-4,y+cellSide-3)], fill=index_to_rgb(lut_pin1, pin1[2,i,j], pin1_range), width=3)
			draw.line([(x+3,y+4),(x+3,y+cellSide-4)], fill=index_to_rgb(lut_pin1, pin1[3,i,j], pin1_range), width=3)

			# Write auxin concentration
			#ImageDraw.Draw(im).text((x+20,y+20), str(round(auxin[i,j],1)), fill=(255, 255, 0))

			# Draw CUC
			# outline=None or outline with same value as fill does not work. To solve this I use the same color as the auxin citoplasm background.	
			draw.ellipse([(x+15,y+15),(x+cellSide-15,y+cellSide-15)], outline=index_to_rgb(lut_auxin, auxin[i,j], auxin_range), fill=index_to_rgba(lut_cuc, cuc[i,j], cuc_range))

			#fill=(0, 0, 255, 25)
			
			x = x + cellSide
		
		y = y + cellSide

	x = x_origin
	y = y_origin


	for i in range(matrix_shape[0]):
		
		x = x_origin
		for j in range(matrix_shape[1]):

			# Draw auxin diffussion vector
			draw.line([(x+25,y+25),(x+25+vector_mag*array_af[8,i,j],y+25+vector_mag*array_af[9,i,j])], fill='white', width=2)
			
			x = x + cellSide
		
		y = y + cellSide

	x = x_origin
	y = y_origin
	
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



# Maps an integer representing the amount of a magnitude (e.g. [auxin]) and translates it to the corresponding RGBA triplet of the selected LUT
def index_to_rgba(lut, level, range):

	# Clip values out of range
	if level < range[0]: level = range[0]
	if level > range[1]: level = range[1]

	# Rescale range to 0-255 (this is typical lut range)
	rescaled_level = int( ( level / range[1] ) * 255 )
	
	# Create RGBA tuple
	rgba = (lut[1,rescaled_level],lut[2,rescaled_level],lut[3,rescaled_level],lut[0,rescaled_level])
	#print rgba
	return rgba



# Create a heatmap using matplotlib's imshow()
def create_heatmap(data):
	
	fig = plt.imshow(data, cmap="plasma", vmin=0, vmax=1, interpolation='none')
	plt.savefig('images/test/image' + str(iteration) + '.png')
	plt.close()

