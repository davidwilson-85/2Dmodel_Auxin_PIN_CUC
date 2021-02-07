#!/usr/bin/env python

import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os, shutil
import cv2

import inputs as ip


def create_cell_plot(matrix_shape, auxin, auxin_range, lut_auxin, pin1, pin1_range, lut_pin1, cuc, cuc_range, lut_cuc, iteration, array_afd, array_afp, img_dest_folder):

	#
	# Create cell plot using PIL
	#
	
	draw_auxin = True
	draw_pin = True
	draw_cuc = False
	draw_values_text = True
	draw_vectors_diff = True
	draw_vectors_pin1 = True

	# Vector magnification factor (only changes visualization)
	vector_mag = 50

	#img = Image.new('RGB', (100, 100))
	#draw = ImageDraw.Draw(img, 'RGBA')
	#drw.polygon([(50, 0), (100, 100), (0, 100)], (255, 0, 0, 125))

	x_origin = 0
	y_origin = 0
	cellSide = 50

	# Calculate required image dimensions
	height = x_origin + ip.tissue_rows * cellSide
	width = y_origin + ip.tissue_columns * cellSide
	
	im = Image.new('RGB', size=(width,height))
	draw = ImageDraw.Draw(im, 'RGBA')
	
	y = y_origin

	for i in range(ip.tissue_rows):
		
		x = x_origin
		for j in range(ip.tissue_columns):
			
			# Draw cell outline and auxin fill
			if draw_auxin == True:
				draw.polygon([(x,y),(x+cellSide,y),(x+cellSide,y+cellSide),(x,y+cellSide)], outline=(50,50,50,255), fill=index_to_rgb(lut_auxin, auxin[i,j], auxin_range))

			# Draw PIN1
			if draw_pin == True:
				draw.line([(x+4,y+3),(x+cellSide-4,y+3)], fill=index_to_rgb(lut_pin1, pin1[0,i,j], pin1_range), width=3)
				draw.line([(x+cellSide-3,y+4),(x+cellSide-3,y+cellSide-4)], fill=index_to_rgb(lut_pin1, pin1[1,i,j], pin1_range), width=3)
				draw.line([(x+4,y+cellSide-3),(x+cellSide-4,y+cellSide-3)], fill=index_to_rgb(lut_pin1, pin1[2,i,j], pin1_range), width=3)
				draw.line([(x+3,y+4),(x+3,y+cellSide-4)], fill=index_to_rgb(lut_pin1, pin1[3,i,j], pin1_range), width=3)
			
			if draw_cuc == True:
				# Draw CUC
				# outline=None or outline with same value as fill does not work. To solve this I use the same color as the auxin citoplasm background.	
				draw.ellipse([(x+15,y+15),(x+cellSide-15,y+cellSide-15)], outline=index_to_rgb(lut_auxin, auxin[i,j], auxin_range), fill=index_to_rgba(lut_cuc, cuc[i,j], cuc_range))
			
			if draw_values_text == True:
				# Write auxin concentration (magenta)
				ImageDraw.Draw(im).text((x+10,y+5), str(round(auxin[i,j],1)), fill=(255, 0, 255))
				# Write PIN1 total concentration (yellow)
				#ImageDraw.Draw(im).text((x+10,y+20), str(round( (pin1[0,i,j]+pin1[1,i,j]+pin1[2,i,j]+pin1[3,i,j]) ,1)), fill=(255, 255, 0))
				# Write CUC concentration (white)
				#ImageDraw.Draw(im).text((x+10,y+35), str(round(cuc[i,j],1)), fill=(255, 255, 255))
			
			x = x + cellSide
		
		y = y + cellSide

	x = x_origin
	y = y_origin

	if draw_vectors_diff == True:
		for i in range(matrix_shape[0]):
			
			x = x_origin
			for j in range(matrix_shape[1]):
	
				# Draw auxin diffussion vector
				draw.line([(x+25,y+25),(x+25+vector_mag*array_afd[8,i,j],y+25+vector_mag*array_afd[9,i,j])], fill='white', width=2)
				
				x = x + cellSide
			
			y = y + cellSide
	
		x = x_origin
		y = y_origin
	
	if draw_vectors_pin1 == True:
		for i in range(matrix_shape[0]):
			
			x = x_origin
			for j in range(matrix_shape[1]):
	
				# Draw PIN1 transport vector
				draw.line([(x+25,y+25),(x+25+vector_mag*array_afp[8,i,j],y+25+vector_mag*array_afp[9,i,j])], fill='yellow', width=3)
				
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


# Create video from images
def create_video():

	'''
	Check also:
	https://stackoverflow.com/questions/44947505/how-to-make-a-movie-out-of-images-in-python
	http://tsaith.github.io/combine-images-into-a-video-with-python-3-and-opencv-3.html
	https://theailearner.com/2018/10/15/creating-video-from-images-using-opencv-python/

	To call this function:
	$ python3 func_graph.py
	'''

	image_folder = 'images/test'
	video_name = 'videos/video7.avi'

	images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
	images.sort()

	frame = cv2.imread(os.path.join(image_folder, images[0]))
	height, width, layers = frame.shape
	print(height, width)

	video = cv2.VideoWriter(video_name, 0, 10, (width,height))

	for image in images:
		video.write(cv2.imread(os.path.join(image_folder, image)))

	cv2.destroyAllWindows()
	video.release()

if __name__ == '__main__':
	create_video()


