#!/usr/bin/env python

import importlib, sys
import numpy as np
from PIL import Image
pr = importlib.import_module(sys.argv[1].split('.')[0], package=None)
#import params as pr

'''
# === LOAD TEMPLATE DATA

# Templates

template = '2D'

auxin_template = 'templates/' + template + '/template_auxin'
pin1_template = 'templates/' + template + '/template_pin1'
cuc_template = 'templates/' + template + '/template_cuc'
middle_domain_template = 'templates/' + template + '/template_middle_domain'
adab_domain_template = 'templates/' + template + '/template_adab_domain'

#auxin = np.loadtxt(auxin_template, delimiter=',', unpack=False)
auxin = pr.template_auxin
auxin = auxin * 10

auxin_matrix_shape = auxin.shape
tissue_rows, tissue_columns = auxin.shape[0], auxin.shape[1]

pin1 = np.loadtxt(pin1_template, delimiter=',', unpack=False).reshape((4,tissue_rows,tissue_columns)) # Format is [z,y,x]
pin1_matrix_shape = pin1.shape
cuc = np.loadtxt(cuc_template, delimiter=',', unpack=False)
middle_domain = np.loadtxt(middle_domain_template, delimiter=',', unpack=False)
adab_domain = np.loadtxt(adab_domain_template, delimiter=',', unpack=False)

'''

# Copy arrays to new variables to not update the params namespace (critical when running series of simulations)
auxin = pr.template_auxin * 10
auxin_matrix_shape = auxin.shape
tissue_rows, tissue_columns = auxin.shape[0], auxin.shape[1]
pin1 = pr.template_pin1 # Format is [z,y,x] = [cell face (T, R, B, L), cell row, cell column]
pin1_matrix_shape = pin1.shape
cuc = pr.template_cuc
middle_domain = pr.template_middle_domain
proximodistal_axis = pr.template_proximodistal_axis
#adab_domain = np.loadtxt(adab_domain_template, delimiter=',', unpack=False)

# Auxin fluxes: number of auxin molecules that cross between 2 cells in a given simultation step
# auxin_fluxes_diffusion: fluxes via diffusion
# auxin_fluxes_pin1: fluxes through PIN1 transporters
auxin_fluxes_diffusion = np.zeros(shape=(10,tissue_rows,tissue_columns)) # Z: T_out, T_in, R_out, R_in ... netXvector, netYvector 
auxin_fluxes_pin1 = np.zeros(shape=(12,tissue_rows,tissue_columns), dtype=(float)) # 3D array = (attribute, column, row)
# Fluxes history buffers changes in fluxes and gives stability to the model (see Bayer 2009)
auxin_fluxes_history = np.zeros(shape=(4,tissue_rows,tissue_columns))

#array_auxin_net_fluxes = np.zeros(shape=(2,tissue_rows,tissue_columns)) # where z[0] => dx and z[1] => dy

# Arrays to store values and track simulation
auxin_auxiliary = auxin * 0
auxin_allcells_historic = []
auxin_incr_allcells_historic = []

pin1_auxiliary = pin1 * 0
pin1_allcells_historic = []
pin1_incr_allcells_historic = []

cuc_auxiliary = cuc * 0
cuc_allcells_historic = []
cuc_incr_allcells_historic = []

# LUTs
lut_auxin = np.loadtxt('luts/lut_red.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_pin1 = np.loadtxt('luts/lut_green.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)
lut_cuc = np.loadtxt('luts/lut_magenta.csv', delimiter=',', unpack=True, dtype=('int'), skiprows=1)

# Images/art
im_arrow_25 = Image.open('art/arrow_white_17x17_25.png')
im_arrow_50 = Image.open('art/arrow_white_17x17_50.png')
im_arrow_75 = Image.open('art/arrow_white_17x17_75.png')
im_arrow_100 = Image.open('art/arrow_white_17x17_100.png')

'''
FOR REFERENCE:

matrix_shape[0] -> number of rows = y
matrix_shape[1] -> number of columns = x

auxin[6,11] = 5 refers to 7th row, 12th column
pin1[0,0,5] = 5 refers to top wall of cell in 1st row, 6th column

'''