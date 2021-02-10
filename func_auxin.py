#!/usr/bin/env python

import numpy as np
import math

import params as pr
import inputs as ip


def auxin_homeostasis(iton):
	
	'''
	Here implement: basal synthesis and turnover, possible effect of CUC on YUC,
	local modifications like exogenous application, etc
	A' = h * ( Synth + custom_synth_degr + C*k(CY) - A*k(Cdecay) )

	Params:
	* iton: Iteration of the simulation. This is used for local auxin synth/degr
	
	'''

	# Simplify var names
	h = pr.euler_h
	k_auxin_synth = pr.k_auxin_synth
	k_auxin_degr = pr.k_auxin_degr
	th_cuc_yuc1 = pr.th_cuc_yuc1
	k_cuc_yuc1 = pr.k_cuc_yuc1
	k_cuc_yuc4 = pr.k_cuc_yuc4
	
	for y in range(ip.tissue_rows):
		for x in range(ip.tissue_columns):
			
			# Get auxin and cuc values for current cell
			auxin_cell = ip.auxin[y,x]
			cuc_cell = ip.cuc[y,x]

			# If current cell has local/custom auxin synth/degr...
			current_cell = (y,x)

			if current_cell in pr.auxin_custom_synth['cells'] and iton in pr.auxin_custom_synth['iterations']:
				local_synth = pr.auxin_custom_synth['value']
			else:
				local_synth = 0
			
			if current_cell in pr.auxin_custom_degr['cells'] and iton in pr.auxin_custom_synth['iterations']:
				local_degr = pr.auxin_custom_degr['value']
			else:
				local_degr = 0
			
			# Test: Effect of CUC on YUC1 follows a step function
			if cuc_cell >= th_cuc_yuc1:
				synth_yuc1 = cuc_cell * k_cuc_yuc1
			else:
				synth_yuc1 = 0
			
			# Calculate change in auxin concentration
			auxin_cell_updated = auxin_cell + h * ( \
				k_auxin_synth + \
				local_synth - \
				local_degr + \
				synth_yuc1 + \
				cuc_cell * k_cuc_yuc4 - \
				auxin_cell * k_auxin_degr \
			)
			
			if auxin_cell_updated < 0:
				auxin_cell_updated = 0
				
			ip.auxin[y,x] = auxin_cell_updated


def auxin_diffusion():
	
	#
	# [auxin](i,j) = [auxin](i,j) - out_diff + in_diff_T + in_diff_R + in_diff_B + in_diff_L
	#
	# Rate of diffusion from cell i to j: dD(i->j)/dt = h * ( [auxin(i)] * k )
	#
	# k = diffusion factor constant
	#

	#fluxes = np.zeros(shape=(8,tissue_rows,tissue_columns)) # Z: T_out, T_in, R_out, R_in...

	# Simplify var names
	h = pr.euler_h
	D = pr.k_auxin_diffusion
	auxin = ip.auxin
	fluxes = ip.auxin_fluxes_difusion
	tissue_rows = ip.tissue_rows
	tissue_columns = ip.tissue_columns

	for y in range(tissue_rows):
		for x in range(tissue_columns):
			
			# Calculate in an out fluxes for each cell face (molecules / face and cycle)
			# (calculating fluxes is necessary for WTF polarization)

			# Top face: out (fluxes[0,y,x]); in (fluxes[1,y,x])
			if y > 0:
				fluxes[0,y,x] = h * ( auxin[y,x] * D )
				fluxes[1,y,x]  = h * ( auxin[y-1,x] * D )
			else:
				fluxes[0,y,x], fluxes[1,y,x] = 0, 0

			# Right face: out (fluxes[2,y,x]); in (fluxes[3,y,x])
			if x < tissue_columns-1:
				fluxes[2,y,x] = h * ( auxin[y,x] * D )
				fluxes[3,y,x]  = h * ( auxin[y,x+1] * D )
			else:
				fluxes[2,y,x], fluxes[3,y,x] = 0, 0

			# Bottom face: out (fluxes[4,y,x]); in (fluxes[5,y,x])
			if y < tissue_rows-1:
				fluxes[4,y,x] = h * ( auxin[y,x] * D )
				fluxes[5,y,x]  = h * ( auxin[y+1,x] * D )
			else:
				fluxes[4,y,x], fluxes[5,y,x] = 0, 0

			# Left face: out (fluxes[6,y,x]); in (fluxes[7,y,x])
			if x > 0:
				fluxes[6,y,x] = h * ( auxin[y,x] * D )
				fluxes[7,y,x]  = h * ( auxin[y,x-1] * D )
			else:
				fluxes[6,y,x], fluxes[7,y,x] = 0, 0

			# Calculate vector components and angle of net flux to draw line

			# Calculate net fluxes (outbound, out - in) to draw vectors
			T_net_flux = fluxes[0,y,x] - fluxes[1,y,x]
			R_net_flux = fluxes[2,y,x] - fluxes[3,y,x]
			B_net_flux = fluxes[4,y,x] - fluxes[5,y,x]
			L_net_flux = fluxes[6,y,x] - fluxes[7,y,x]
			# Vector X component
			fluxes[8,y,x] = R_net_flux - L_net_flux
			# Vector Y component
			fluxes[9,y,x] = B_net_flux - T_net_flux

	# Update the auxin concentrations after calculating all the fluxes to avoid polarity effect of looping through numpy array
	# This could go inside auxin_homeostasis()
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			fluxes_out = fluxes[0,y,x] + fluxes[2,y,x] + fluxes[4,y,x] + fluxes[6,y,x]
			fluxes_in = fluxes[1,y,x] + fluxes[3,y,x] + fluxes[5,y,x] + fluxes[7,y,x]
			
			auxin[y,x] = auxin[y,x] - fluxes_out + fluxes_in


def pin_on_auxin(k_pin1_transp):

	#
	# PIN1-MEDIATED AUXIN TRANSPORT
	# Apply PIN1 transport to auxin concentration values
	#
	# PIN1_Tr = Nbr auxin molecules / ( PIN1 molecule * cycle )
	#
	# Transport rate = h * ( [auxin] * [PIN1] * k )
	#

	# Simplify var names
	h = pr.euler_h
	K = k_pin1_transp
	auxin = ip.auxin
	pin1 = ip.pin1
	#k_pin1_transp = pr.k_pin1_transp
	tissue_rows = ip.tissue_rows
	tissue_columns = ip.tissue_columns
	fluxes_pin1 = ip.auxin_fluxes_pin1
	
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			# Certain combinations of [PIN1] and K could result in that the sum of auxin molecules effluxed (through all the cell faces) exceeds the num of auxin molecules in the cell. To detect that, first calculate the efflux in the whole cell and check if it is higher than the num of auxin molecules available.
			total_pin1 = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
			transported_molecules_total = h * auxin[y,x] * total_pin1 * K 

			# Manually limit transport if it exceeds the amount of auxin molecules
			# Later on: calculate excess ratio and then use to reduce the value of the vector below
			if transported_molecules_total > auxin[y,x]:
				print('warning: transported molecules had to be manually adjusted in cell ' + str(y), str(x))

			# Calculate in an out fluxes for each cell face (molecules / face and cycle)
			# (calculating fluxes is necessary for WTF polarization)

			# Top face: out (fluxes_pin1[0,y,x]); in (fluxes_pin1[1,y,x])
			if y > 0:
				fluxes_pin1[0,y,x] = h * ( auxin[y,x] * pin1[0,y,x] * K )
				fluxes_pin1[1,y,x]  = h * ( auxin[y-1,x] * pin1[2,y-1,x] * K )
			else:
				fluxes_pin1[0,y,x], fluxes_pin1[1,y,x] = 0, 0
			
			# Right face: out (fluxes_pin1[2,y,x]); in (fluxes_pin1[3,y,x])
			if x < tissue_columns - 1:
				fluxes_pin1[2,y,x] = h * ( auxin[y,x] * pin1[1,y,x] * K )
				fluxes_pin1[3,y,x]  = h * ( auxin[y,x+1] * pin1[3,y,x+1] * K )
			else:
				fluxes_pin1[2,y,x], fluxes_pin1[3,y,x] = 0, 0
			
			# Bottom face: out (fluxes_pin1[4,y,x]); in (fluxes_pin1[5,y,x])
			if y < tissue_rows - 1:
				fluxes_pin1[4,y,x] = h * ( auxin[y,x] * pin1[2,y,x] * K )
				fluxes_pin1[5,y,x]  = h * ( auxin[y+1,x] * pin1[0,y+1,x] * K )
			else:
				fluxes_pin1[4,y,x], fluxes_pin1[5,y,x] = 0, 0
			
			# Left face: out (fluxes_pin1[6,y,x]); in (fluxes_pin1[7,y,x])
			if x > 0:
				fluxes_pin1[6,y,x] = h * ( auxin[y,x] * pin1[3,y,x] * K )
				fluxes_pin1[7,y,x]  = h * ( auxin[y,x-1] * pin1[1,y,x-1] * K )
			else:
				fluxes_pin1[6,y,x], fluxes_pin1[7,y,x] = 0, 0
			
			# Calculate vector components and angle of net flux to draw line

			# Calculate net fluxes (outbound, out - in) to draw vectors
			T_net_flux = fluxes_pin1[0,y,x] - fluxes_pin1[1,y,x]
			R_net_flux = fluxes_pin1[2,y,x] - fluxes_pin1[3,y,x]
			B_net_flux = fluxes_pin1[4,y,x] - fluxes_pin1[5,y,x]
			L_net_flux = fluxes_pin1[6,y,x] - fluxes_pin1[7,y,x]
			# Vector X component
			vector_x = fluxes_pin1[8,y,x] = R_net_flux - L_net_flux
			# Vector Y component
			vector_y = fluxes_pin1[9,y,x] = B_net_flux - T_net_flux
			# Calculate sine of angle described by vector
			vector_hyp = math.sqrt(vector_x**2 + vector_y**2)
			if vector_y != 0 and vector_hyp != 0:
				sine = vector_y / vector_hyp
				# Calculate degres of vector (convert sine to radians and then to degrees)
				fluxes_pin1[10,y,x] = math.degrees(math.asin(sine))
			else:
				fluxes_pin1[10,y,x] = 366.0

			#print(vector_y, vector_hyp, fluxes_pin1[10,y,x])

	# Update the auxin concentrations after calculating all the fluxes to avoid polarity effect of looping through numpy array
	# This could go inside auxin_homeostasis()
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			fluxes_pin1_out = fluxes_pin1[0,y,x] + fluxes_pin1[2,y,x] + fluxes_pin1[4,y,x] + fluxes_pin1[6,y,x]
			fluxes_pin1_in = fluxes_pin1[1,y,x] + fluxes_pin1[3,y,x] + fluxes_pin1[5,y,x] + fluxes_pin1[7,y,x]

			auxin[y,x] = auxin[y,x] - fluxes_pin1_out + fluxes_pin1_in


def pin_on_auxin_old(k_pin1_transp):

	#
	# PIN1-MEDIATED AUXIN TRANSPORT
	# Apply PIN1 transport to auxin concentration values
	#
	# PIN1_Tr = Nbr auxin molecules / ( PIN1 molecule * cycle )
	#
	# Transport rate = h * ( [auxin] * [PIN1] * k )
	#

	# Simplify var names
	h = pr.euler_h
	auxin = ip.auxin
	pin1 = ip.pin1
	#k_pin1_transp = pr.k_pin1_transp
	tissue_rows = ip.tissue_rows
	tissue_columns = ip.tissue_columns
	transpVectors = ip.auxin_fluxes_pin1
	
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			auxin_molecules = auxin[y,x]
			total_pin1 = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
			transported_molecules_total = auxin_molecules * total_pin1 * k_pin1_transp
			
			#excess_ratio = transported_molecules_total / 

			# Manually limit transport if it exceeds the amount of auxin molecules
			# Later on: calculate excess ratio and then use to reduce the value of the vector below
			if transported_molecules_total > auxin_molecules:
				transported_molecules_total = auxin_molecules
				print('warning: transported molecules had to be manually adjusted in cell ' + str(y), str(x))

			# Go through cell matrix and calculate all transport vectors (all are efflux)

			# To top (y,x -> y-1,x)
			if y > 0:
				transpVectors[0,y,x] = h * ( auxin_molecules * pin1[0,y,x] * k_pin1_transp )
			else:
				transpVectors[0,y,x] = 0

			# To right (y,x -> y,x+1)
			if x < tissue_columns - 1:
				transpVectors[1,y,x] = h * ( auxin_molecules * pin1[1,y,x] * k_pin1_transp )
			else:
				transpVectors[1,y,x] = 0

			# To bottom (y,x -> y+1,x)
			if y < tissue_rows - 1:
				transpVectors[2,y,x] = h * ( auxin_molecules * pin1[2,y,x] * k_pin1_transp )
			else:
				transpVectors[2,y,x] = 0

			# To left (y,x -> y,x-1)
			if x > 0:
				transpVectors[3,y,x] = h * ( auxin_molecules * pin1[3,y,x] * k_pin1_transp )
			else:
				transpVectors[3,y,x] = 0

	for y in range(tissue_rows):
		for x in range(tissue_columns):

			# From top (y,x <- y-1,x)
			if y > 0:
				transpFromTop = transpVectors[2,y-1,x]
			else:
				transpFromTop = 0

			# From right (y,x <- y,x+1)
			if x < tissue_columns - 1:
				transpFromRight = transpVectors[3,y,x+1]
			else:
				transpFromRight = 0 

			# From bottom (y,x <- y+1,x)
			if y < tissue_rows - 1:
				transpFromBottom = transpVectors[0,y+1,x]
			else:
				transpFromBottom = 0

			# From left (y,x <- y,x-1)
			if x > 0:
				transpFromLeft = transpVectors[1,y,x-1]
			else:
				transpFromLeft = 0

			# Update the auxin concentration in each cell
			total_efflux = transpVectors[0,y,x] + transpVectors[1,y,x] + transpVectors[2,y,x] + transpVectors[3,y,x]
			auxin[y,x] = auxin[y,x] - total_efflux + transpFromTop + transpFromRight + transpFromBottom + transpFromLeft


if __name__ == '__main__':
    pass




