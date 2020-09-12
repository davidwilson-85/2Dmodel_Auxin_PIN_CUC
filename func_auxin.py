#!/usr/bin/env python

import numpy as np



def auxin_homeostasis(auxin, cuc, k_auxin_synth, k_cuc_yuc, k_auxin_decay):
	
	#
	# Here implement: basal synthesis and turnover, possible effect of CUC on YUC,
	# local modifications like exogenous application, etc
	#
	# A' = Synth + C*k(CA) - A*k(Cdecay)
	#
	
	for y in range(auxin.shape[0]):
		for x in range(auxin.shape[1]):
			
			auxin_cell = auxin[y,x]
			cuc_cell = cuc[y,x]
			
			auxin_cell_updated = auxin_cell + k_auxin_synth + cuc_cell * k_cuc_yuc - auxin_cell * k_auxin_decay
			
			if auxin_cell_updated < 0:
				auxin_cell_updated = 0
				
			auxin[y,x] = auxin_cell_updated



def auxin_diffusion(k_auxin_diffusion, gridShape, tissue_columns, tissue_rows, auxin, array_af, iteration):
	
	#
	# [auxin](i,j) = [auxin](i,j) - out_diff + in_diff_T + in_diff_R + in_diff_B + in_diff_L
	#
	# Rate of diffusion from cell i to j: dD(i->j)/dt = [auxin(i)] * k
	#
	# k = diffusion factor constant
	#

	#array_af = np.zeros(shape=(8,tissue_rows,tissue_columns)) # Z: T_out, T_in, R_out, R_in...

	for y in range(tissue_rows):
		for x in range(tissue_columns):
			
			# Top face
			if y > 0:
				T_out = auxin[y,x] * k_auxin_diffusion
				T_in = auxin[y-1,x] * k_auxin_diffusion
			else:
				T_out, T_in = 0, 0

			# Right face
			if x < tissue_columns-1:
				R_out = auxin[y,x] * k_auxin_diffusion
				R_in = auxin[y,x+1] * k_auxin_diffusion
			else:
				R_out, R_in = 0, 0

			# Bottom face
			if y < tissue_rows-1:
				B_out = auxin[y,x] * k_auxin_diffusion
				B_in = auxin[y+1,x] * k_auxin_diffusion
			else:
				B_out, B_in = 0, 0

			# Left face
			if x > 0:
				L_out = auxin[y,x] * k_auxin_diffusion
				L_in = auxin[y,x-1] * k_auxin_diffusion
			else:
				L_out, L_in = 0, 0

			array_af[0,y,x] = T_out
			array_af[1,y,x] = T_in
			array_af[2,y,x] = R_out
			array_af[3,y,x] = R_in
			array_af[4,y,x] = B_out
			array_af[5,y,x] = B_in
			array_af[6,y,x] = L_out
			array_af[7,y,x] = L_in


			# Calculate net fluxes (outbound) to draw vectors
			T_net_flux = T_out - T_in
			R_net_flux = R_out - R_in
			B_net_flux = B_out - B_in
			L_net_flux = L_out - L_in

			# Change sign of T and L net fluxes
			T_net_flux = - T_net_flux
			L_net_flux = - L_net_flux
			
			vector_x_component = L_net_flux + R_net_flux
			vector_y_component = T_net_flux + B_net_flux

			array_af[8,y,x] = vector_x_component
			array_af[9,y,x] = vector_y_component

	#print array_af[0:8,:,:]


	# Update the auxin concentrations after calculating all the fluxes to avoid polarity effect of looping through numpy array
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			#auxin[x,y] = auxin[x,y] - (T_out + R_out + B_out + L_out) + (T_in + R_in + B_in + L_in)

			auxin[y,x] = auxin[y,x] - (array_af[0,y,x] + array_af[2,y,x] + array_af[4,y,x] + array_af[6,y,x]) + (array_af[1,y,x] + array_af[3,y,x] + array_af[5,y,x] + array_af[7,y,x])





def pin_on_auxin(auxin, pin1, k_pin1_transp, tissue_rows, tissue_columns, pin1_matrix_shape):

	#
	# PIN1-MEDIATED AUXIN TRANSPORT
	# Apply PIN1 transport to auxin concentration values
	#
	# PIN1_Tr = Nbr auxin molecules / ( PIN1 molecule * cycle )
	#
	# Transport rate = [auxin] * [PIN1] * k
	#

	# Create absolute efflux transport values for each cell
	transpVectors = np.zeros(pin1_matrix_shape, dtype=(float,1)) # 3D array = (cell_face, column, row)
	
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
				print 'warning: transported molecules had to be manually adjusted in cell ' + str(y), str(x)

			# To top (y,x -> y-1,x)
			if y > 0:
				transpVectors[0,y,x] = auxin_molecules * pin1[0,y,x] * k_pin1_transp
			else:
				transpVectors[0,y,x] = 0

			# To right (y,x -> y,x+1)
			if x < tissue_columns - 1:
				transpVectors[1,y,x] = auxin_molecules * pin1[1,y,x] * k_pin1_transp
			else:
				transpVectors[1,y,x] = 0

			# To bottom (y,x -> y+1,x)
			if y < tissue_rows - 1:
				transpVectors[2,y,x] = auxin_molecules * pin1[2,y,x] * k_pin1_transp
			else:
				transpVectors[2,y,x] = 0

			# To left (y,x -> y,x-1)
			if x > 0:
				transpVectors[3,y,x] = auxin_molecules * pin1[3,y,x] * k_pin1_transp
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




