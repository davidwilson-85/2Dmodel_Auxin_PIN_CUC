#!/usr/bin/env python

import random, importlib, sys, math
import numpy as np
from scipy.integrate import odeint

#import params as pr
pr = importlib.import_module(sys.argv[1].split('.')[0], package=None)
import inputs as ip


def auxin_custom_manipulation(iteration, sim_time):
	
	'''
	This function implements some changes in auxin:
	- Local exogenous application, etc
	- Noise: random variation in concentration

	Inputs and outputs:
	- Function does not have parameter inputs. It reads pr.template_auxin
	- Function does not return any objects. It writes to pr.template_auxin
	
	Params:
	- sim_time: Iteration of the simulation. This is used for custom (local) auxin synth/degr
	
	'''

	# Define simpler aliases
	h = pr.euler_h
		
	for y in range(ip.tissue_rows):
		for x in range(ip.tissue_columns):
			
			# Simplify var names
			auxin_cell = ip.auxin[y,x]

			current_cell = (y,x)
			custom_synth, custom_degr, noise = 0, 0, 0

			# Perfect sources and sinks (overwrites noise and synthesis/degradation)
			if pr.auxin_perfect_sources['active'] == True:
				if current_cell in pr.auxin_perfect_sources['cells'] \
				and sim_time >= pr.auxin_perfect_sources['time_interval'][0] \
				and sim_time < pr.auxin_perfect_sources['time_interval'][1]:
					auxin_cell = pr.auxin_perfect_sources['value']
			if pr.auxin_perfect_sinks['active'] == True:
				if current_cell in pr.auxin_perfect_sinks['cells'] \
				and sim_time >= pr.auxin_perfect_sinks['time_interval'][0] \
				and sim_time < pr.auxin_perfect_sinks['time_interval'][1]:
					auxin_cell = pr.auxin_perfect_sinks['value']
			
			# Local/custom auxin synth/degr
			if pr.auxin_custom_synth['value'] > 0:
				if current_cell in pr.auxin_custom_synth['cells'] \
				and sim_time >= pr.auxin_custom_synth['time_interval'][0] \
				and sim_time < pr.auxin_custom_synth['time_interval'][1]:
					custom_synth = pr.auxin_custom_synth['value']
			if pr.auxin_custom_degr['value'] > 0:
				if current_cell in pr.auxin_custom_degr['cells'] \
				and sim_time >= pr.auxin_custom_degr['time_interval'][0] \
				and sim_time < pr.auxin_custom_degr['time_interval'][1]:
					custom_degr = pr.auxin_custom_degr['value'] * auxin_cell
			
			# Noise
			if pr.auxin_noise['limit'] > 0 \
			and iteration >= pr.auxin_noise['iteration_interval'][0] \
			and iteration < pr.auxin_noise['iteration_interval'][1]:
				# Convert % to absolute limits for current cell, then take random number within limits 
				noise_lim_cell = auxin_cell * ( pr.auxin_noise['limit'] / 100 )
				noise = random.uniform(-noise_lim_cell, noise_lim_cell)
				# This is old system using predefined absolute limits
				#noise = random.uniform(-pr.auxin_noise['limit'], pr.auxin_noise['limit'])
				#print(noise)

			# Calculate change in auxin concentration and correct negative values if necessary
			auxin_cell_updated = auxin_cell + h * ( custom_synth - custom_degr ) + noise
			#if math.isnan(auxin_cell_updated) == True: print(auxin_cell_updated)
			if auxin_cell_updated < 0:
				auxin_cell_updated = float(1E-6)
			
			# Update auxin array
			ip.auxin[y,x] = auxin_cell_updated


def model_auxin_movement(init_values, t):

	'''
	This model is solved by ODEint.
	It only defines changes that can be described as differential equations and that occur within each individual cell.
	Changes not included: those that involve movement of auxin between cells, PIN1 polarization 
	'''

	A, y, x = init_values

	At = ip.auxin[ y-1 , x   ]
	Ar = ip.auxin[ y   , x+1 ]
	Ab = ip.auxin[ y+1 , x   ]
	Al = ip.auxin[ y   , x-1 ]

	#E = pr.k_pin1_effi_basal + pr.k_pin1_effi_cuc * ip.cuc[y,x]
	#Et = pr.k_pin1_effi_basal + pr.k_pin1_effi_cuc * ip.cuc[ y-1 , x   ]
	#Er = pr.k_pin1_effi_basal + pr.k_pin1_effi_cuc * ip.cuc[ y   , x+1 ]
	#Eb = pr.k_pin1_effi_basal + pr.k_pin1_effi_cuc * ip.cuc[ y+1 , x   ]
	#El = pr.k_pin1_effi_basal + pr.k_pin1_effi_cuc * ip.cuc[ y   , x-1 ]

	#Pt = ip.pin1[0, y, x]
	#Pr = ip.pin1[1, y, x]
	#Pb = ip.pin1[2, y, x]
	#Pl = ip.pin1[3, y, x]
	
	#Ptb = ip.pin1[ y-1 , x   ]
	#Prl = ip.pin1[ y   , x+1 ]
	#Pbt = ip.pin1[ y+1 , x   ]
	#Plr = ip.pin1[ y   , x-1 ]
	
	D = pr.k_auxin_diffusion

	#PIN1_net = Ptb*Et*At + Prl*Er*Ar + Pbt*Eb*Ab + Plr*El*Al - (Pt+Pr+Pb+Pl)*E*A
	#diff_net = (At+Ar+Ab+Al)*D - A*D
	#dA_dt = PIN1_net + diff_net
	
	dA_dt = (At+Ar+Ab+Al)*D - A*D
	y_dummy = 0
	x_dummy = 0

	return [dA_dt, y_dummy, x_dummy]


def solve_am_model():
		
	"""
	Function to calculate auxin movement (fluxes) via diffusion and active transport. Active transport is for now only PIN1 export.
	"""
	
	# Solve in a cell-by-cell basis
	for y in range(ip.tissue_rows):
		for x in range(ip.tissue_columns):
			
			# Gather initial values for ODEint
			model_init_values = [
				ip.auxin[y,x],
				y,
				x
			]

			# Solve
			cell_solution = odeint(model_auxin_movement, model_init_values, np.linspace(0, pr.euler_h, 2))
			
			# Update current cell in data arrays with solution output
			ip.auxin[y,x] = cell_solution[-1,0]
			

def auxin_diffusion():
	
	#
	# [auxin](i,j) = [auxin](i,j) - out_diff + in_diff_T + in_diff_R + in_diff_B + in_diff_L
	#
	# Rate of diffusion from cell i to j: dD(i->j)/dt = h * ( [auxin(i)] * k )
	#
	# k = diffusion constant
	#

	#fluxes = np.zeros(shape=(8,tissue_rows,tissue_columns)) # Z: T_out, T_in, R_out, R_in...

	# Simplify var names
	h = pr.euler_h
	D = pr.k_auxin_diffusion
	auxin = ip.auxin
	fluxes = ip.auxin_fluxes_diffusion
	tissue_rows = ip.tissue_rows
	tissue_columns = ip.tissue_columns

	for y in range(tissue_rows):
		for x in range(tissue_columns):
			
			# Calculate in an out fluxes for each cell face (molecules / face and cycle)
			# (calculating fluxes is necessary for WTF polarization)

			# Test
			#D = D * math.fabs(9 - (ip.middle_domain[x]) / 10)

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

			# Calculate vector components to draw line and to orient image with arrow direction

			# Calculate net fluxes (outbound, out - in) to draw vectors
			T_net_flux = fluxes[0,y,x] - fluxes[1,y,x]
			R_net_flux = fluxes[2,y,x] - fluxes[3,y,x]
			B_net_flux = fluxes[4,y,x] - fluxes[5,y,x]
			L_net_flux = fluxes[6,y,x] - fluxes[7,y,x]
			# Vector X component
			fluxes[8,y,x] = R_net_flux - L_net_flux
			# Vector Y component
			fluxes[9,y,x] = B_net_flux - T_net_flux

			# Angle of net flux is calculated in functions func_graph.create_cell_plot() and func_graph.vector_to_degrees()

	# Update the auxin concentrations after calculating all the fluxes to avoid polarity effect of looping through numpy array
	# This could go inside auxin_homeostasis()
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			fluxes_out = fluxes[0,y,x] + fluxes[2,y,x] + fluxes[4,y,x] + fluxes[6,y,x]
			fluxes_in = fluxes[1,y,x] + fluxes[3,y,x] + fluxes[5,y,x] + fluxes[7,y,x]
			
			auxin[y,x] = auxin[y,x] - fluxes_out + fluxes_in


def pin_on_auxin_ALTERNATIVE():

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
	cuc = ip.cuc
	Kp = pr.k_pin1_effi_basal
	Kp_p = pr.k_pin1_effi_pho # Efficiency of phosphorylated PIN1
	Kcp = pr.k_pin1_effi_cuc
	tissue_rows = ip.tissue_rows
	tissue_columns = ip.tissue_columns
	fluxes_pin1 = ip.auxin_fluxes_pin1

	# Calculate efficiencies of each cell:
	# Create array with same shape as tissue to store the PIN1 efficiency value of each cell
	# If model is UTGeff: Calculate proportion of phosphorylated PIN1 (high efficiency) in range 0-1, for cell y,x and the neighbours. Then calculate the combined PIN1 efficiency of cell y,x and the neighbours.
	
	effs = np.zeros(shape=(tissue_rows, tissue_columns))
	
	if pr.pin1_polarity == 'dual':
		effs.fill(Kp)
	if pr.pin1_polarity == 'utg_smith2006':
		c = 0.5 # Factor to use in function to calculate CUC-dependent proportion of phosphorylated PIN1
		for y in range(tissue_rows):
			for x in range(tissue_columns):
				pin1_p = 10*cuc[y,x] / (10*cuc[y,x] + c)
				effs[y,x] = Kp * (1-pin1_p) + Kp_p * pin1_p
				print(cuc[y,x], effs[y,x])
	
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			# Assign PIN1 efficiency of cell y,x and its neighbours to explicit var names
			K_out = effs[y,x]
			if y > 0:
				K_fromT = effs[y-1,x]
			if x < tissue_columns - 1:
				K_fromR = effs[y,x+1]
			if y < tissue_rows - 1:
				K_fromB = effs[y+1,x]
			if x > 0:
				K_fromL = effs[y,x-1]
		
			# Certain combinations of [PIN1] and K could result in that the sum of auxin molecules effluxed (through all the cell faces) exceeds the num of auxin molecules in the cell. To detect that, first calculate the efflux in the whole cell and check if it is higher than the num of auxin molecules available.
			total_pin1 = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
			transported_molecules_total = h * auxin[y,x] * total_pin1 * K_out 

			# Manually limit transport if it exceeds the amount of auxin molecules
			# Later on: calculate excess ratio and then use to reduce the value of the vector below
			# 2022.01.17: This situation probably means that paramenters are very unrealistic (it is not likely that in a very small fraction of time all auxin molecules in the cell are transported by PIN1 molecules...). Therefore, increase the auxin / PIN1 ratio ot or reduce k_pin1_transp
			if transported_molecules_total > auxin[y,x]:
				print('warning: transported molecules had to be manually adjusted in cell ' + str(y), str(x))

			# Calculate in an out fluxes for each cell face (molecules / face and cycle)
			# (calculating fluxes is necessary for WTF polarization)

			# Top face: out (fluxes_pin1[0,y,x]); in (fluxes_pin1[1,y,x])
			if y > 0:
				fluxes_pin1[0,y,x] = h * ( auxin[y,x] * pin1[0,y,x] * K_out )
				fluxes_pin1[1,y,x]  = h * ( auxin[y-1,x] * pin1[2,y-1,x] * K_fromT )
			else:
				fluxes_pin1[0,y,x], fluxes_pin1[1,y,x] = 0, 0
			
			# Right face: out (fluxes_pin1[2,y,x]); in (fluxes_pin1[3,y,x])
			if x < tissue_columns - 1:
				fluxes_pin1[2,y,x] = h * ( auxin[y,x] * pin1[1,y,x] * K_out )
				fluxes_pin1[3,y,x]  = h * ( auxin[y,x+1] * pin1[3,y,x+1] * K_fromR )
			else:
				fluxes_pin1[2,y,x], fluxes_pin1[3,y,x] = 0, 0
			
			# Bottom face: out (fluxes_pin1[4,y,x]); in (fluxes_pin1[5,y,x])
			if y < tissue_rows - 1:
				fluxes_pin1[4,y,x] = h * ( auxin[y,x] * pin1[2,y,x] * K_out )
				fluxes_pin1[5,y,x]  = h * ( auxin[y+1,x] * pin1[0,y+1,x] * K_fromB )
			else:
				fluxes_pin1[4,y,x], fluxes_pin1[5,y,x] = 0, 0
			
			# Left face: out (fluxes_pin1[6,y,x]); in (fluxes_pin1[7,y,x])
			if x > 0:
				fluxes_pin1[6,y,x] = h * ( auxin[y,x] * pin1[3,y,x] * K_out )
				fluxes_pin1[7,y,x]  = h * ( auxin[y,x-1] * pin1[1,y,x-1] * K_fromL )
			else:
				fluxes_pin1[6,y,x], fluxes_pin1[7,y,x] = 0, 0
			
			# Calculate vector components to draw line and to orient image with arrow direction

			# Calculate net fluxes (outbound, out - in) to draw vectors
			T_net_flux = fluxes_pin1[0,y,x] - fluxes_pin1[1,y,x]
			R_net_flux = fluxes_pin1[2,y,x] - fluxes_pin1[3,y,x]
			B_net_flux = fluxes_pin1[4,y,x] - fluxes_pin1[5,y,x]
			L_net_flux = fluxes_pin1[6,y,x] - fluxes_pin1[7,y,x]
			# Vector X component
			vector_x = fluxes_pin1[8,y,x] = R_net_flux - L_net_flux
			# Vector Y component
			vector_y = fluxes_pin1[9,y,x] = B_net_flux - T_net_flux

	# Update the auxin concentrations after calculating all the fluxes to avoid polarity effect of looping through numpy array
	# This could go inside auxin_homeostasis()
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			fluxes_pin1_out = fluxes_pin1[0,y,x] + fluxes_pin1[2,y,x] + fluxes_pin1[4,y,x] + fluxes_pin1[6,y,x]
			fluxes_pin1_in = fluxes_pin1[1,y,x] + fluxes_pin1[3,y,x] + fluxes_pin1[5,y,x] + fluxes_pin1[7,y,x]

			auxin[y,x] = auxin[y,x] - fluxes_pin1_out + fluxes_pin1_in


def pin_on_auxin():

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
	cuc = ip.cuc
	Kb = pr.k_pin1_effi_basal
	Kp = pr.k_pin1_effi_pho
	Kcp = pr.pin1_pho_k05
	hc = pr.pin1_pho_hc
	tissue_rows = ip.tissue_rows
	tissue_columns = ip.tissue_columns
	fluxes_pin1 = ip.auxin_fluxes_pin1

	# Calculate efficiencies of each cell:
	# Create array with same shape as tissue to store the PIN1 efficiency value of each cell
	# The PIN1 efficiency of a cell depends on the proportion of unphosphorylated and phosporylated PIN1 molecules, which have efficiencies Kb (basal) and Kp (phospho-activated) 
	
	effs = np.zeros(shape=(tissue_rows, tissue_columns))
	
	'''
	# Older method (linear function with a cap)
	if pr.pin1_polarity == 'dual':
		effs.fill(Kb)
	if pr.pin1_polarity == 'utg_smith2006':

		Kcp = .01

		for y in range(tissue_rows):
			for x in range(tissue_columns):
				K = Kb + cuc[y,x] * Kcp
				if K <= Kp:
					effs[y,x] = K
				else:
					effs[y,x] = Kp
				print(cuc[y,x], effs[y,x])
	'''
	
	if pr.pin1_polarity == 'utg_smith2006':
		for y in range(tissue_rows):
			for x in range(tissue_columns):
				# Calculate proportion of phosphorylated and unphosphorylated PIN molecules in the cell
				fraction_pin1_p = cuc[y,x]**hc / (Kcp**hc + cuc[y,x]**hc)
				fraction_pin1_u = 1 - fraction_pin1_p
				# Compute the efficiency of the whole cell
				effs[y,x] = ( Kb * fraction_pin1_u ) + ( Kp * fraction_pin1_p )
				#print(cuc[y,x], effs[y,x])
	else:
		effs.fill(Kb)
	
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			# Calculate K as K basal + effect of CUC
			# Kp = basal efficiency
			# Kcp = cuc (phosphorilation) strength to increase efficiency
			K_out = effs[y,x]
			if y > 0:
				K_fromT = effs[y-1,x]
			if x < tissue_columns - 1:
				K_fromR = effs[y,x+1]
			if y < tissue_rows - 1:
				K_fromB = effs[y+1,x]
			if x > 0:
				K_fromL = effs[y,x-1]

			# Certain combinations of [PIN1] and K could result in that the sum of auxin molecules effluxed (through all the cell faces) exceeds the num of auxin molecules in the cell. To detect that, first calculate the efflux in the whole cell and check if it is higher than the num of auxin molecules available.
			total_pin1 = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
			transported_molecules_total = h * auxin[y,x] * total_pin1 * K_out 

			# Manually limit transport if it exceeds the amount of auxin molecules
			# Later on: calculate excess ratio and then use to reduce the value of the vector below
			# 2022.01.17: This situation probably means that paramenters are very unrealistic (it is not likely that in a very small fraction of time all auxin molecules in the cell are transported by PIN1 molecules...). Therefore, increase the auxin / PIN1 ratio ot or reduce k_pin1_transp
			if transported_molecules_total > auxin[y,x]:
				print('warning: transported molecules had to be manually adjusted in cell ' + str(y), str(x))

			# Calculate in an out fluxes for each cell face (molecules / face and cycle)
			# (calculating fluxes is necessary for WTF polarization)

			# Top face: out (fluxes_pin1[0,y,x]); in (fluxes_pin1[1,y,x])
			if y > 0:
				fluxes_pin1[0,y,x] = h * ( auxin[y,x] * pin1[0,y,x] * K_out )
				fluxes_pin1[1,y,x]  = h * ( auxin[y-1,x] * pin1[2,y-1,x] * K_fromT )
			else:
				fluxes_pin1[0,y,x], fluxes_pin1[1,y,x] = 0, 0
			
			# Right face: out (fluxes_pin1[2,y,x]); in (fluxes_pin1[3,y,x])
			if x < tissue_columns - 1:
				fluxes_pin1[2,y,x] = h * ( auxin[y,x] * pin1[1,y,x] * K_out )
				fluxes_pin1[3,y,x]  = h * ( auxin[y,x+1] * pin1[3,y,x+1] * K_fromR )
			else:
				fluxes_pin1[2,y,x], fluxes_pin1[3,y,x] = 0, 0
			
			# Bottom face: out (fluxes_pin1[4,y,x]); in (fluxes_pin1[5,y,x])
			if y < tissue_rows - 1:
				fluxes_pin1[4,y,x] = h * ( auxin[y,x] * pin1[2,y,x] * K_out )
				fluxes_pin1[5,y,x]  = h * ( auxin[y+1,x] * pin1[0,y+1,x] * K_fromB )
			else:
				fluxes_pin1[4,y,x], fluxes_pin1[5,y,x] = 0, 0
			
			# Left face: out (fluxes_pin1[6,y,x]); in (fluxes_pin1[7,y,x])
			if x > 0:
				fluxes_pin1[6,y,x] = h * ( auxin[y,x] * pin1[3,y,x] * K_out )
				fluxes_pin1[7,y,x]  = h * ( auxin[y,x-1] * pin1[1,y,x-1] * K_fromL )
			else:
				fluxes_pin1[6,y,x], fluxes_pin1[7,y,x] = 0, 0
			
			# Calculate vector components to draw line and to orient image with arrow direction

			# Calculate net fluxes (outbound, out - in) to draw vectors
			T_net_flux = fluxes_pin1[0,y,x] - fluxes_pin1[1,y,x]
			R_net_flux = fluxes_pin1[2,y,x] - fluxes_pin1[3,y,x]
			B_net_flux = fluxes_pin1[4,y,x] - fluxes_pin1[5,y,x]
			L_net_flux = fluxes_pin1[6,y,x] - fluxes_pin1[7,y,x]
			# Vector X component
			vector_x = fluxes_pin1[8,y,x] = R_net_flux - L_net_flux
			# Vector Y component
			vector_y = fluxes_pin1[9,y,x] = B_net_flux - T_net_flux

			# Angle of net flux is calculated in func_graph module, in functions create_cell_plot() and vector_to_degrees()
			
			'''
			# Calculate sine of angle described by vector
			vector_hyp = math.sqrt(vector_x**2 + vector_y**2)
			if vector_y != 0 and vector_hyp != 0:
				sine = - vector_y / vector_hyp
				# Calculate degres of vector (convert sine to radians and then to degrees)
				fluxes_pin1[10,y,x] = math.degrees(math.asin(sine))
			else:
				fluxes_pin1[10,y,x] = 366.0

			#print(vector_y, vector_hyp, fluxes_pin1[10,y,x])
			'''

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


def auxin_homeostasis_old(iteration, sim_time):
	
	'''
	This function implements all changes in auxin concentration in which there
	is not communication between cells. These processes can be:
	- Basal (global) de novo synthesis
	- Basal (global) turnover
	- Local synthesis (e.g. CUC effect on YUC)
	- Local degradation
	- Local exogenous application, etc
	- Noise: random variation in concentration 

	A' = h * ( k(synth) - A*k(decay) + custom_synth_degr + CUC*k(CY) ) + noise

	Params:
	* sim_time: Iteration of the simulation. This is used for local auxin synth/degr
	
	'''

	# Define simpler aliases
	h = pr.euler_h
	k_auxin_synth = pr.k_auxin_synth
	k_auxin_degr = pr.k_auxin_degr
	k_cuc_auxin_synth = pr.k_cuc_auxin_synth
	k_md_auxin_synth = pr.k_md_auxin_synth
		
	for y in range(ip.tissue_rows):
		for x in range(ip.tissue_columns):
			
			# Simplify var names
			auxin_cell = ip.auxin[y,x]
			cuc_cell = ip.cuc[y,x]
			md_cell = ip.middle_domain[x]
			
			# Local/custom auxin synth/degr and noise...
			current_cell = (y,x)
			custom_synth, custom_degr, noise = 0, 0, 0

			if pr.auxin_custom_synth['value'] > 0:
				if current_cell in pr.auxin_custom_synth['cells'] \
				and sim_time >= pr.auxin_custom_synth['time_interval'][0] \
				and sim_time < pr.auxin_custom_synth['time_interval'][1]:
					custom_synth = pr.auxin_custom_synth['value']
			
			if pr.auxin_custom_degr['value'] > 0:
				if current_cell in pr.auxin_custom_degr['cells'] \
				and sim_time >= pr.auxin_custom_degr['time_interval'][0] \
				and sim_time < pr.auxin_custom_degr['time_interval'][1]:
					custom_degr = pr.auxin_custom_degr['value'] * auxin_cell
			
			# Noise
			if pr.auxin_noise['limit'] > 0 \
			and iteration >= pr.auxin_noise['iteration_interval'][0] \
			and iteration < pr.auxin_noise['iteration_interval'][1]:
				noise = auxin_cell * ( 1 + random.uniform(-pr.auxin_noise['limit'], pr.auxin_noise['limit']) )
			
			# Calculate change in auxin concentration
			auxin_cell_updated = auxin_cell + \
				h * ( \
					k_auxin_synth - \
					k_auxin_degr * auxin_cell + \
					k_cuc_auxin_synth * cuc_cell + \
					k_md_auxin_synth * md_cell + \
					custom_synth - \
					custom_degr \
				) \
				+ noise
			
			if auxin_cell_updated < 0:
				auxin_cell_updated = float(1E-6)
				
			ip.auxin[y,x] = auxin_cell_updated


def auxin_custom_manipulation_old(iteration, sim_time):
	
	'''
	This function implements some changes in auxin:
	- Local exogenous application, etc
	- Noise: random variation in concentration

	Inputs and outputs:
	- Function does not have parameter inputs. It reads pr.template_auxin
	- Function does not return any objects. It writes to pr.template_auxin
	
	Params:
	- sim_time: Iteration of the simulation. This is used for custom (local) auxin synth/degr
	
	'''

	# Define simpler aliases
	h = pr.euler_h
		
	for y in range(ip.tissue_rows):
		for x in range(ip.tissue_columns):
			
			# Simplify var names
			auxin_cell = ip.auxin[y,x]
			cuc_cell = ip.cuc[y,x]
			md_cell = ip.middle_domain[x]

			current_cell = (y,x)
			custom_synth, custom_degr, noise = 0, 0, 0
			
			# Local/custom auxin synth/degr and noise...
			if pr.auxin_custom_synth['value'] > 0:
				if current_cell in pr.auxin_custom_synth['cells'] \
				and sim_time >= pr.auxin_custom_synth['time_interval'][0] \
				and sim_time < pr.auxin_custom_synth['time_interval'][1]:
					custom_synth = pr.auxin_custom_synth['value']
			if pr.auxin_custom_degr['value'] > 0:
				if current_cell in pr.auxin_custom_degr['cells'] \
				and sim_time >= pr.auxin_custom_degr['time_interval'][0] \
				and sim_time < pr.auxin_custom_degr['time_interval'][1]:
					custom_degr = pr.auxin_custom_degr['value'] * auxin_cell
			
			# Noise
			if pr.auxin_noise['limit'] > 0 \
			and iteration >= pr.auxin_noise['iteration_interval'][0] \
			and iteration < pr.auxin_noise['iteration_interval'][1]:
				
				# Convert % to absolute limits for current cell, then take random number within limits 
				noise_lim_cell = auxin_cell * ( pr.auxin_noise['limit'] / 100 )
				noise = random.uniform(-noise_lim_cell, noise_lim_cell)
				
				# This is old system using predefined absolute limits
				#noise = random.uniform(-pr.auxin_noise['limit'], pr.auxin_noise['limit'])
				
				#print(noise)

			# Calculate change in auxin concentration
			auxin_cell_updated = auxin_cell + h * ( custom_synth - custom_degr ) + noise

			# Perfect sources and sinks (overwrites noise and synthesis/degradation)
			if pr.auxin_perfect_sources['active'] == True:
				if current_cell in pr.auxin_perfect_sources['cells'] \
				and sim_time >= pr.auxin_perfect_sources['time_interval'][0] \
				and sim_time < pr.auxin_perfect_sources['time_interval'][1]:
					auxin_cell_updated = pr.auxin_perfect_sources['value']
			if pr.auxin_perfect_sinks['active'] == True:
				if current_cell in pr.auxin_perfect_sinks['cells'] \
				and sim_time >= pr.auxin_perfect_sinks['time_interval'][0] \
				and sim_time < pr.auxin_perfect_sinks['time_interval'][1]:
					auxin_cell_updated = pr.auxin_perfect_sinks['value']
			
			if auxin_cell_updated < 0:
				auxin_cell_updated = float(1E-6)
			



			
			ip.auxin[y,x] = auxin_cell_updated


if __name__ == '__main__':
    pass




