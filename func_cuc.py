#!/usr/bin/env python

import importlib, sys, random
import numpy as np
#import params as pr
pr = importlib.import_module(sys.argv[1].split('.')[0], package=None)
import inputs as ip


def cuc_custom_manipulation(iteration):
	
	'''
	This function implements some changes in cuc:
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
			cuc_cell = ip.cuc[y,x]

			custom_synth, custom_degr, noise = 0, 0, 0
			
			# Noise
			if pr.cuc_noise['limit'] > 0 \
			and iteration >= pr.cuc_noise['iteration_interval'][0] \
			and iteration < pr.cuc_noise['iteration_interval'][1]:
				# Convert % to absolute limits for current cell, then take random number within limits 
				noise_lim_cell = cuc_cell * ( pr.cuc_noise['limit'] / 100 )
				noise = random.uniform(-noise_lim_cell, noise_lim_cell)
				#print(noise)

			# Calculate change in CUC concentration and correct negative values if necessary
			cuc_cell_updated = cuc_cell + noise
			#if math.isnan(cuc_cell_updated) == True: print(auxin_cell_updated)
			if cuc_cell_updated < 0:
				cuc_cell_updated = float(1E-6)
			
			# Update CUC array
			ip.cuc[y,x] = cuc_cell_updated
	
def cuc_expression():

	'''
	DEPRECATED. THIS FUNCTIONALITY IS NOW INTEGRATED IN regulatory_network.py.

	CUC expression is produced at a constant rate, (promoted in the middle domain), repressed by side (adaxial and abaxial) domain, repressed by auxin, and decays at a constant rate.

	C' = h * [ k(C) - S*k(SC) - A*k(AC) - C*k(Cdecay) ]
	
	'''

	# Rename parameters
	h = pr.euler_h
	k_cuc = pr.k_cuc
	k_md_cuc = pr.k_md_cuc
	#k_adab_cuc = pr.k_adab_cuc
	k_auxin_cuc = pr.k_auxin_cuc
	k_cuc_decay = pr.k_cuc_decay

	for y in range(ip.tissue_rows):
		for x in range(ip.tissue_columns):
		
			cuc_cell = ip.cuc[y,x]
			auxin_cell = ip.auxin[y,x]
			#adab_cell = ip.adab_domain[x]
			md_cell = ip.middle_domain[x]
			#md_cell = ip.middle_domain[y,x] For 1-D middle domain
			
			cuc_cell_updated = cuc_cell + h * ( k_cuc + md_cell * k_md_cuc - auxin_cell * k_auxin_cuc - cuc_cell * k_cuc_decay )
			
			if cuc_cell_updated < 0:
				cuc_cell_updated = 0
				
			ip.cuc[y,x] = cuc_cell_updated



if __name__ == '__main__':
    pass


