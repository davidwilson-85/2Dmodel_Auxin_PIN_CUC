#!/usr/bin/env python

import numpy as np


	
def cuc_expression(middle_domain, auxin, cuc):

	#
	# CUC expression is promoted in the middle domain, repressed by auxin, and decays at a constant rate.
	#
	# C' = M*K(MC) - A*K(AC) - C*K(Cdecay)
	#
	#

	for y in range(auxin.shape[0]):
		for x in range(auxin.shape[1]):
		
			cuc_cell = cuc[y,x]
			auxin_cell = auxin[y,x]
			md_cell = cuc[y,x]
			
			k_md_cuc = 0.01
			k_auxin_cuc = 0
			k_cuc_decay = 0
			
			cuc_cell_updated = cuc_cell + md_cell * k_md_cuc - auxin_cell * k_auxin_cuc - cuc_cell * k_cuc_decay
			
			if cuc_cell_updated < 0:
				cuc_cell_updated = 0
				
			cuc[y,x] = cuc_cell_updated





