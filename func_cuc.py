#!/usr/bin/env python

import numpy as np


	
def cuc_expression(middle_domain, auxin, cuc, k_md_cuc, k_auxin_cuc, k_cuc_decay):

	#
	# CUC expression is promoted in the middle domain, repressed by auxin, and decays at a constant rate. [Maybe add also turnover]
	#
	# C' = h * ( M*k(MC) - A*k(AC) - C*k(Cdecay) )
	#
	#

	for y in range(auxin.shape[0]):
		for x in range(auxin.shape[1]):
		
			cuc_cell = cuc[y,x]
			auxin_cell = auxin[y,x]
			md_cell = middle_domain[x]
			
			cuc_cell_updated = cuc_cell + h * ( md_cell * k_md_cuc - auxin_cell * k_auxin_cuc - cuc_cell * k_cuc_decay )
			
			if cuc_cell_updated < 0:
				cuc_cell_updated = 0
				
			cuc[y,x] = cuc_cell_updated



if __name__ == '__main__':
    pass


