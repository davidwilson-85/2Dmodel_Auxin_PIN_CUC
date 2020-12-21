#!/usr/bin/env python

import numpy as np
import params as pr
import inputs as ip

	
def cuc_expression():

	#
	# CUC expression is promoted in the middle domain, repressed by auxin, and decays at a constant rate.
	#
	# C' = h * ( M*k(MC) - A*k(AC) - C*k(Cdecay) )
	#
	#

	# Rename parameters
	h = pr.euler_h
	k_md_cuc = pr.k_md_cuc
	k_auxin_cuc = pr.k_auxin_cuc
	k_cuc_decay = pr.k_cuc_decay

	for y in range(ip.tissue_rows):
		for x in range(ip.tissue_columns):
		
			cuc_cell = ip.cuc[y,x]
			auxin_cell = ip.auxin[y,x]
			md_cell = ip.middle_domain[y,x]
			#md_cell = ip.middle_domain[y,x] For 1-D midlle domain
			
			cuc_cell_updated = cuc_cell + h * ( md_cell * k_md_cuc - auxin_cell * k_auxin_cuc - cuc_cell * k_cuc_decay )
			
			if cuc_cell_updated < 0:
				cuc_cell_updated = 0
				
			ip.cuc[y,x] = cuc_cell_updated



if __name__ == '__main__':
    pass


