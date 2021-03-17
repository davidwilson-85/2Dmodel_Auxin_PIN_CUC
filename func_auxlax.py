#!/usr/bin/env python

import numpy as np

import params as pr
import inputs as ip


def aux1_expression():

	'''
	* Auxin promotes PIN1 expression.
	* CUC promotes PIN1 expression.
	
	To model this, assume that auxin and CUCs increase PIN1 expression and PIN1 has a constant turnover/decay
	
	P' = h * ( A*P*K(AP) + C*P*K(CP) - P*K(Pdecay) )

	'''

	# Rename parameters
	pin1 = ip.pin1
	auxin = ip.auxin
	cuc = ip.cuc
	k_auxin_pin1 = pr.k_auxin_pin1
	k_cuc_pin1 = pr.k_cuc_pin1
	k_pin1_decay = pr.k_pin1_decay
	h = pr.euler_h
	
	# For each cell, calculate updated PIN1 expression value
	for y in range(auxin.shape[0]):
		for x in range(auxin.shape[1]):
		
			pin1_cell = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
			auxin_cell = auxin[y,x]
			cuc_cell = cuc[y,x]
			
			pin1_cell_updated = pin1_cell + h * ( auxin_cell * k_auxin_pin1 + cuc_cell * k_cuc_pin1 - pin1_cell * k_pin1_decay )

			#pin1_cell_updated = pin1_cell + h * ( auxin_cell * pin1_cell * k_auxin_pin1 + cuc_cell * pin1_cell * k_cuc_pin1 - pin1_cell * k_pin1_decay )
			
			if pin1_cell_updated < 0:
				pin1_cell_updated = 0
			
			pin1_ratio = float(pin1_cell_updated / pin1_cell)
			
			pin1[0,y,x] = pin1[0,y,x] * pin1_ratio
			pin1[1,y,x] = pin1[1,y,x] * pin1_ratio
			pin1[2,y,x] = pin1[2,y,x] * pin1_ratio
			pin1[3,y,x] = pin1[3,y,x] * pin1_ratio



def lax1_expression():

	'''
	Kasprzewska et al 2015: LAX1::GUS is expressed in same places as DR5.

	'''

	# Rename parameters
	pin1 = ip.pin1
	auxin = ip.auxin
	cuc = ip.cuc
	k_auxin_pin1 = pr.k_auxin_pin1
	k_cuc_pin1 = pr.k_cuc_pin1
	k_pin1_decay = pr.k_pin1_decay
	h = pr.euler_h
	
	# For each cell, calculate updated PIN1 expression value
	for y in range(auxin.shape[0]):
		for x in range(auxin.shape[1]):
		
			pin1_cell = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
			auxin_cell = auxin[y,x]
			cuc_cell = cuc[y,x]
			
			pin1_cell_updated = pin1_cell + h * ( auxin_cell * k_auxin_pin1 + cuc_cell * k_cuc_pin1 - pin1_cell * k_pin1_decay )

			#pin1_cell_updated = pin1_cell + h * ( auxin_cell * pin1_cell * k_auxin_pin1 + cuc_cell * pin1_cell * k_cuc_pin1 - pin1_cell * k_pin1_decay )
			
			if pin1_cell_updated < 0:
				pin1_cell_updated = 0
			
			pin1_ratio = float(pin1_cell_updated / pin1_cell)
			
			pin1[0,y,x] = pin1[0,y,x] * pin1_ratio
			pin1[1,y,x] = pin1[1,y,x] * pin1_ratio
			pin1[2,y,x] = pin1[2,y,x] * pin1_ratio
			pin1[3,y,x] = pin1[3,y,x] * pin1_ratio



if __name__ == '__main__':
    pass

