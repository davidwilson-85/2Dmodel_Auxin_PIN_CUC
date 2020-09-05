#!/usr/bin/env python

import numpy as np
import math



def auxin_on_pin_expression():

	#
	# Auxin promotes PIN1 expression. To model this, assume that auxin increases PIN1 expression and PIN1 has a constant turnover/decay
	#
	#

	pass


def cuc_on_pin_expression():

	#
	# The same as for auxin
	#
	#

	pass


def auxin_on_pin_polarity(auxin, pin1, k_pin1_UTGresponsiveness, tissue_rows, tissue_columns):

	#
	# Auxin promotes PIN1 exppression
	# It also affect its subcellular localization (up-the-gradient model = UTG)
	# UTG: PIN1 accumulates at membrane abutting cells with higher auxin concentration 
	# 
	#

	# Calculate updated PIN1 allocation
	for y in range(tissue_rows):
		for x in range(tissue_columns):

			# Top	
			if y > 0:
				utg_auxinRatioT = math.sqrt(auxin[y-1,x] / auxin[y,x])
			else:
				utg_auxinRatioT = 1.0
			
			# Right
			if x < tissue_columns - 1:
				utg_auxinRatioR = math.sqrt(auxin[y,x+1] / auxin[y,x])
			else:
				utg_auxinRatioR = 1.0
				
			# Bottom
			if y < tissue_rows - 1:
				utg_auxinRatioB = math.sqrt(auxin[y+1,x] / auxin[y,x])
			else:
				utg_auxinRatioB = 1.0
				
			# Left
			if x > 0:
				utg_auxinRatioL = math.sqrt(auxin[y,x-1] / auxin[y,x])
			else:
				utg_auxinRatioL = 1.0
			
			# Calculate updated normalized PIN1 concentration values per cell face
			total_pin1 = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]

			updated_pin1_T = pin1[0,y,x] * utg_auxinRatioT
			updated_pin1_R = pin1[1,y,x] * utg_auxinRatioR
			updated_pin1_B = pin1[2,y,x] * utg_auxinRatioB
			updated_pin1_L = pin1[3,y,x] * utg_auxinRatioL

			total_pin1_updated = updated_pin1_T + updated_pin1_R + updated_pin1_B + updated_pin1_L

			normalization_ratio = total_pin1 / total_pin1_updated

			# Normalize (force the sum of PIN1 in all faces in the cell to remain unchanged)
			pin1[0,y,x] = updated_pin1_T * normalization_ratio
			pin1[1,y,x] = updated_pin1_R * normalization_ratio
			pin1[2,y,x] = updated_pin1_B * normalization_ratio
			pin1[3,y,x] = updated_pin1_L * normalization_ratio

			#pin1[0,0,5] = 9

			#print pin1[0,y,x], pin1[1,y,x], pin1[2,y,x], pin1[3,y,x]

			#print utg_auxinRatioT, utg_auxinRatioR, utg_auxinRatioB, utg_auxinRatioL





def cuc_on_pin_polarity():

	# 
	# CUC genes affect PIN1 subcellular localization
	# 
	# It is not clear how. Test different hypotheses (WTF, reversal, non-polar, DTCG, UTG dampening, etc...) 
	# 
	# 

	pass


