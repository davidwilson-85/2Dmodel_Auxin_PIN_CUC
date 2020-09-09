#!/usr/bin/env python

import numpy as np
import math



def pin_expression(pin1, auxin, cuc):

	#
	# Auxin promotes PIN1 expression. To model this, assume that auxin increases PIN1 expression and PIN1 has a constant turnover/decay
	#
	# P' = A*P*K(AP) + C*P*K(CP) - P*K(Pdecay)
	#
	#
	
	for y in range(auxin.shape[0]):
		for x in range(auxin.shape[1]):
		
			pin1_cell = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
			auxin_cell = auxin[y,x]
			cuc_cell = cuc[y,x]
			
			k_auxin_pin1 = 0 #0.0001
			k_cuc_pin1 = 0 #0.0001
			k_pin1_decay = 0 # 0.004
			
			pin1_cell_updated = pin1_cell + auxin_cell * pin1_cell * k_auxin_pin1 + cuc_cell * pin1_cell * k_cuc_pin1 - pin1_cell * k_pin1_decay
			
			if pin1_cell_updated < 0:
				pin1_cell_updated = 0
			
			pin1_ratio = float(pin1_cell_updated / pin1_cell)
			
			pin1[0,y,x] = pin1[0,y,x] * pin1_ratio
			pin1[1,y,x] = pin1[1,y,x] * pin1_ratio
			pin1[2,y,x] = pin1[2,y,x] * pin1_ratio
			pin1[3,y,x] = pin1[3,y,x] * pin1_ratio




def pin_utg_smith2006(auxin, pin1, k_UTG, cuc, cuc_threshold_pin1):

	#
	# Auxin affect PIN1 subcellular localization (up-the-gradient model = UTG)
	# UTG: PIN1 accumulates at membrane abutting cells with higher auxin concentration
	# 
	# I use formula from Smith 2006 and Bilsborough 2011:
	#
	#                                    b^A[i]
	# PIN[ij] (potential) = PIN[i] * ---------------
	#                                 SUM[k] b^A[k] 
	#
	
	# Base of exponential function to tweak with UTG responsiveness
	b = k_UTG
	
	
	for y in range(auxin.shape[0]):
		for x in range(auxin.shape[1]):
		
			# Current PIN1 total amount in the cell
			total_pin1 = pin1[0,y,x] + pin1[1,y,x] + pin1[2,y,x] + pin1[3,y,x]
	
			# Correct boundary effect
			# Top
			if y > 0:
				auxin_top = auxin[y-1,x]
			else:
				auxin_top = auxin[y,x]
			
			# Right
			if x < auxin.shape[1] - 1:
				auxin_right = auxin[y,x+1]
			else:
				auxin_right = auxin[y,x]
			
			# Bottom
			if y < auxin.shape[0] - 1:
				auxin_bottom = auxin[y+1,x]
			else:
				auxin_bottom = auxin[y,x]
			
			# Left
			if x > 0:
				auxin_left = auxin[y,x-1]
			else:
				auxin_left = auxin[y,x]
	
			# Calculate normalization factor (eq. denominator)
			norm_factor = b**auxin_top + b**auxin_right + b**auxin_bottom + b**auxin_left
	
			pin1[0,y,x] = total_pin1 * ( b**auxin_top / norm_factor )
			pin1[1,y,x] = total_pin1 * ( b**auxin_right / norm_factor )
			pin1[2,y,x] = total_pin1 * ( b**auxin_bottom / norm_factor )
			pin1[3,y,x] = total_pin1 * ( b**auxin_left / norm_factor )


			#print pin1[0,y,x], pin1[1,y,x], pin1[2,y,x], pin1[3,y,x]

			#print utg_auxinRatioT, utg_auxinRatioR, utg_auxinRatioB, utg_auxinRatioL
			
			# CUC effect on PIN1 polarity
			# For now, simply reverse the values in the X and Y axes. This will be in favour of the auxin gradient.
			if cuc[y,x] > cuc_threshold_pin1:
				
				pin1[0,y,x] = pin1[2,y,x]
				pin1[1,y,x] = pin1[3,y,x]
				pin1[2,y,x] = pin1[0,y,x]
				pin1[3,y,x] = pin1[1,y,x]





def pin_utg_ratio(auxin, pin1, k_UTG, cuc, cuc_threshold_pin1):

	#
	# Auxin affect PIN1 subcellular localization (up-the-gradient model = UTG)
	# UTG: PIN1 accumulates at membrane abutting cells with higher auxin concentration
	# 
	#

	# Calculate updated PIN1 allocation
	for y in range(auxin.shape[0]):
		for x in range(auxin.shape[1]):

			# Top	
			if y > 0:
				utg_auxinRatioT = math.sqrt(auxin[y-1,x] / auxin[y,x])
			else:
				utg_auxinRatioT = 1.0
			
			# Right
			if x < auxin.shape[1] - 1:
				utg_auxinRatioR = math.sqrt(auxin[y,x+1] / auxin[y,x])
			else:
				utg_auxinRatioR = 1.0
				
			# Bottom
			if y < auxin.shape[0] - 1:
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
			
			# CUC effect on PIN1 polarity
			# For now, simply reverse the values in the X and Y axes. This will be in favour of the auxin gradient.
			if cuc[y,x] > cuc_threshold_pin1:
				
				pin1[0,y,x] = pin1[2,y,x]
				pin1[1,y,x] = pin1[3,y,x]
				pin1[2,y,x] = pin1[0,y,x]
				pin1[3,y,x] = pin1[1,y,x]





def cuc_on_pin_polarity():

	# 
	# CUC genes affect PIN1 subcellular localization
	# 
	# It is not clear how. Test different hypotheses (WTF, reversal, non-polar, DTCG, UTG dampening, etc...) 
	# 
	# 

	pass



if __name__ == '__main__':
    pass

