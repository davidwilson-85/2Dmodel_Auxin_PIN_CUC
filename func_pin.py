#!/usr/bin/env python

import importlib, sys, math
import numpy as np

#import params as pr
pr = importlib.import_module(sys.argv[1].split('.')[0], package=None)
import inputs as ip

from scipy.integrate import odeint


def pin_expression():

	'''
	* Auxin promotes PIN1 expression.
	* CUC promotes PIN1 expression.
	
	To model this, assume that auxin and CUCs increase PIN1 expression and PIN1 has a constant turnover/decay
	
	P' = h * ( A*K(AP) + C*K(CP) - P*K(Pdecay) )

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


def pin_polarity():
	
	"""
	This function determines, for each cell in the tissue, the PIN1 polarity mode.
	Then it calls the corresponding function [e.g. pin_utg_smith2006()] that will update the PIN polarity of a single cell
	
	CUC genes affect PIN1 subcellular localization. It is not clear how. Test different
	hypotheses (WTF, reversal, non-polar, DTCG, UTG dampening, etc...)
	
	"""

	for y in range(ip.auxin.shape[0]):
		for x in range(ip.auxin.shape[1]):

			if pr.pin1_polarity == 'utg_smith2006':
				pin_utg_smith2006(y, x, ip.auxin, ip.pin1[:,y,x], pr.k_UTG)

			if pr.pin1_polarity == 'wtf_abley2016':
				pin_wtf_abley2016_odeint_solver(y, x)
			
			if pr.pin1_polarity == 'dual':
				
				if int(ip.cuc[y,x]) >= pr.cuc_threshold_pin1:
					pin_wtf_abley2016_odeint_solver(y, x)
				else:
					pin_utg_smith2006(y, x, ip.auxin, ip.pin1[:,y,x], pr.k_UTG)


def pin_utg_smith2006(y, x, auxin, pin1, k_UTG):

	'''
	Auxin affects PIN1 subcellular localization (up-the-gradient model = UTG)
	UTG: PIN1 tends to accumulates at membrane abutting neighbour cell with highest auxin concentration
	 
	I use formula from Smith 2006 (also used in Bilsborough 2011):
	
	                                   b^A[i]
	PIN[ij] (potential) = PIN[i] * ---------------
	                                SUM[k] b^A[k] 
	
	Current problem: As it is now, function does not take into account current PIN1 distribution, so it erases any initial state defined in the template
	
	Parameter pin1 here does not refer to whole tissue pin1 array, but to cell y,x.

	'''
	
	# Base of exponential function to regulate UTG responsiveness
	b = k_UTG

	# Regulation of k_UTG by m-l axis (suppression of polarity convergences by abaxial or adaxial identity, like KAN1 an KAN2)
	# Map range of MD values to (1 - k_UTG) and overwrite b accorging to the position of the cell in the m-l axis
	# slope = (output_end - output_start) / (input_end - input_start)
	# output = output_start + slope * (input - input_start)
	if pr.md_on_pin1_UTG == True:
		slope = (b - 1) / np.amax(ip.middle_domain) - np.amin(ip.middle_domain)
		b = 1 + slope * (ip.middle_domain[x] - np.amin(ip.middle_domain))
	
	# Current PIN1 total amount in the cell
	total_pin1 = pin1[0] + pin1[1] + pin1[2] + pin1[3]

	# Calculate auxin in neighbours (correcting in boundary cells)
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

	# Calculate PIN1 allocation at each cell face
	pin1[0] = total_pin1 * ( b**auxin_top / norm_factor )
	pin1[1] = total_pin1 * ( b**auxin_right / norm_factor )
	pin1[2] = total_pin1 * ( b**auxin_bottom / norm_factor )
	pin1[3] = total_pin1 * ( b**auxin_left / norm_factor )

	#print pin1[0,y,x], pin1[1,y,x], pin1[2,y,x], pin1[3,y,x]	
	#print utg_auxinRatioT, utg_auxinRatioR, utg_auxinRatioB, utg_auxinRatioL


def pin_wtf_abley2016_odeint_solver(y, x):

	# Calculate net flux at each cell face (out = positive; in = negative)
	# To express it as molecules / time step, I divide by the step size (Euler h)
	h = pr.euler_h
	fd = ip.auxin_fluxes_diffusion
	fp = ip.auxin_fluxes_pin1
	net_flux_t = ( fd[0,y,x] - fd[1,y,x] + fp[0,y,x] - fp[1,y,x] ) / h
	net_flux_r = ( fd[2,y,x] - fd[3,y,x] + fp[2,y,x] - fp[3,y,x] ) / h
	net_flux_b = ( fd[4,y,x] - fd[5,y,x] + fp[4,y,x] - fp[5,y,x] ) / h
	net_flux_l = ( fd[6,y,x] - fd[7,y,x] + fp[6,y,x] - fp[7,y,x] ) / h
	# If net outflux is negative, there is no effect on PIN1 allocation to membrane, so flux is considered as if it was 0
	if net_flux_t < 0: net_flux_t = 0
	if net_flux_r < 0: net_flux_r = 0
	if net_flux_b < 0: net_flux_b = 0
	if net_flux_l < 0: net_flux_l = 0

	# Gather initial values for ODEint in a list
	model_init_values = [
		ip.pin1[0,y,x],
		ip.pin1[1,y,x],
		ip.pin1[2,y,x],
		ip.pin1[3,y,x],
		net_flux_t,
		net_flux_r,
		net_flux_b,
		net_flux_l
	]

	# Solve
	cell_solution = odeint(pin_wtf_abley2016_odeint_model, model_init_values, np.linspace(0, pr.euler_h, 2))
	
	# Update current cell in PIN array with solution output
	# Abley2016: If [PIN1](ij) reaches a threshold [], no more PIN1 can be allocated to membrane ij
	pm = pr.k_WTF_pin1_max
	ip.pin1[0,y,x] = cell_solution[-1,0] if cell_solution[-1,0] <= pm else pm
	ip.pin1[1,y,x] = cell_solution[-1,1] if cell_solution[-1,1] <= pm else pm
	ip.pin1[2,y,x] = cell_solution[-1,2] if cell_solution[-1,2] <= pm else pm
	ip.pin1[3,y,x] = cell_solution[-1,3] if cell_solution[-1,3] <= pm else pm


def pin_wtf_abley2016_odeint_model(init_values, t):

	a = pr.k_WTF_a
	b = pr.k_WTF_b

	Pt, Pr, Pb, Pl, Ft, Fr, Fb, Fl = init_values

	# Linear effect of flux on PIN1
	dPt_dt = a * Ft - b * Pt
	dPr_dt = a * Fr - b * Pr
	dPb_dt = a * Fb - b * Pb
	dPl_dt = a * Fl - b * Pl
	dFt_dt, dFr_dt, dFb_dt, dFl_dt = 0, 0, 0, 0

	return [dPt_dt, dPr_dt, dPb_dt, dPl_dt, dFt_dt, dFr_dt, dFb_dt, dFl_dt]


def pin_wtf_p(y, x, auxin_fluxes, pin1, k_WTF):

	'''
	Calculation similar to UTG (Smith2006) but using passive (p) net flux to allocate PIN1 instead of nighbour auxin concentrations.
	Try using the net flux VS only the efflux
	
	Parameter pin1 here does not refer to whole tissue pin1 array, but to cell y,x
	
	'''

	# Coefficient to tweak WTF responsiveness
	b = k_WTF
		
	# Current PIN1 total amount in the cell
	total_pin1 = pin1[0] + pin1[1] + pin1[2] + pin1[3]

	# Calculate net flux at each cell face
	netflux_top = auxin_fluxes[0,y,x] - auxin_fluxes[1,y,x]
	netflux_right = auxin_fluxes[2,y,x] - auxin_fluxes[3,y,x]
	netflux_bottom = auxin_fluxes[4,y,x] - auxin_fluxes[5,y,x]
	netflux_left = auxin_fluxes[6,y,x] - auxin_fluxes[7,y,x]
	
	# Calculate normalization factor (eq. denominator)
	norm_factor = b**netflux_top + b**netflux_right + b**netflux_bottom + b**netflux_left
	
	# Calculate new PIN amount at each cell face
	pin1[0] = total_pin1 * ( b**netflux_top / norm_factor )
	pin1[1] = total_pin1 * ( b**netflux_right / norm_factor )
	pin1[2] = total_pin1 * ( b**netflux_bottom / norm_factor )
	pin1[3] = total_pin1 * ( b**netflux_left / norm_factor )


def pin_utg_ratio(auxin, pin1, k_UTG, cuc, cuc_threshold_pin1):

	#
	# Auxin affect PIN1 subcellular localization (up-the-gradient model = UTG)
	# UTG: PIN1 accumulates at membrane abutting cells with higher auxin concentration
	# 
	# WARNING: This function has not been adapted yet to be called from pin_polarity()
	# 

	print('WARNING: This function has not been adapted yet to be called from pin_polarity()')

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


def pin_utg_smith2006_mod(y, x, auxin, pin1, k_UTG):

	'''
	Auxin affects PIN1 subcellular localization (up-the-gradient model = UTG)
	UTG: PIN1 accumulates at membrane abutting cells with higher auxin concentration
	 
	I use formula from Smith 2006 (also used in Bilsborough 2011):
	
	                                   b^A[i]
	PIN[ij] (potential) = PIN[i] * ---------------
	                                SUM[k] b^A[k] 
	
	Current problem: As it is now, function does not take into account current PIN1 distribution, so it erases any initial state defined in the template
	
	Parameter pin1 here does not refer to whole tissue pin1 array, but to cell y,x.
	'''
	cuc = ip.cuc
	
	# Base of exponential function to tweak with UTG responsiveness
	b = k_UTG
	
	# Current PIN1 total amount in the cell
	total_pin1 = pin1[0] + pin1[1] + pin1[2] + pin1[3]

	# Calculate auxin in neighbours (correcting in boundary cells)
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

	# Calculate PIN1 alocation at each cell face
	pin1[0] = total_pin1 * ( b**auxin_top / norm_factor )
	pin1[1] = total_pin1 * ( b**auxin_right / norm_factor )
	pin1[2] = total_pin1 * ( b**auxin_bottom / norm_factor )
	pin1[3] = total_pin1 * ( b**auxin_left / norm_factor )

	if cuc[y,x] > pr.cuc_threshold_pin1:
		pin1[0] = 9
		pin1[1] = 9
		pin1[2] = 9
		pin1[3] = 9

	#print pin1[0,y,x], pin1[1,y,x], pin1[2,y,x], pin1[3,y,x]	
	#print utg_auxinRatioT, utg_auxinRatioR, utg_auxinRatioB, utg_auxinRatioL


def pin_wtf_abley2016_old(y, x):

	'''
	Based on Abley et al 2016 (Coen lab). Flux = diffusion + PIN1 transport + import.
	For now I implement only diffusion + PIN1 transport.

	Eqs. from on Abley et al 2016:

	dPIN(ij) / dt = 
		if flux(i->j) >= 0: ka * flux(i->j) - kb * PIN(ij)
		if flux(i->j) <  0: - kb * PIN(ij)

	flux(i->j) = D * [ A(i) - A(j) ] + T * [ (PIN(ij) * A(i) - PIN(ji) * A(j) ]

	'''

	# Aliases var names
	h = pr.euler_h
	a = pr.k_WTF_a
	b = pr.k_WTF_b
	wtf_pin1_max = pr.k_WTF_pin1_max
	pin1 = ip.pin1
	flux_diff = ip.auxin_fluxes_diffusion
	flux_pin1 = ip.auxin_fluxes_pin1
	
	# Calculate net flux at each cell face (out = positive; in = negative)
	# To express it as molecules / time step, I divide by the step size (Euler h)
	net_flux_t = ( flux_diff[0,y,x] - flux_diff[1,y,x] + flux_pin1[0,y,x] - flux_pin1[1,y,x] ) / h
	net_flux_r = ( flux_diff[2,y,x] - flux_diff[3,y,x] + flux_pin1[2,y,x] - flux_pin1[3,y,x] ) / h
	net_flux_b = ( flux_diff[4,y,x] - flux_diff[5,y,x] + flux_pin1[4,y,x] - flux_pin1[5,y,x] ) / h
	net_flux_l = ( flux_diff[6,y,x] - flux_diff[7,y,x] + flux_pin1[6,y,x] - flux_pin1[7,y,x] ) / h

	# Calculate new PIN amount at each cell face
	
	# If net outflux is negative, there is no effect on PIN1 allocation to membrane
	if net_flux_t < 0: net_flux_t = 0
	if net_flux_r < 0: net_flux_r = 0
	if net_flux_b < 0: net_flux_b = 0
	if net_flux_l < 0: net_flux_l = 0

	# Linear effect of flux on PIN1
	pin1[0,y,x] = pin1[0,y,x] + h * ( (a * net_flux_t) - (b * pin1[0,y,x]) )
	pin1[1,y,x] = pin1[1,y,x] + h * ( (a * net_flux_r) - (b * pin1[1,y,x]) )
	pin1[2,y,x] = pin1[2,y,x] + h * ( (a * net_flux_b) - (b * pin1[2,y,x]) )
	pin1[3,y,x] = pin1[3,y,x] + h * ( (a * net_flux_l) - (b * pin1[3,y,x]) )

	# Quadratic effect of flux on PIN1
	#pin1[0,y,x] = pin1[0,y,x] + (net_flux_t**a) - (b * pin1[0,y,x])
	#pin1[1,y,x] = pin1[1,y,x] + (net_flux_r**a) - (b * pin1[1,y,x])
	#pin1[2,y,x] = pin1[2,y,x] + (net_flux_b**a) - (b * pin1[2,y,x])
	#pin1[3,y,x] = pin1[3,y,x] + (net_flux_l**a) - (b * pin1[3,y,x])

	# Abley2016: If [PIN1](ij) reaches a threshold [], no more PIN1 can be allocated to membrane ij
	if pin1[0,y,x] >= wtf_pin1_max: pin1[0,y,x] = wtf_pin1_max
	if pin1[1,y,x] >= wtf_pin1_max: pin1[1,y,x] = wtf_pin1_max
	if pin1[2,y,x] >= wtf_pin1_max: pin1[2,y,x] = wtf_pin1_max
	if pin1[3,y,x] >= wtf_pin1_max: pin1[3,y,x] = wtf_pin1_max




if __name__ == '__main__':
    pass

