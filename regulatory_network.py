#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint

import params as pr
import inputs as ip


def model_regulatory_network(init_values, t):

    '''
    This model is solved by ODEint.
    It only defines changes that can be described as differential equations and that occur within each individual cell.
    Changes not included: those that involve movement of auxin between cells, PIN1 polarization 
    '''

    A, C, Pt, Pr, Pb, Pl, MD = init_values
    
    k_AS = pr.k_auxin_synth
    k_CA = pr.k_cuc_auxin_synth
    k_MDA = pr.k_md_auxin_synth
    k_AT = pr.k_auxin_degr

    k_CS = pr.k_cuc_synth
    k_MDC = pr.k_md_cuc
    k_AC = pr.k_auxin_cuc
    k_CT = pr.k_cuc_decay
    
    k_PS = pr.k_pin1_synth
    k_AP = pr.k_auxin_pin1
    k_CP = pr.k_cuc_pin1
    k_PT = pr.k_pin1_decay

    dA_dt = k_AS + C * k_CA + MD * k_MDA - A * k_AT
    #dA_dt = k_AS + .0001 * A**2 + C * k_CA + MD * k_MDA - A * k_AT
    dC_dt = k_CS + MD * k_MDC - A * C * k_AC - C * k_CT
    dPt_dt = k_PS + A * k_AP + C * k_CP - Pt * k_PT
    dPr_dt = k_PS + A * k_AP + C * k_CP - Pr * k_PT
    dPb_dt = k_PS + A * k_AP + C * k_CP - Pb * k_PT
    dPl_dt = k_PS + A * k_AP + C * k_CP - Pl * k_PT
    dMD_dt = 0

    return [dA_dt, dC_dt, dPt_dt, dPr_dt, dPb_dt, dPl_dt, dMD_dt]


def solve_model():
    
    # Solve in a cell-by-cell basis
    for y in range(ip.tissue_rows):
        for x in range(ip.tissue_columns):
            # Gather initial values for ODEint
            model_init_values = [
				ip.auxin[y,x],
				ip.cuc[y,x],
				ip.pin1[0,y,x],
				ip.pin1[1,y,x],
				ip.pin1[2,y,x],
				ip.pin1[3,y,x],
				ip.middle_domain[x]
			]

			# Solve
            cell_solution = odeint(model_regulatory_network, model_init_values, np.linspace(0, pr.euler_h, 2))
            
            # Update current cell in data arrays with solution output
            ip.auxin[y,x] = cell_solution[-1,0]
            ip.cuc[y,x] = cell_solution[-1,1]
            ip.pin1[0,y,x] = cell_solution[-1,2]
            ip.pin1[1,y,x] = cell_solution[-1,3]
            ip.pin1[2,y,x] = cell_solution[-1,4]
            ip.pin1[3,y,x] = cell_solution[-1,5]


def model_v2(init_values, t):

    k_AD = 0.05
    k_AS = 0
    k_AT = 0

    k_CS = 0.01
    k_CT = 0
    k_AC = 0.01
    k_CA = 0

    k_WTF1 = 1
    k_WTF2 = 0.005

    k_PA = 0.0005

    # Auxin in current cell and its neighbours
    A = init_values[0]
    A_T, A_R, A_B, A_L = init_values[1], init_values[2], init_values[3], init_values[4]
    C = init_values[5]
    P_toT, P_toR, P_toB, P_toL = init_values[6], init_values[7], init_values[8], init_values[9]
    P_T_toB, P_R_toL, P_B_toT, P_L_toR = init_values[10], init_values[11], init_values[12], init_values[13]
    
    der_A = k_AS \
          + C * k_CA \
          - k_AD * (4 * A - A_T - A_R - A_B - A_L) \
          - k_PA * ( A * (P_toT + P_toR + P_toB + P_toL) - A_T * P_T_toB - A_R * P_R_toL - A_B + P_B_toT - A_L * P_L_toR ) \
          - A * k_AT

    der_A_flux_T = k_AD * (A_T - A)
    der_A_flux_R = k_AD * (A_R - A)
    der_A_flux_B = k_AD * (A_B - A)
    der_A_flux_L = k_AD * (A_L - A)

    der_C = k_CS - C * k_CT - A * k_AC

    # PIN1 WTF
    net_flux_toT = k_AD * (A - A_T) + k_PA * (P_toT * A - P_T_toB * A_T)
    der_P_toT = k_WTF1 * net_flux_toT - k_WTF2 * P_toT
    net_flux_toR = k_AD * (A - A_R) + k_PA * (P_toR * A - P_R_toL * A_R)
    der_P_toR = k_WTF1 * net_flux_toR - k_WTF2 * P_toR
    net_flux_toB = k_AD * (A - A_B) + k_PA * (P_toB * A - P_B_toT * A_B)
    der_P_toB = k_WTF1 * net_flux_toB - k_WTF2 * P_toB
    net_flux_toL = k_AD * (A - A_L) + k_PA * (P_toL * A - P_L_toR * A_L)
    der_P_toL = k_WTF1 * net_flux_toL - k_WTF2 * P_toL

    der_P_T_toB, der_P_R_toL, der_P_B_toT, der_P_L_toR = 0, 0, 0, 0
    
    return [
        der_A,
        der_A_flux_T,
        der_A_flux_R,
        der_A_flux_B,
        der_A_flux_L,
        der_C,
        der_P_toT,
        der_P_toR,
        der_P_toB,
        der_P_toL,
        der_P_T_toB,
        der_P_R_toL,
        der_P_B_toT,
        der_P_L_toR
    ]

