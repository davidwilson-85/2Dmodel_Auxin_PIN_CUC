#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
import params_v3 as pr


def model_regulatory_network(init_values, t):

    '''
    This model is solved by ODEint.
    It only defines changes that can be described as differential equations and that occur within each individual cell.
    Changes not included: those resulting from movement of auxin between cells, PIN1 polarization 
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
    
    k_PS = 0.01
    k_AP = pr.k_auxin_pin1
    k_CP = pr.k_cuc_pin1
    k_PT = pr.k_pin1_decay

    dA_dt = k_AS + C * k_CA + MD * k_MDA - A * k_AT
    dC_dt = k_CS + MD * k_MDC - A * k_AC - C * k_CT
    dPt_dt = k_PS + A * k_AP * C * k_CP - Pt * k_PT
    dPr_dt = k_PS + A * k_AP * C * k_CP - Pr * k_PT
    dPb_dt = k_PS + A * k_AP * C * k_CP - Pb * k_PT
    dPl_dt = k_PS + A * k_AP * C * k_CP - Pl * k_PT
    dMD_dt = 0

    return [dA_dt, dC_dt, dPt_dt, dPr_dt, dPb_dt, dPl_dt, dMD_dt]


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

###########

'''
time_points = np.linspace(0, 1, 2)
solution = odeint(model, cell_A, time_points)

t = np.linspace(0, 1000, 100)
a_t0 = np.zeros(4)
sol = odeint(model, a_t0, t)

print(sol[: , :][-1])


t_ini = 0
t_fin = 499
nbr_points = 500000

timepoints = np.linspace(t_ini, t_fin, nbr_points) # From A to B including X points
interval_length = (t_fin - t_ini) / (nbr_points - 1)
a_seq = np.zeros(4)
stacked_result = []

for i in timepoints:
    t = np.linspace(i, i+interval_length, 2)
    sol_seq = odeint(model, a_seq, t)
    a_seq = sol_seq[1 , :].flatten()
    results = sol_seq[: , :]
    stacked_result.append(results)

print(stacked_result[-1][-1])
'''
