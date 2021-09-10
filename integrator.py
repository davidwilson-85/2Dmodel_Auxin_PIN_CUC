#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint

def model(init_values, t):

    k_degr = 0.02
    k_AD = 0.05
    k_AS = 0
    k_AT = 0

    k_AC = 0
    k_CA = 0

    # Auxin in current cell and its neighbours
    A_i = init_values[0]
    A_T = init_values[1]
    A_R = init_values[2]
    A_B = init_values[3]
    A_L = init_values[4]

    # CUC
    C = 0
    # PIN1
    P = 0
    
    der_A_i = k_AS + C * k_CA + k_AD * (A_T + A_R + A_B + A_L - 4 * A_i) - A_i * k_AT

    der_A_flux_T = k_AD * (A_T - A_i)
    der_A_flux_R = k_AD * (A_R - A_i)
    der_A_flux_B = k_AD * (A_B - A_i)
    der_A_flux_L = k_AD * (A_L - A_i)

    der_C = k_CS - C * k_CT - A * k_AC
    
    return [
        der_A_i,
        der_A_flux_T,
        der_A_flux_R,
        der_A_flux_B,
        der_A_flux_L
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
