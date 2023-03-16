#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model(a,t):
    k_synth = 0.05
    k_diff = 0.15
    k_degr = 0.1
    a = a.reshape(2,2)
    da_dt = [0,0,0,0]
    # Top left cell
    da_dt[0] = k_diff * (a[0,1] + a[1,0] - 2 * a[0,0]) + k_synth
    # Top right
    da_dt[1] = k_diff * (a[0,0] + a[1,1] - 2 * a[0,1])
    # Bottom left
    da_dt[2] = k_diff * (a[0,0] + a[1,1] - 2 * a[1,0]) - k_degr * a[1,0]
    # Bottom right
    da_dt[3] = k_diff * (a[0,1] + a[1,0] - 2 * a[1,1])
    #a2 = np.array(da_dt).flatten()
    #return a2
    return da_dt

###########

t = np.linspace(0, 1000, 100)
a_t0 = np.zeros(4)
sol = odeint(model, a_t0, t)

plt.plot(t, sol[: , 0], label='0,0')
plt.plot(t, sol[: , 1], label='0,1')
plt.plot(t, sol[: , 2], label='1,0')
plt.plot(t, sol[: , 3], label='1,1')
plt.xlabel('simulation time (AU)')
plt.ylabel('[auxin] (AU)')
plt.legend()
plt.savefig('odeint/test.png', bbox_inches='tight')

print(sol[: , :][-1])

###########

'''
Sequential approach. I run the solver for a subinterval (e.g. t0 to t5) and store the result. This result is used as initial values for a subsequent run of the solver (e.g. t5 to t10), and so on.
'''

#timepoints = [0,50,100,150,200,250,300,350,400,450]

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

