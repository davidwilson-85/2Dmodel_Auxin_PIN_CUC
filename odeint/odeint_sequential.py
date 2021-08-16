#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model(a,t):
    k_sy = 0.05
    k_di = 0.15
    k_de = 0.1
    a = a.reshape(2,2)
    da_dt = [0,0,0,0]
    da_dt[0] = k_di * (a[0,1] + a[1,0] - 2 * a[0,0]) + k_sy
    da_dt[1] = k_di * (a[0,0] + a[1,1] - 2 * a[0,1])
    da_dt[2] = k_di * (a[0,0] + a[1,1] - 2 * a[1,0]) - k_de * a[1,0]
    da_dt[3] = k_di * (a[0,1] + a[1,0] - 2 * a[1,1])
    #a2 = np.array(da_dt).flatten()
    #return a2
    return da_dt

###########
'''
t = np.linspace(0, 500, 500)
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

print(sol[: , :])
'''
###########

#timepoints = [0,50,100,150,200,250,300,350,400,450]
timepoints = np.linspace(0, 499, 500)
a_seq = np.zeros(4)
stacked_result = []

for i in timepoints:
    t = np.linspace(i, i+50, 5)
    sol_seq = odeint(model, a_seq, t)
    print(sol_seq)
    a_seq = sol_seq[4 , :].flatten()
    results = sol_seq[: , :]
    stacked_result.append(results)

print(stacked_result)

