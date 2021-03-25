#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

t = np.linspace(0, 500, 500)

a_t0 = 50
a_t0_grid = np.zeros(4)

def model(a,t):

    k = 0.25
    da_dt = a * k
    return da_dt


def model_grid(a,t):
    k_sy = 0.05
    k_di = 0.15
    k_de = 0.1
    a = a.reshape(2,2)
    #def func(x):
    #    return - x * k_de
    #da_dt = [func(x) for x in auxin]
    da_dt = [0,0,0,0]
    da_dt[0] = k_di * (a[0,1] + a[1,0] - 2 * a[0,0]) + k_sy
    da_dt[1] = k_di * (a[0,0] + a[1,1] - 2 * a[0,1])
    da_dt[2] = k_di * (a[0,0] + a[1,1] - 2 * a[1,0]) - k_de * a[1,0]
    da_dt[3] = k_di * (a[0,1] + a[1,0] - 2 * a[1,1])
    #a2 = np.array(da_dt).flatten()
    #return a2
    return da_dt

###########

#sol = odeint(model, a_t0, t)

sol_grid = odeint(model_grid, a_t0_grid, t)

plt.plot(t, sol_grid[: , 0], label='0,0')
plt.plot(t, sol_grid[: , 1], label='0,1')
plt.plot(t, sol_grid[: , 2], label='1,0')
plt.plot(t, sol_grid[: , 3], label='1,1')
plt.legend()
plt.savefig('tests/test.png', bbox_inches='tight')

print(sol_grid[: , :])