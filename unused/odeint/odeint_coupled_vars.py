#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def hiv(x,t):

    kr1 = 1e5
    kr2 = 0.1
    kr3 = 2e-7
    kr4 = 0.5
    kr5 = 5
    kr6 = 100

    h = x[0]
    i = x[1]
    v = x[2]

    p = kr3 * h * v

    dhdt = kr1 - kr2 * h - p
    didt = p - kr4 * i 
    dvdt = - p - kr5 * v + kr6 * i

    return [dhdt, didt, dvdt]


def model_auxin_AUX1(x, t):

    k1 = 1
    k2 = 0.1
    k3 = 0.05
    k4 = 0.1
    k5 = 0.2

    auxin = x[0] 
    aux1 = x[1]
    
    dA_dt = k1 - k2 * auxin - auxin * aux1 * k3
    dAUX1_dt = auxin * k4 - aux1 * k5
    
    return [dA_dt, dAUX1_dt]

x_t0 = [10,2]
t = np.linspace(0, 100, 1000)
sim = odeint(model_auxin_AUX1, x_t0, t)

plt.plot(t, sim[: , 0])
plt.plot(t, sim[: , 1])
plt.savefig('auxin.png', bbox_inches='tight')

quit()

#print(hiv([1e6, 0, 100],0))

x0 = [1e6, 0, 100]
t = np.linspace(0, 15, 1000)
x = odeint(hiv, x0, t)

h = x[:, 0]
i = x[:, 1]
v = x[:, 2]

plt.semilogy(t, h)
plt.semilogy(t, i)
plt.semilogy(t, v)

plt.savefig('test.png', dpi=400, bbox_inches='tight', facecolor='w')