#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model_diff(a,t):

    k_s = 0.05
    k_d = 0.02
    k_t = 0.005

    a1 = a[0]
    a2 = a[1]

    da1_dt = k_s + k_d * a2 - k_d * a1
    da2_dt = k_d * a1 - k_d * a2 - k_t * a2

    return [da1_dt, da2_dt]

a_t0 = [10,0]

t_sim1 = np.linspace(0, 1000, 10000)
sim1 = odeint(model_diff, a_t0, t_sim1)

plt.plot(t_sim1, sim1[: , 0])
plt.plot(t_sim1, sim1[: , 1])
plt.savefig('diff.png', bbox_inches='tight')