#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

t = np.linspace(0, 30, 30)

a_t0 = 50


def model(a,t):

    k = 0.25
    da_dt = a * k
    return da_dt

sol = odeint(model, a_t0, t)

plt.plot(t, sol[: , 0])
plt.savefig('test.png', bbox_inches='tight')

print(sol[: , 0])