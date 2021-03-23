#!/usr/bin/env python

import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

'''
This is tha basis for testing whether my approach (a loop with sequential
update of each variable in the model) can give comparable results to using
an integrator like ODEint.
'''

def model(a,t):

    k_sy = 0.05
    k_di = 0.01
    k_de = 0.025

    a = a.reshape(2,2)

    #def func(x):
    #    return - x * k_de

    #da_dt = [func(x) for x in auxin]
    da_dt = [0,0,0,0]

    da_dt[0] = k_di * (a[0,1] + a[1,0] - 2 * a[0,0]) + k_sy
    da_dt[1] = k_di * (a[0,0] + a[1,1] - 2 * a[0,1])
    da_dt[2] = k_di * (a[0,0] + a[1,1] - 2 * a[1,0]) + k_de * a[1,0]
    da_dt[3] = k_di * (a[0,1] + a[1,0] - 2 * a[1,1])

    #a2 = np.array(da_dt).flatten()
    #return a2

    return da_dt

'''
a_a = np.ones(4)
mo = model(a_a,0)
print(mo)
for i in mo:
    print(i)

quit()
'''

# Create a heatmap using matplotlib's imshow()
def create_heatmap(data, filenum):
    fig = plt.imshow(data, cmap="plasma", vmin=0, vmax=3, interpolation='none')
    plt.savefig('test_imgs/img_' + str(filenum) + '.png', bbox_inches='tight')
    plt.close()

a_t0 = np.zeros(4)

t = np.linspace(0, 20, 20)
sim = odeint(model, a_t0, t)

print(sim)
quit()

plt.plot(t, sim[: , 0])
plt.plot(t, sim[: , 1])
plt.plot(t, sim[: , 2])
plt.plot(t, sim[: , 3])
plt.savefig('test_imgs/diff.png', bbox_inches='tight')

for i,item in enumerate(sim[:,:]):
    create_heatmap(item.reshape(2,2), i)

print('Done.')