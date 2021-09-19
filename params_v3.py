#!/usr/bin/env python

import numpy as np

# General
euler_h = .01              # Euler step size = h; (0 - 1]
#odeint_timepoints = np.linspace(0, euler_h, 2)
simulation_time = 24       # Arbitrary Units (AU) (let's assume it is hours of development)
img_dest_folder = 'images'
cell_plot_frequency = 10     # To do: change this to images / hour
create_video = True
create_gif = False
simulation_description = 'Testing UTG and WTF in combination to reproduce leaf margin observations'

# Heatmap ranges
auxin_range = (0, 95)       # This is only to map variable values to heatmap values
pin1_range = (0, 10)
cuc_range = (0, 9)
middle_domain = (0, 9)
adab_domain = (0, 9)

# Auxin diffusion
k_auxin_diffusion = 0.3 #0.3 0.12  # Relative amount of molecules that cross between two adjacent cells per cycle

# Auxin homeostasis
k_auxin_synth = 0 # Basal absolute amount of molecules synthesized per unit of time
k_auxin_degr = 0 #0.02 #0.01 #0.06 # Relative amount of molecules degraded per unit of time
k_cuc_auxin_synth = 1 #1 #0.5 #0.6 #0.3
k_md_auxin_synth = 0 #0.25 #0

# Auxin noise ( interval is [) and refers to iterations)
auxin_noise = {
  'limit': 0, #0.0025,
  'iteration_interval': (0,1)
}

# Auxin custom local synthesis (treated as absolute rate; interval is [) and refers to simulation time)
auxin_custom_synth = {
  "cells": ((5,5),(0,100)),
  "time_interval": (0,1000),
  "value": 0 #20
}
# Auxin custom local degradation (treated as relative rate; interval is [) and refers to simulation time)
auxin_custom_degr = {
  "cells": ((0,4),(100,100)),
  "time_interval": (0,1000),
  "value": 0 #0.05
}

# PIN1 localization
pin1_polarity = 'wtf_abley2016'   # 'multi' OR 'smith2006' OR 'wtf_abley2016'
k_UTG = 1.1 # 6 (6 in Bilsborough 2011, Smith 2006)
k_WTF_a = 10 #1500 #1 in Abley 2016 ('linear WTF')
k_WTF_b = .18 #0.2 #.005 in Abley 2016 ('linear WTF')
k_WTF_pin1_max = 9 #100000000000 #9
cuc_threshold_pin1 = 5

# PIN1 expression and activity
k_auxin_pin1 = 0.005 #0.003 #0.0001
k_cuc_pin1 = 0 #0.01
k_pin1_synth = 0
k_pin1_decay = 0.02 #0.1 # 0.004
k_pin1_transp = 0.02 #0.02 #0.01   # = auxin molecules transported / ( PIN1 molecule * cycle )

# CUC expression
k_cuc_synth = 0.01 #0.35
k_md_cuc = 0 #0.01
k_adab_cuc = 0
k_auxin_cuc = 0.005 #0.01
k_cuc_decay = 0

'''
Local synthesis or degradation (absolute or relative)
Here define a list of elements, each specifying the cell coordinates, synth/degr, abs/rel, cycles, etc

Example: marginSynthesis = new Local(type='synth', mode='abs', coords=(2,4,0), cycles=(0-4)) 
'''