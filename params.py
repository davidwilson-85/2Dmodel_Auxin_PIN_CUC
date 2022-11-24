#!/usr/bin/env python

import numpy as np

# General
euler_h = .1              # Euler step size = h; (0 - 1]
#odeint_timepoints = np.linspace(0, euler_h, 2)
simulation_time = 2       # Arbitrary Units (AU) (let's assume it is hours of development)
img_dest_folder = 'images'
cell_plot_frequency = 1     # To do: change this to images / hour
create_video = True
create_gif = False
simulation_description = ''

# Heatmap ranges
auxin_range = (0, 235)       # This is only to map variable values to heatmap values
pin1_range = (0, 16)
cuc_range = (0, 9)
middle_domain = (0, 9)
adab_domain = (0, 9)

# Auxin homeostasis
k_auxin_diffusion = .4 #.3 .12  # Rel. amount of molecules that cross between two adjacent cells per cycle
k_auxin_synth = 0 # Basal abs. amount of molecules synthesized per unit of time
k_auxin_degr = 0 #0.2 #0.02 #0.01 #0.06 # Rel. amount of molecules degraded per unit of time
k_cuc_auxin_synth = .2 #.75 #.75 #1 #.5 #.6 #.3
k_md_auxin_synth = 0 #.05 #.25 #0

# Auxin noise ( interval is [) and refers to iterations)
auxin_noise = {
  'limit': 0, #.0025,
  'iteration_interval': (0,1)
}

# Auxin custom local synthesis (treated as absolute rate; interval is [) and refers to simulation time)
auxin_custom_synth = {
  "cells": ((12,0), (12,1), (12,2), (12,3), (12,4), (12,5), (12,6), (12,7), (12,8), (12,9), (12,10)),
  "time_interval": (0,100000000),
  "value": 5 #2.5 #20
}
# Auxin custom local degradation (treated as relative rate; interval is [) and refers to simulation time)
auxin_custom_degr = {
  "cells": ((1,5),(0,5)),
  "time_interval": (0,10000),
  "value": .99 #.2 #.05
}

# PIN1 localization
pin1_polarity = 'wtf_abley2016'   # 'multi' OR 'smith2006' OR 'wtf_abley2016'
k_UTG = 2 #1.1 #1.3 # 6 (6 in Bilsborough 2011, Smith 2006)
k_WTF_a = 1 #10 #1500 #1 in Abley 2016 ('linear WTF')
k_WTF_b = .5 #.18 #.2 #.005 in Abley 2016 ('linear WTF')
k_WTF_pin1_max = 10 #100000000000 #9
cuc_threshold_pin1 = 5 # For dual PIN1 polarization

# PIN1 expression and activity
k_auxin_pin1 = 0 #.005 #.003 #.0001
k_cuc_pin1 = .02 #.1 #.01
k_pin1_synth = 0 #1
k_pin1_decay = 0 #.02 #.1 # .004
k_pin1_effi_basal = .02 #.02 #.01   # = auxin mol. transp. / ( PIN1 molecule * cycle )
k_pin1_effi_cuc = 0 #.04 #.1 # CUC effect on PIN1 efficiency [0 = no effect]

# CUC expression
k_cuc_synth = 0 #.01 #.35
k_md_cuc = 0 #.01
k_adab_cuc = 0
k_auxin_cuc = 0 #.0005 #.01
k_cuc_decay = 0

'''
Local synthesis or degradation (absolute or relative)
Here define a list of elements, each specifying the cell coordinates, synth/degr, abs/rel, cycles, etc

Example: marginSynthesis = new Local(type='synth', mode='abs', coords=(2,4,0), cycles=(0-4)) 
'''