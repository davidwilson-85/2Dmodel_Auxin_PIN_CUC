#!/usr/bin/env python

import func_auxin

# General
euler_h = 0.1               # Euler step size = h; (0 - 1]
nbr_iterations = 400
img_dest_folder = 'images/test'
cell_plot_frequency = 5

# Heatmap ranges
auxin_range = (0, 99)       # This is only to map variable values to heatmap values
pin1_range = (0, 10)
cuc_range = (0, 9)
middle_domain = (0, 9)

# Auxin diffusion
k_auxin_diffusion = 0.12 #0.3 0.12    		 # Relative amount of molecules that cross between two adjacent cells per cycle

# Auxin homeostasis
k_auxin_synth = 0			  	     	 # Basal absolute amount of molecules synthesized per cycle
k_cuc_auxin_synth = 0.1 #0.3
k_auxin_degr = 0 #0.06
# Auxin - other params
auxin_noise_factor = 0.025

# Auxin custom local synthesis
auxin_custom_synth = {
  "cells": ((6,0),(6,1)),
  "iterations": range(0,1000),
  "value": 50
}
# Auxin custom local degradation
auxin_custom_degr = {
  "cells": (),
  "iterations": range(0),
  "value": 0
}

# PIN1 localization
pin1_polarity = 'wtf_abley2016'   # 'multi' OR 'smith2006' OR 'wtf_abley2016'
k_UTG = 1.2             # 6 (Bilsborough 2011, Smith 2006), 1.3
k_WTF_a = 1 # 4E-3  1 in Abley 2016
k_WTF_b = 0.005
k_WTF_pin1_max = 9

# PIN1 expression and activity
k_auxin_pin1 = 0 #0.003 #0.0001
k_cuc_pin1 = 0 #0.01
k_pin1_decay = 0 #0.1 # 0.004
k_pin1_transp = 0.01 #0.01   # = auxin molecules transported / ( PIN1 molecule * cycle )

# CUC activity
cuc_on_pin1Pol = 0          # Maybe this is a dead parameter
cuc_threshold_pin1 = 5

# CUC expression
k_md_cuc = 0.11
k_auxin_cuc = 0.03
k_cuc_decay = 0.03

'''
Local synthesis or degradation (absolute or relative)
Here define a list of elements, each specifying the cell coordinates, synth/degr, abs/rel, cycles, etc

Example: marginSynthesis = new Local(type='synth', mode='abs', coords=(2,4,0), cycles=(0-4)) 
'''
