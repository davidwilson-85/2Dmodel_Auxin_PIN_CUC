#!/usr/bin/env python

import func_auxin

# General
euler_h = 0.1                 # Euler step size = h; (0 - 1]
nbr_iterations = 250
img_dest_folder = 'images/test'
cell_plot_frequency = 1

# Heatmap ranges
auxin_range = (0, 99)       # This is only to map variable values to heatmap values
pin1_range = (0, 9)
cuc_range = (0, 9)
middle_domain = (0, 9)

# Switches
PIN1_UTG = True
PIN1_WTF = 'linear' # linear, cuadratic...
AUX_LAX_transport = False
CUC = False

# Auxin diffusion
k_auxin_diffusion = 1 #0.12    		 # Relative amount of molecules that cross between two adjacent cells per cycle
# Auxin homeostasis
k_auxin_synth = 0			  	     	 # Basal absolute amount of molecules synthesized per cycle
k_cuc_yuc = 0 #0.1
k_auxin_degr = 0.05
# Auxin - other params
auxin_noise_factor = 0
# Auxin custom local synthesis
auxin_custom_synth = {
  "cells": ((3,0),(3,1)),
  "iterations": (),
  "value": 10
}
# Auxin custom local degradation
auxin_custom_degr = {
  "cells": (),
  "iterations": (),
  "value": 0
}

# PIN1 localization/activity
pin1_polarity = 'multi'   # 'multi' OR 'smith2006' OR 'ratio' OR 'wtf'
k_UTG = 1.5             # 6 (Bilsborough 2011, Smith 2006), 1.3
k_WTF = 200000
k_pin1_transp = 0.035        # = Nbr auxin molecules transported / ( PIN1 molecule * cycle ); used values=0.01
# PIN1 expression
k_auxin_pin1 = 0.003 #0.0001
k_cuc_pin1 = 0 #0.0002
k_pin1_decay = 0.1 # 0.004

# CUC activity
cuc_on_pin1Pol = 0          # Maybe this is a dead parameter
cuc_threshold_pin1 = 1000
# CUC expression
k_md_cuc = 0.1
k_auxin_cuc = 0.03
k_cuc_decay = 0.03

# Local synthesis or degradation (absolute or relative)
	# Here define a list of elements, each specifying the cell coordinates, synth/degr, abs/rel, cycles, etc

	# Example: marginSynthesis = new Local(type='synth', mode='abs', coords=(2,4,0), cycles=(0-4)) 

