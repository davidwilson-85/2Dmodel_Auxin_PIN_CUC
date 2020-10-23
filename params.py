#!/usr/bin/env python

# General
euler_h = 0.1                 # Euler step size = h; (0 - 1]
nbr_iterations = 1000
img_dest_folder = 'images/test'
cell_plot_frequency = 10

# Heatmap ranges
auxin_range = (0, 99)       # This is only to map variable values to heatmap values
pin1_range = (0, 9)
cuc_range = (0, 9)

# Switches
PIN1_UTG = True
PIN1_WTF = 'linear' # linear, cuadratic...
AUX_LAX_transport = False
CUC = False

# Auxin diffusion
k_auxin_diffusion = 0.2 #0.12    		 # Relative amount of molecules that cross between two adjacent cells per cycle
# Auxin homeostasis
k_auxin_synth = 0			  	     	 # Basal absolute amount of molecules synthesized per cycle
k_cuc_yuc = 0 #0.01           WHAT IS THIS????????????
k_auxin_decay = 0.1
# Auxin - other params
auxin_noise_factor = 0

# PIN1 localization/activity
pin1_pol_mode = 'wtf'   # 'smith2006' OR 'ratio' OR 'wtf'
k_UTG = 1.3             # 6 (Bilsborough 2011, Smith 2006), 1.3
k_WTF = 2000
k_pin1_transp = 0.035        # = Nbr auxin molecules transported / ( PIN1 molecule * cycle ); used values=0.01
# PIN1 expression
k_auxin_pin1 = 0 #0.00006 #0.0001
k_cuc_pin1 = 0 #0.0005 #0.0002
k_pin1_decay = 0 #0.003 # 0.004

# CUC activity
cuc_on_pin1Pol = 0
cuc_threshold_pin1 = 500
# CUC expression
k_md_cuc = 0 #0.01
k_auxin_cuc = 0 #0.0019
k_cuc_decay = 0 #0.01