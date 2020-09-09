#!/usr/bin/env python


# General
nbr_iterations = 500
img_dest_folder = 'images/test'
cell_plot_frequency = 5

# Switches
PIN1_UTG = True
PIN1_WTF = 'linear' # linear, cuadratic...
AUX_LAX_transport = False
CUC = False

# Auxin
auxin_range = (0, 99)       # This is only to map variable values to heatmap values
k_auxin_diffusion = 0.12     # Relative amount of molecules that cross between two adjacent cells per cycle
auxin_synthesis = 0.1       # Absolute amount of molecules synthesized per cycle
auxin_destruction = 0.1     # Absolute amount of molecules destroyed per cycle
auxin_noise_factor = 0

# PIN1
utg_function = 'ratio'   # 'smith2006' OR 'ratio'
pin1_range = (0, 9)
k_UTG = 1.3 #6			   	 # Relative amount of molecules that can change cell face per cycle
k_pin1_transp = 0.002        # = Nbr auxin molecules transported / ( PIN1 molecule * cycle ); used values=0.01

# CUC
cuc_range = (0, 9)          # This is only to map variable values to heatmap values
auxin_on_cuc = 0
cuc_on_pin1Pol = 0
cuc_threshold_pin1 = 500


