
import numpy as np

# General
euler_h = .1             # Euler step size = h; (0 - 1]
simulation_time_total = 200     # Arbitrary Units (AU) (let's assume it is hours)
img_dest_folder = 'images'
create_cell_plots = 'multiple' # False / 'multiple' / 'only_last'
cell_plot_frequency = 1     # Simulation time units between plots
create_video = True
create_gif = False
simulation_description = 'UTG model used in the paper'

# Heatmap ranges
auxin_range = (0, 250)       # This is only to map variable values to heatmap values
pin1_range = (0, 21)
cuc_range = (0, 6)
middle_domain = (0, 9)
adab_domain = (0, 9)

# Auxin homeostasis
k_auxin_diffusion = .3 #.3 .12 # Rel. amount of molecules that cross between two adjacent cells per cycle
k_auxin_synth = 1 #1 #0.8 # Basal abs. amount of molecules synthesized per cell per unit of time
k_auxin_degr = .025 #.02 #0.2 #0.02 #0.01 #0.06 # Rel. amount of molecules degraded per unit of time
k_cuc_auxin_synth = 0 #0.42 #0.25 #0.5 #.75 #.75 #1 #.5 #.6 #.3
k_md_auxin_synth = 0 #.05 #.25 #0
AA_Vmax = 0      # Auxin on auxin creation: Maximum of Hill Function. Use 0 to turn off
AA_n = 0 #10        # Auxin on auxin creation: Hill coefficient, Steepness of the Hill function
AA_k_05 = 0 #50     # Auxin on auxin creation: [auxin] where 1/2 og Hill function Vmax is reached

# Auxin noise (limit is in %; iteration_interval is [) and refers to iterations)
auxin_noise = {
  'limit': 0,
  'iteration_interval': (0,1000000)
}
# Auxin perfect sources (interval is [) and refers to simulation time)
auxin_perfect_sources = {
    "active": True,
    "cells": ((-1,-1), (1,5)),
    "time_interval": (0,100000000),
    "value": 300
}
# Auxin perfect sinks (interval is [) and refers to simulation time)
auxin_perfect_sinks = {
    "active": False,
    "cells": ((10,5), (9,5)),
    "time_interval": (0,100000000),
    "value": 0
}
# Auxin custom local synthesis (treated as absolute rate; interval is [) and refers to simulation time)
auxin_custom_synth = {
    "cells": ((-1,-1), (1,5)),
    "time_interval": (0,100000000),
    "value": 0 #50 #2.5 #20
}
# Auxin custom local degradation (treated as relative rate; interval is [) and refers to simulation time)
auxin_custom_degr = {
    "cells": ((11,7), (1000,5)),
    "time_interval": (0,10000),
    "value": 0 #0.5 #.2 #.05
}

# PIN1 localization
pin1_polarity = 'utg_smith2006' # 'dual' OR 'smith2006' OR 'wtf_abley2016'
k_UTG = 1.2 #1.1 #1.3 # 6 (6 in Bilsborough 2011, Smith 2006)
k_WTF = 2
k_WTF_a = 1 #10 #1500 #1 in Abley 2016 ('linear WTF')
k_WTF_b = .005 #.18 #.2 #.005 in Abley 2016 ('linear WTF')
k_WTF_pin1_max = 15 #29 #9
cuc_threshold_pin1 = 5 # For dual PIN1 polarization # POSSIBLY DEPRECATED
md_on_pin1_UTG = True # template_middle_domain modulates k_UTG
k_netflux_on_fluxhistory = 1 # For function pin_wtf_instant()
k_fluxhistory_decay = .1 # For function pin_wtf_instant()
wtf_instant_exponential = True # False / True

# PIN1 expression and activity
k_pin1_synth = .1 #1
k_pin1_decay = .025 #.02 #.1 #.004
k_auxin_pin1 = .0001 #.003 #.0001
k_cuc_pin1 = .01 #.2 #.1 #.01
k_md_pin1 = 0 #.05
k_pin1_effi_basal = .007 #.02 #.01   # = auxin mol. transp. / ( PIN1 molecule * cycle )
k_pin1_effi_pho = .07   # Transport efficiency of PIN1-P molecules
pin1_pho_k05 = 5 #.04 #.1 # [CUC] where effect half PIN1 molecules are phosphorylated. Use very high value to turn off
pin1_pho_hc = 1 # Hill coefficient. Determines the steepness and shape of the Hill function

# CUC expression
k_cuc_synth = 0 #.01 #.35
k_md_cuc = .08 #0.1
k_pd_cuc = .1
k_adab_cuc = 0
k_auxin_cuc = .003 #.0005 #.01
k_cuc_decay = .1
cuc_noise = { # CUC noise (limit is in %; iteration_interval is [) and refers to iterations)
  'limit': 0,
  'iteration_interval': (0,1000000)
}

# Series of simulations for parameter value exploration
is_series = False # Specifies whether simulation is a single run or a series
series_param_a = { # If is_series = True, this overrides value of the chosen parameter
    'name': 'dummy',
	'min': 0,
	'max': 1,
	'num_points': 10
}

# Batch of simulations for automation
is_batch = True # Specifies whether simulation is a batch of simulations
batch_params = { # If is_series = True, this overrides value of the chosen parameter
    'wt': {'k_md_cuc': .08},
    'cuc': {'k_md_cuc': 0}
}

## Templates with initial values

#template_middle_domain = np.array([0,0,4,6,9,9,9,6,4,0,0], dtype=float)
template_middle_domain = np.array([0,0,4,6,8,9,8,6,4,0,0], dtype=float)

template_proximodistal_axis = np.array([0,0,0,0,0,0,1,3,5,8,9,9,9,9,9], dtype=float)

template_auxin = np.array(
    [
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0]
    ],
    dtype=float
)

template_cuc = np.array(
    [
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0]
    ],
    dtype=float
)

template_pin1 = np.array(
    [
        [
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4]
        ],
        [
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4]
        ],
        [
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4]
        ],
        [
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4],
            [4,4,4,4,4,4,4,4,4,4,4]
        ],
    ],
    dtype=float
)