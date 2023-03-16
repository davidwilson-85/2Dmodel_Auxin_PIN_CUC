#!/usr/bin/env python

import os, shutil
import numpy as np
import params as pr
import TEST_main
#from sim_logs import params_2022_12_20_22_53_59 as pr # To re-run logged simulation


# Cleanup destination folder (remove and create)
#shutil.rmtree(pr.img_dest_folder) 
#os.mkdir(pr.img_dest_folder)

if pr.is_series == False:

	# Run single simulation
	TEST_main.run()

if pr.is_series == True:

	# Retrieve parameter that varies and range of values that it takes
	#param_a_space = np.linspace(pr.series_param_a['min'], pr.series_param_a['max'], pr.series_param_a['num_points'])
	param_a_space = [0.777, 1, 0.111, 0.777, 1]
	# Calculate total number of simulations in the series and set parameter
	series_num_total = len(param_a_space)

	print(param_a_space)
	print(series_num_total)
	for i in param_a_space: print(i)

	# Run simulation series
	for key, val in enumerate(param_a_space):
		exec('pr.' + pr.series_param_a['name'] + '=' + str(val)) # exec() converts str to code
		print('series_num: ' + str(key) + '; ' + pr.series_param_a['name'] + ' = ' + str(val))
		#main.run(key, series_num_total)
		os.system(f"python3 TEST_main.py {key} {val}")
