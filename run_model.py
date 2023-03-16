#!/usr/bin/env python

import time, os, shutil, datetime
import numpy as np

import importlib

from main import run
import params as pr
import inputs as ip

import tests.check as check

# Do checks
check.check_dirs()


def run_single(batch_sim_params = None):
	'''
	Runs a single simulation with the parameters specified in params.py
	'''
	
	# Cleanup destination folder (remove and create)
	shutil.rmtree(pr.img_dest_folder) 
	os.mkdir(pr.img_dest_folder)

	# Force reload the module inputs to reset all values to those in the parameters file (this is critical).
	# Reassign parameters defined inside dict pr.batch_params and the parameter being overwritten in each simulation of the series.
	if batch_sim_params is not None:
		importlib.reload(pr)
		importlib.reload(ip)
		for param, val in batch_sim_params.items():
			exec('pr.' + str(param) + '=' + str(val)) # exec() converts str to code
	
	# Run single simulation
	run()


def run_series(batch_sim_params = None):
	'''
	Runs a series of simulations
	Parameters are defined in params.py
	Each simulation will share all parameters except for one
	The changing parameter and its value in each simulation of the series is defined in params.series_param_a.
	'''
	
	# Cleanup destination folder (remove and create)
	shutil.rmtree(pr.img_dest_folder) 
	os.mkdir(pr.img_dest_folder)

	# Retrieve parameter that varies and range of values that it takes
	param_a_space = np.linspace(pr.series_param_a['min'], pr.series_param_a['max'], pr.series_param_a['num_points'])
	series_num_total = len(param_a_space)

	# Run simulation series
	for key, series_val in enumerate(param_a_space):
		#val = float(val)

		# Force reload the module inputs to reset all values to those in the parameters file (this is critical).
		# Reassign parameters defined inside dict pr.batch_params and the parameter being overwritten in each simulation of the series.
		importlib.reload(pr)
		importlib.reload(ip)
		if batch_sim_params is not None:
			for param, val in batch_sim_params.items():
				exec('pr.' + str(param) + '=' + str(val)) # exec() converts str to code
		exec('pr.' + pr.series_param_a['name'] + '=' + str(series_val)) # exec() converts str to code
		
		# Run
		print('series_num: ' + str(key) + '; ' + pr.series_param_a['name'] + ' = ' + str(series_val))
		run(series_num_total, series_num = key)


def run_batch():
	'''
	Runs a batch of simulations
	A batch consists of a set of single or simulations or series of simulations
	The parameters for each simulations in the batch are defined in params.batch_params
	This function calls run_single() or run_series() and collects the output in folder out_batch for review at the end of the batch

	Steps:
	1. Define parameters to inject in the model, overriding those specified indide dict pr.batch_params for single simulations.
	2. Run model.
	3. Copy and rename desired output files to keep them for later.
	'''

	# Define batch unique id and create folder to contain all simulations
	batch_id = str(datetime.datetime.now())[:19].replace(':','-').replace(' ','_')
	batch_dir = 'out_batch/batch_' + str(batch_id)
	os.mkdir(batch_dir)

	for sim_id, sim in pr.batch_params.items():
		
		importlib.reload(pr)
		importlib.reload(ip)

		# Make subdir for simulation
		sim_dir = 'out_batch/batch_' + str(batch_id) + '/' + str(sim_id)
		os.mkdir(sim_dir)

		'''# Update parameters defined inside dict pr.batch_params
		for param, val in sim.items():
			exec('pr.' + str(param) + '=' + str(val)) # exec() converts str to code'''

		print('\nbatch_num: ' + str(sim_id) + '; params: ' + str(sim))

		# Run model and collect the output
		
		if pr.is_series == False:
			
			run_single(sim)

			# Copy output to batch/simulation folder
			shutil.copytree(pr.img_dest_folder, sim_dir, dirs_exist_ok=True)
			shutil.copy('graphs/auxin_profile.png', batch_dir + '/auxin_profile_' + str(sim_id) + '.png')
			shutil.copy('graphs/auxin_profile.csv', batch_dir + '/auxin_profile_' + str(sim_id) + '.csv')
			shutil.copy('graphs/levels.png', batch_dir + '/levels_' + str(sim_id) + '.png')

		if pr.is_series == True:
			
			run_series(sim)

			# Copy output to batch/simulation folder
			shutil.copytree(pr.img_dest_folder, sim_dir, dirs_exist_ok=True)
			shutil.copy('graphs/auxin_profile_multiple.png', batch_dir + '/auxin_profile_multiple_' + str(sim_id) + '.png')
			shutil.copy('graphs/auxin_profile_multiple.csv', batch_dir + '/auxin_profile_multiple_' + str(sim_id) + '.csv')



if __name__ == '__main__':

	global_start_time = time.time()

	if pr.is_batch == False:
		
		if pr.is_series == False:
			run_single()

		if pr.is_series == True:
			run_series()

	if pr.is_batch == True:

		run_batch()
	
	print("\nTotal time: %s minutes" % ( (time.time() - global_start_time) / 60) )
		


'''
TO DO:
* Try staggered cells
* Try PD growth?
* Try to move more things (like custom auxin synth and degradation) to ODEint model
* Check if the order in which the functions are called has an effect on output of simulations

'''
