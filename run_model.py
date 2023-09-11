#!/usr/bin/env python

import time, os, shutil, datetime, sys
import numpy as np

import importlib

from main import main
#import params as pr
pr = importlib.import_module(sys.argv[1].split('.')[0], package=None)
import inputs as ip

import tests.check as check

# Do checks
check.check_dirs()

# Create unique identifier of simulation based on date and time
timestamp = str(datetime.datetime.now())[:19].replace(':','-').replace(' ','_')

# Create simulation dirs
sim_dir = 'output/sim_' + str(timestamp)
os.mkdir(sim_dir)
os.mkdir(sim_dir + '/images')

# Copy params file inside simualation dir
params_file = sys.argv[1]
shutil.copy(params_file, sim_dir)


def run_single(batch_sim_params=None):
	'''
	Runs a single simulation with the parameters specified in params.py
	'''

	# Force reload the module inputs to reset all values to those in the parameters file (this is critical).
	# Reassign parameters defined inside dict pr.batch_params and the parameter being overwritten in each simulation of the series.
	# The inputs (ip) module has to be reloaded after params have been overwritte
	if batch_sim_params is not None:
		importlib.reload(pr)
		for param, val in batch_sim_params.items():
			exec('pr.' + str(param) + '=' + str(val)) # exec() converts str to code
		importlib.reload(ip)
	
	# Run single simulation
	main(timestamp, sim_dir)


def run_series(batch_sim_params=None):
	'''
	Runs a series of simulations
	Parameters are defined in params.py
	Each simulation will share all parameters except for one
	The changing parameter and its value in each simulation of the series is defined in params.series_param_a.
	'''
	
	# Cleanup destination folder (remove and create)
	#shutil.rmtree(pr.img_dest_folder) 
	#os.mkdir(pr.img_dest_folder)

	# Retrieve parameter that varies and range of values that it takes
	param_a_space = np.linspace(pr.series_param_a['min'], pr.series_param_a['max'], pr.series_param_a['num_points'])
	series_num_total = len(param_a_space)

	# Run simulation series
	for key, series_val in enumerate(param_a_space):

		# Force reload the module inputs to reset all values to those in the parameters file (this is critical).
		# Reassign parameters defined inside dict pr.batch_params and the parameter being overwritten in each simulation of the series.
		# The inputs (ip) module has to be reloaded after params have been overwritte
		importlib.reload(pr)
		if batch_sim_params is not None:
			for param, val in batch_sim_params.items():
				exec('pr.' + str(param) + '=' + str(val)) # exec() converts str to code
		exec('pr.' + pr.series_param_a['name'] + '=' + str(series_val)) # exec() converts str to code
		importlib.reload(ip)
		
		# Run single simulation part of a series
		print('series_num: ' + str(key) + '; ' + pr.series_param_a['name'] + ' = ' + str(series_val))
		main(timestamp, sim_dir, series_num_total, series_num = key)


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

	for sim_id, sim in pr.batch_params.items():
		
		importlib.reload(pr)
		importlib.reload(ip)

		# Make subdir for simulation images
		sim_imgs_subdir = sim_dir + '/' + str(sim_id)
		os.mkdir(sim_imgs_subdir)

		print('\nbatch element: ' + str(sim_id) + '; params: ' + str(sim))

		# Run model and collect the output
		
		if pr.is_series == False:

			# Cleanup destination folder (remove and create)
			shutil.rmtree(sim_dir + '/images') 
			os.mkdir(sim_dir + '/images')
			
			run_single(sim)

			# Move and rename files
			shutil.copytree(sim_dir + '/images', sim_imgs_subdir, dirs_exist_ok=True)
			os.rename(sim_dir + '/auxin_profile.csv', sim_dir + '/auxin_profile_' + str(sim_id) + '.csv')
			os.rename(sim_dir + '/auxin_profile.png', sim_dir + '/auxin_profile_' + str(sim_id) + '.png')
			os.rename(sim_dir + '/levels.png', sim_dir + '/levels_' + str(sim_id) + '.png')
			os.rename(sim_dir + '/vid_' + timestamp + '.mp4', sim_dir + '/vid_' + str(sim_id) + '.mp4')

		if pr.is_series == True:
			
			run_series(sim)

			# Move and rename files
			shutil.copytree(sim_dir + '/images', sim_imgs_subdir, dirs_exist_ok=True)
			os.rename(sim_dir + '/auxin_profile_multiple.csv', sim_dir + '/auxin_profile_multiple_' + str(sim_id) + '.csv')
			os.rename(sim_dir + '/auxin_profile_multiple.png', sim_dir + '/auxin_profile_multiple_' + str(sim_id) + '.png')
	
	# Remove generic image folder
	shutil.rmtree(sim_dir + '/images')


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
* Check if the order in which the functions are called has an effect on output of simulations

'''
