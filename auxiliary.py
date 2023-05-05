#!/usr/bin/env python

'''
Auxiliary functions go here
'''

import shutil, importlib, sys

import inputs as ip
#import params as pr
pr = importlib.import_module(sys.argv[1].split('.')[0], package=None)

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import seaborn as sns
from PIL import Image, ImageDraw


def track_simulation(iteration, nbr_iterations):
    """
    Creates a graph x=simtime y=total level of auxin, etc. This can be useful to detect bugs (for example, if there was not synth nor degr of auxin in a simulation, total value should remain constant).
    
    Detects when simulation has reached stationary state. This can be known by comparing every sim step with the previous one. Values to compare are the levels of auxin/PIN1/CUC. Comparisons are done cell wise, changes are considered as absolute, and the changes in all cells in the grid are added together. If combined absolute changes are less than a certain threshold value, stationary state has been reached. Simulation can be stopped then.
    """

    # Total auxin in system
    auxin_allcells = ip.auxin.sum()
    ip.auxin_allcells_historic.append(auxin_allcells)
    # Total PIN1 in system
    pin1_allcells = ip.pin1.sum()
    ip.pin1_allcells_historic.append(pin1_allcells)
    # Total CUC in system
    cuc_allcells = ip.cuc.sum()
    ip.cuc_allcells_historic.append(cuc_allcells)

    if iteration == 0:
        ip.auxin_auxiliary = ip.auxin.copy()
        ip.pin1_auxiliary = ip.pin1.copy()
        ip.cuc_auxiliary = ip.cuc.copy()
    
    if iteration > 0:
        auxin_increment = np.absolute(ip.auxin - ip.auxin_auxiliary)
        auxin_increment_allcells = auxin_increment.sum()
        ip.auxin_incr_allcells_historic.append(auxin_increment_allcells)
        ip.auxin_auxiliary = ip.auxin.copy()

        pin1_increment = np.absolute(ip.pin1 - ip.pin1_auxiliary)
        pin1_increment_allcells = pin1_increment.sum()
        ip.pin1_incr_allcells_historic.append(pin1_increment_allcells)
        ip.pin1_auxiliary = ip.pin1.copy()

        cuc_increment = np.absolute(ip.cuc - ip.cuc_auxiliary)
        cuc_increment_allcells = cuc_increment.sum()
        ip.cuc_incr_allcells_historic.append(cuc_increment_allcells)
        ip.cuc_auxiliary = ip.cuc.copy()
    
    if iteration == nbr_iterations:
        fig1 = plt.figure()
        fig1, axes = plt.subplots(6, figsize=(5,10), sharex=True)
        axes[0].plot(ip.auxin_allcells_historic)
        axes[1].plot(ip.auxin_incr_allcells_historic)
        axes[2].plot(ip.pin1_allcells_historic)
        axes[3].plot(ip.pin1_incr_allcells_historic)
        axes[4].plot(ip.cuc_allcells_historic)
        axes[5].plot(ip.cuc_incr_allcells_historic)

        axes[0].set(ylabel='Total auxin')
        axes[1].set(ylabel='Auxin abs. change')
        axes[2].set(ylabel='Total PIN1')
        axes[3].set(ylabel='PIN1 abs. change')
        axes[4].set(ylabel='Total CUC')
        axes[5].set(ylabel='CUC abs. change')

        for ax in axes:
            ax.ticklabel_format(useOffset=False, style='plain')
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            #ax.set_ylim(bottom=0)

        plt.xlabel('Simulation iteration')
        fig1.savefig('graphs/levels.png', bbox_inches='tight')
        plt.close()

        print('Created plot: levels')


def make_kymograph(iteration, nbr_iterations):
    '''docstring'''

    if iteration == 0:
        global data_kymograph
        data_kymograph = []

    # Create the y values of graph
    data_kymograph.append(ip.auxin[:,5].copy())

    if iteration == nbr_iterations:

        # make fig
        width = 5000
        height = 700
        patch_w = 1
        patch_h = 40
        im = Image.new('RGB', size=(width,height))
        draw = ImageDraw.Draw(im, 'RGBA')

        for t, profile in enumerate(data_kymograph):
            for cell, auxin in enumerate(profile):
                c = int(auxin)
                x_ref = t * patch_w
                y_ref = cell * patch_h
                draw.polygon([
                    (x_ref,y_ref),
                    (x_ref+patch_w,y_ref),
                    (x_ref+patch_w,y_ref+patch_h),
                    (x_ref,y_ref+patch_h)],
                    fill=(c, c, c)
                )
        
        # Save image
        im.save('graphs/kymograph.png')


def create_line_plot_single(timestamp, series_num, series_num_total):

    """
    Creates a line plot representing the auxin profile in a chosen column (or row) of cells

    Params:
        * timestamp: unique identifier of the simulation based on the date and time at which the simulation is run, and to be printed on the image to be able to recover the parameters used
        * series_num: 0-based number of the current simulation in the series
        * series_num_total: total number of simulations in the series
        * data arrays of the model
    """

    if series_num == 0:
        global auxin_series_historic
        auxin_series_historic = []

    # Create y axis points of graph
    tissue_rows = ip.auxin.shape[0]
    x = np.linspace(1, tissue_rows, tissue_rows)
    
    values = ip.auxin[:,5].copy()

    # Convert to Pandas dataframe and save as a csv file
    df = pd.DataFrame(values)
    df.to_csv("graphs/auxin_profile.csv", index=False)

    # Make plot
    if series_num == series_num_total - 1:
        fig2 = plt.figure()
        fig2, ax = plt.subplots(1, figsize=(3.5,5))
        # Choices that worked: 'Spectral', 'coolwarm', 'plasma', 'viridis', 'jet', 'brg'
        ax.set_prop_cycle('color', plt.cm.viridis(np.linspace(0,1,series_num_total)))
        ax.plot(values,x, c='gray')
        ax.invert_yaxis()
        plt.xlabel('[Auxin] A.U.')
        plt.ylabel('Cell row')
        ax.annotate(timestamp, xy=(2, 1), xytext=(0.01, .99), textcoords='figure fraction', va='top', ha='left')
        fig2.savefig('graphs/auxin_profile.png', bbox_inches='tight')
        #plt.colorbar(label="param value", orientation="vertical")
        plt.close()

        print('Created plot: auxin profile')


def create_line_plot_multi(timestamp, series_num, series_num_total):

    """
    Function is called in the last iteration of each simulation and gathers the auxin values in the desired cells
    In the last series of the simulation, it creates a line plot with representing the auxin profile in a chosen column (or row) of cells

    Params:
        series_num: 0-based number of the current simulation in the series
        series_num_total: total number of simulations in the series
    """

    if series_num == 0:
        global auxin_series_historic
        auxin_series_historic = []

    # Create the x values of graph
    values = ip.auxin[:,5].copy()
    auxin_series_historic.append(values)

    if series_num == series_num_total - 1:

        # Convert to Pandas dataframe, transpose, and save as a csv file
        df = pd.DataFrame(auxin_series_historic).T
        df.to_csv("graphs/auxin_profile_multiple.csv", index=False)

        # Create the y values of graph
        tissue_rows = ip.auxin.shape[0]
        x = np.linspace(1, tissue_rows, tissue_rows)
        
        # Define colormap
        cmap = plt.get_cmap('viridis', series_num_total)
        
        fig2, ax = plt.subplots(1, 1, figsize=(5, 5))
        
        for n, y in enumerate(auxin_series_historic):
            ax.plot(y, x, c=cmap(n))
        ax.invert_yaxis()
        ax.annotate(timestamp, xy=(2, 1), xytext=(0.01, .99), textcoords='figure fraction', va='top', ha='left')
        plt.xlabel('[Auxin] A.U.')
        plt.ylabel('Cell row')
        
        # Colorbar
        # Normalizer
        norm = mpl.colors.Normalize(vmin=pr.series_param_a['min'], vmax=pr.series_param_a['max'])
        # creating ScalarMappable
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ticks=np.linspace(pr.series_param_a['min'], pr.series_param_a['max'], pr.series_param_a['num_points']), label=pr.series_param_a['name'])

        fig2.savefig('graphs/auxin_profile_multiple.png', bbox_inches='tight')

        print('Created plot: auxin profile')


def write_to_log(timestamp):
    
    ''''
    Write contents of inputs and templates files to log file
    '''

    timestamp = str(timestamp).replace('-','_')
    name_log_file = 'sim_logs/params_' + timestamp + '.py'
    params_file = sys.argv[1]
    shutil.copy(params_file, name_log_file)


def write_to_log_old(timestamp):
    
    ''''
    Write contents of inputs and templates files to log file
    '''

    #timestamp = str(datetime.datetime.now())[:19].replace(':','-').replace(' ','_')
    name_log_file = 'sim_logs/' + timestamp + '_params'
    shutil.copy('params.py', name_log_file)

    with open('templates/2D/template_auxin', mode='r') as template_auxin:
        template_auxin_contents = template_auxin.read()

    with open('templates/2D/template_pin1', mode='r') as template_pin1:
        template_pin1_contents = template_pin1.read()

    with open('templates/2D/template_cuc', mode='r') as template_cuc:
        template_cuc_contents = template_cuc.read()

    with open('templates/2D/template_middle_domain', mode='r') as template_middle_domain:
        template_middle_domain_contents = template_middle_domain.read()
    
    with open('templates/2D/template_adab_domain', mode='r') as template_adab_domain:
        template_adab_domain_contents = template_adab_domain.read()

    with open(name_log_file, 'a') as log_file:
        
        log_file.write('\n\n\n***** template_auxin *****\n')
        log_file.write(template_auxin_contents)
        log_file.write('\n\n')

        log_file.write('\n\n***** template_pin1 *****\n')
        log_file.write(template_pin1_contents)
        log_file.write('\n\n')

        log_file.write('\n\n***** template_cuc *****\n')
        log_file.write(template_cuc_contents)
        log_file.write('\n\n')

        log_file.write('\n\n***** template_middle_domain *****\n')
        log_file.write(template_middle_domain_contents)
        log_file.write('\n\n')

        log_file.write('\n\n***** template_adaxial/abaxial_domain *****\n')
        log_file.write(template_adab_domain_contents)
        log_file.write('\n\n')


def save_ndarray():

    '''
    The standard way with Numbpy gives me an error:
    'Cannot load file containing pickled data when allow_pickle=False'
    And makes things complicated in other aspects too
    '''

    #with open('test.npy', 'wb') as f:
    #    np.save(f, np.array([1, 2]))
    
    #with open('test.npy', 'rb') as f:
    # a = np.load(f)

    with open('templates/2D/template_auxin_1', 'a') as file:

        for y in range(ip.tissue_rows):

            row_values = []
            
            for x in range(ip.tissue_columns):
                row_values.append(str(ip.auxin[y,x]/10)[:4] + ',')
            
            file.write(''.join(row_values)[:-1])
            
            if y < ip.tissue_rows - 1:
                file.write('\n')
    

if __name__ == '__main__':
    #write_to_log('test_timestamp')
    create_line_plot_multi(0, 1)