#!/usr/bin/env python

'''
Auxiliary functions go here
'''

import shutil, datetime
import inputs as ip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns


def track_simulation(iteration, nbr_iterations):
    """
    Creates a graph x=simtime y=total level of auxin etc. This can be useful to detect bugs (for example, if there is not synth nor degr of auxin, total value has to remain constant).
    
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

    if iteration == 0 :
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


def track_series(series_num, series_num_total):

    """
    Creates a line plot with representing the auxin profile in a chosen column (or row) of cells

    Params:
        series_num: 0-based number of the current simulation in the series
        series_num_total: total number of simulations in the series
    """
    
    values = ip.auxin[:,5].copy()
    ip.auxin_series_historic.append(values)

    if series_num == series_num_total - 1:
        fig2 = plt.figure()
        fig2, ax = plt.subplots(1, figsize=(5,5))
        ax.set_prop_cycle('color', plt.cm.Spectral(np.linspace(0,1,6)))
        for i in ip.auxin_series_historic:
            ax.plot(i)
        fig2.savefig('graphs/series.png', bbox_inches='tight')
        plt.close()


def write_to_log(timestamp):
    
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
    write_to_log()