#!/usr/bin/env python

'''
Auxiliary functions go here
'''

import shutil, datetime
import inputs_v3 as ip
import numpy as np
import matplotlib.pyplot as plt


def track_simulation(iteration, nbr_iterations):
    """
    Creates a graph x=simtime y=total level of auxin. This can be useful to detect bugs (for example, if there is not synth nor degr of auxin, total value has to remain constant).
    
    Detects when simulation has reached stationary state. This can be known by comparing every sim step with the previous one. values to compare are the levels of auxin/PIN1/CUC. Comparisons are done cell wise, changes are considered as absolute, and the changes in all cells in the grid are added together. If combined absolute changes are less than a certain threshold value, stationary state has been reached. Simulation can be stoped then.
    """

    if iteration == 0 :
        ip.auxin_auxiliary = ip.auxin.copy()
    
    if iteration > 0:
        auxin_diff_abs = np.absolute(ip.auxin - ip.auxin_auxiliary)
        auxin_diff_abs_sum = auxin_diff_abs.sum()
        print(auxin_diff_abs_sum)
        ip.auxin_sum_per_step.append(auxin_diff_abs_sum)
        ip.auxin_auxiliary = ip.auxin.copy()
    
    if iteration == nbr_iterations:
        plt.plot(ip.auxin_sum_per_step)
        plt.savefig('test.png')


def write_to_log(timestamp):
    
    ''''
    Write contents of inputs and templates files to log file
    '''

    #timestamp = str(datetime.datetime.now())[:19].replace(':','-').replace(' ','_')
    name_log_file = 'sim_logs/' + timestamp + '_params'
    shutil.copy('params_v3.py', name_log_file)

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