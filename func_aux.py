#!/usr/bin/env python

'''
Auxiliary functions go here
'''

import shutil, datetime
import inputs as ip
import numpy as np


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