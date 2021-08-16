#!/usr/bin/env python

import os

def check_dirs():

    """
    Check and create folders for images, videos, sim_logs and templates
    
    This function is executed every time the program runs but will only tipically do stuff after cloning from Github or after accidentally deleting any of the folders listed above.
    """

    if not os.path.exists('images'):
        os.makedirs('images')

    if not os.path.exists('videos'):
        os.makedirs('videos')
    
    if not os.path.exists('sim_logs'):
        os.makedirs('sim_logs')
    
    if not os.path.exists('templates/1D'):
        os.makedirs('templates/1D')
    
    if not os.path.exists('templates/2D'):
        os.makedirs('templates/2D')