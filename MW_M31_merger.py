#!/usr/bin/env python
# coding: utf-8

import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '' and os.name != 'nt':
    print('No display found. Using non-interactive Agg backend for matplotlib')
    mpl.use('Agg')


import numpy as np
import matplotlib.pyplot as plt
import random
import math

import inspect
import subprocess
import sys

#finds script path
filename = inspect.getframeinfo(inspect.currentframe()).filename
script_path = os.path.dirname(os.path.abspath(filename))
#add progressbar path to possible paths to import modules
sys.path.insert(1, script_path + '/python-progressbar/')


plot_folder = script_path + '/plots'

if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)
    
gal_folder = script_path + '/galaxies/data'

if not os.path.exists(gal_folder):
    os.makedirs(gal_folder)


from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model
from amuse.couple import bridge


#ignore warnings (AmuseWarning at end of simulation)

import warnings

warnings.filterwarnings('ignore')


#defines parser for terminal usage
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--plot', 
                    help='Allow output plots', 
                    action='store_true')
parser.add_argument('--test', 
                    help='Use test initial conditions', 
                    action='store_true')
parser.add_argument('-f', 
                    help='Foo value, defined for jupyter notebook compatibility')
args = parser.parse_args()

PLOT = args.plot
TEST = args.test

if PLOT:
    print('Plots turned on', flush=True)


#progressbar import

import progressbar as pbar
import progressbar.widgets as pbwg


from galaxies import galaxies as gal


M_galaxy = 1.0e12 | units.MSun
R_galaxy = 10 | units.kpc
t_end = 300 | units.Myr

data_format = '.hdf5'

if not TEST:
    n_bulge = 10000
    n_disk = 10000
    n_halo = 20000
    
    mw_data = script_path + '/galaxies/data/MW_full' 
    m31_data = script_path + '/galaxies/data/M31_full'

    if not os.path.exists(mw_data) and not os.path.exists(m31_data):
        M31, MW, converter = gal.make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk, script_path, test=TEST)
        
    if os.path.exists(mw_data) and os.path.exists(m31_data):
        widgets = ['Loading galaxies data: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
        with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
            converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
            MW = read_set_from_file(mw_data, "hdf5")
            M31 = read_set_from_file(m31_data, 'hdf5')
        
if TEST:
    print('Using test initial conditions', flush=True)
    n_bulge = 750
    n_disk = 750
    n_halo = 1500
    
    mw_data = script_path + '/galaxies/data/MW_test' 
    m31_data = script_path + '/galaxies/data/M31_test'

    if not os.path.exists(mw_data) and not os.path.exists(m31_data):
        M31, MW, converter = gal.make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk, script_path, test=TEST)
        
    if os.path.exists(mw_data) and os.path.exists(m31_data):
        widgets = ['Loading galaxies data: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
        with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
            converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
            MW = read_set_from_file(mw_data, "hdf5")
            M31 = read_set_from_file(m31_data, 'hdf5')

test_disk = gal.test_particles(MW, n_halo, n_bulge, n_disk)
#gal.simulate_merger_with_particles(M31, MW, converter, n_halo, n_bulge, n_disk, t_end, script_path, plot=PLOT)
#stars = gal.nstars(50, 8.27806, 0.10, 230000)
gal.simulate_merger(M31, MW, converter, n_halo, t_end, script_path, plot=PLOT)
#mw_and_stars(galaxy2, stars, converter, n_halo, t_end)


if not TEST:
    from IGM import IGM_homogenous_Gadget2 as igm

    N1 = 5000
    N2 = 1000
    L = 10 | units.kpc
    rho = 1000 | units.MSun / (units.kpc)**3
    u = 1.6e+15 | (units.m)**2 / (units.s)**2

    widgets = ['Building IGM: ', pbwg.AnimatedMarker(), ' ',
               pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        sph_code = igm.setup_sph_code(Gadget2, N1, N2, L, rho, u)

    gal.merger_and_igm(galaxy1, galaxy2, converter, sph_code, n_halo, t_end, script_path, plot=PLOT)




