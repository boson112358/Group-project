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


folder = script_path + '/plots'

if not os.path.exists(folder):
    os.makedirs(folder)


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
parser.add_argument('-f', 
                    help='Foo value, defined for jupyter notebook compatibility')
args = parser.parse_args()

PLOT = args.plot

if PLOT:
    print('Plots turned on')


#progressbar import

import progressbar as pbar
import progressbar.widgets as pbwg


from galaxies import galaxies as gal


M_galaxy = 1.0e12 | units.MSun
R_galaxy = 10 | units.kpc
n_bulge = 10000
n_disk = 10000
n_halo = 20000
t_end = 200 | units.Myr
#t_end = 200

galaxy1, galaxy2, converter = gal.make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk)
stars = gal.nstars(50, 8.27806, 0.10, 230000)
#simulate_merger(galaxy1, galaxy2, converter, n_halo, t_end)
#mw_and_stars(galaxy2, stars, converter, n_halo, t_end)


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




