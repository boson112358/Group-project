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
import sys

#finds script path
filename = inspect.getframeinfo(inspect.currentframe()).filename
SCRIPT_PATH = os.path.dirname(os.path.abspath(filename))
#add progressbar path to possible paths to import modules
sys.path.insert(1, SCRIPT_PATH + '/python-progressbar/')


plot_folders = [SCRIPT_PATH + '/plots', 
                SCRIPT_PATH + '/plots/merger_plots',
                SCRIPT_PATH + '/plots/solar_system_plots',
                SCRIPT_PATH + '/plots/merger_contour/',
                SCRIPT_PATH + '/plots/zoomed_merger_plots/']

for folder in plot_folders:
    if not os.path.exists(folder):
        os.makedirs(folder)
    
gal_folder = SCRIPT_PATH + '/galaxies/data'

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
parser.add_argument('--solar', 
                    help='Add solar system', 
                    action='store_true')
parser.add_argument('--disk', 
                    help='Add test disk to MW', 
                    action='store_true')
parser.add_argument('--igm', 
                    help='Add IGM', 
                    action='store_true')
parser.add_argument('--test', 
                    help='Use test initial conditions', 
                    action='store_true')
parser.add_argument('-f', 
                    help='Foo value, defined for jupyter notebook compatibility')
args = parser.parse_args()

PLOT = args.plot
SOLAR = args.solar
DISK = args.disk
IGM = args.igm
TEST = args.test

if PLOT:
    print('Plots turned on', flush=True)


#progressbar import

import progressbar as pbar
import progressbar.widgets as pbwg


#Galaxies starting conditions
if not TEST:
    n_bulge = 20000
    n_disk = 20000
    n_halo = 40000
elif TEST:
    print('Using test initial conditions', flush=True)
    n_bulge = 10000
    n_disk = 10000
    n_halo = 20000

mw_parameters = {'name': 'mw',
                 'n_halo': n_halo,
                 'halo_outer_radius' : 244.48999 | units.kpc,
                 'disk_number_of_particles' : n_disk,
                 'disk_mass' : 19.66 * 2.33 * 10e9 | units.MSun,
                 'disk_scale_length' : 2.806 | units.kpc,
                 'disk_outer_radius' : 30 | units.kpc, 
                 'disk_scale_height_sech2' : 0.409 | units.kpc,
                 'bulge_scale_radius' : 0.788 | units.kpc,
                 'bulge_number_of_particles' : n_bulge}

m31_parameters = {'name': 'm31_not_displaced',
                  'n_halo': n_halo,
                  'halo_outer_radius' : 201.619995 | units.kpc,
                  'disk_number_of_particles' : n_disk,
                  'disk_mass' : 33.40 * 2.33 * 10e9 | units.MSun,
                  'disk_scale_length' : 5.577 | units.kpc,
                  'disk_outer_radius' : 30 | units.kpc, 
                  'disk_scale_height_sech2' : 0.3 | units.kpc,
                  'bulge_scale_radius' : 1.826 | units.kpc,
                  'bulge_number_of_particles' : n_bulge}

#simulation parameters
scale_mass_galaxy = 1.0e12 | units.MSun
scale_radius_galaxy = 10 | units.kpc
t_end = 8000 | units.Myr

#Solar system starting conditions
n_stars = 10                                       #How many particles we will add
solar_radial_distance = 8.2
solar_position = (-np.sqrt(0.5 * solar_radial_distance**2),
                  np.sqrt(0.5 * solar_radial_distance**2),
                  0) | units.kpc                 #If we displace MW, we can do it through this vector aswell
system_radius = 0.100 | units.kpc                #neighbourhood in which we distribute our stars
solar_tang_velocity = 220 | (units.km/units.s)   #This is roughly the velocity of solarsystem around MW

#Intergal medium starting conditions
N1 = 5000
N2 = 1000
L = 10 | units.kpc
rho = 1000 | units.MSun / (units.kpc)**3
u = 1.6e+15 | (units.m)**2 / (units.s)**2


#M31 displacement
rotation = np.array([[0.7703,  0.3244,  0.5490],
                     [-0.6321, 0.5017,  0.5905],
                     [-0.0839, -0.8019, 0.5915]])

traslation = [-379.2, 612.7, 283.1] | units.kpc
radial_velocity = 117 * np.array([0.4898, -0.7914, 0.3657]) | units.kms
transverse_velocity = 50 * np.array([0.5236, 0.6024, 0.6024]) | units.kms


from galaxies import galaxies as gal
from stars import solar_system as sol
from intergal_medium import IGM_homogenous_Gadget2 as igm


#Galaxy initialization
converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)

if not TEST:
    mw_data_path = SCRIPT_PATH + '/galaxies/data/{}_full'.format(mw_parameters['name'])
    m31_not_displaced_data_path = SCRIPT_PATH + '/galaxies/data/{}_full'.format(m31_parameters['name'])
    
elif TEST:
    mw_data_path = SCRIPT_PATH + '/galaxies/data/{}_test'.format(mw_parameters['name'])
    m31_not_displaced_data_path = SCRIPT_PATH + '/galaxies/data/{}_test'.format(m31_parameters['name']) 

if os.path.exists(mw_data_path) and os.path.exists(m31_not_displaced_data_path):
    widgets = ['Loading galaxies data: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        mw = read_set_from_file(mw_data_path, "hdf5")
        m31_not_displaced = read_set_from_file(m31_not_displaced_data_path, 'hdf5')
else:
    mw = gal.make_galaxy(converter, mw_parameters, SCRIPT_PATH, test=TEST)
    m31_not_displaced = gal.make_galaxy(converter, m31_parameters, SCRIPT_PATH, test=TEST)


if all(value == False for value in [SOLAR, DISK, IGM]):
    m31 = gal.displace_galaxy(m31_not_displaced, rotation, traslation, radial_velocity, transverse_velocity)
    print('Simulating merger with no additional components (t = {} Myr) ...'.format(int(np.round(t_end.value_in(units.Myr), 
                                                                                                 decimals=0))), flush=True)
    gal.simulate_merger(mw, m31, n_halo, t_end, SCRIPT_PATH, plot=PLOT)
    
if DISK:
    print('Simulating merger with disk test particles ...', flush=True)
    test_disk = gal.test_particles(MW, n_halo, n_bulge, n_disk)
    gal.simulate_merger_with_particles(M31, MW, converter, n_halo, n_bulge, n_disk, t_end, SCRIPT_PATH, plot=PLOT)

if SOLAR:
    print('Simulating MW with solar system ...', flush=True)
    MW.position += [100.0, 100.0, 0] | units.kpc
    t_end = 400 | units.Myr
    #MW.rotate(-np.pi/4, -np.pi/4, 0.0)
    stars = sol.make_solar_system(n_stars, solar_position, system_radius, MW_velocity_vector, solar_tang_velocity)
    gal.mw_and_stars(MW, stars, converter, sol.leapfrog_alg, n_halo, t_end, SCRIPT_PATH, plot=PLOT)


if IGM:
    print('Simulating merger with IGM ...', flush=True)
    widgets = ['Building IGM: ', pbwg.AnimatedMarker(), ' ',
               pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        sph_code = igm.setup_sph_code(Gadget2, N1, N2, L, rho, u)

    gal.merger_and_igm(galaxy1, galaxy2, converter, sph_code, n_halo, t_end, SCRIPT_PATH, plot=PLOT)

