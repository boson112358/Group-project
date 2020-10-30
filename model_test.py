###### importing modules ######

import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '' and os.name != 'nt':
    print('No display found. Using non-interactive Agg backend for matplotlib')
    mpl.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import inspect
import sys

#finds script path
SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))

#creates plots folder
plot_folder = SCRIPT_PATH + '/plots/testrun/'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)
    
#creates csv folder
csv_folder = SCRIPT_PATH + '/data/galaxy-csv/'
if not os.path.exists(csv_folder):
    os.makedirs(csv_folder)
    
#ignore warnings (AmuseWarning at end of simulation)
import warnings
warnings.filterwarnings('ignore')

#defines parser for terminal usage
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--nosimulation', 
                    help='Skip galaxy simulation', 
                    action='store_true')
#parser.add_argument('-f', 
#                    help='Foo value, defined for jupyter notebook compatibility')
args = parser.parse_args()

NOSIMULATION = args.nosimulation

if NOSIMULATION:
    print('Skipping galaxy simulation', flush=True)
else:
    print('Simuation test run', flush=True)

#progressbar import
import progressbar as pbar
import progressbar.widgets as pbwg

from amuse.lab import units, Particle, read_set_from_file, nbody_system

#importing data analysis functions
import data_analysis as da
    
###### initial conditions ######

#simulation parameters

t_end = 10 | units.Myr
t_step = 0.1 | units.Myr

#M31 displacement
rotation = np.array([[0.7703,  0.3244,  0.5490],
                     [-0.6321, 0.5017,  0.5905],
                     [-0.0839, -0.8019, 0.5915]])

traslation = [-379.2, 612.7, 283.1] | units.kpc
radial_velocity = 117 * np.array([0.4898, -0.7914, 0.3657]) | units.kms
transverse_velocity = 50 * np.array([0.5236, 0.6024, 0.6024]) | units.kms

#galaxy parameters

scale_mass_galaxy = 1.0e12 | units.MSun
scale_radius_galaxy = 100 | units.kpc
n_bulge = 10000
n_disk = 20000
n_halo = 40000

mw_parameters = {'name': 'mw',
                 #halo parameters
                 'n_halo': n_halo,
                 'halo_scale_length': 12.96 | units.kpc,
                 #disk parameters
                 'disk_number_of_particles' : n_disk,
                 'disk_mass' : 19.66 * 2.33 * 10e9 | units.MSun,
                 'disk_scale_length' : 2.806 | units.kpc,
                 'disk_outer_radius' : 30 | units.kpc, 
                 'disk_scale_height_sech2' : 0.409 | units.kpc,
                 'disk_central_radial_velocity_dispersion': 0.7,
                 #bulge parameters
                 'bulge_scale_radius' : 0.788 | units.kpc,
                 'bulge_number_of_particles' : n_bulge,
                 #unused parameters (unclear effect)
                 "halo_streaming_fraction": 0.,
                 "bulge_streaming_fraction": 0.,
                 #unused parameters (they cause problems)
                 'disk_scale_length_of_sigR2': 2.806 | units.kpc}


###### loading galaxy ######

converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)

mw_data_path = SCRIPT_PATH + '/galaxies/data/testrun/{}_test'.format(mw_parameters['name'])

widgets = ['Found galaxy data, loading: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
    mw = read_set_from_file(mw_data_path, "hdf5")

#print('Correcting particles velocities ...', flush=True)
#mw[n_disk:n_disk+n_bulge].velocity = mw[n_disk:n_disk+n_bulge].velocity/10000


###### model analysis ######

mw_total_mass = da.galaxy_total_mass(mw)
mw_halo, mw_disk, mw_bulge = da.galaxy_structures(mw, n_disk, n_bulge)
print('Halo particles: {}; expected: {}'.format(len(mw_halo), n_halo)) 
print('Disk particles: {}; expected: {}'.format(len(mw_disk), n_disk))
print('Bulge particles: {}; expected: {}'.format(len(mw_bulge), n_bulge))

widgets = ['Plotting galaxy structures: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
    for structure, filename in zip([mw_halo, mw_disk, mw_bulge], ['mw_halo', 'mw_disk', 'mw_bulge']):
        da.plot_galaxy_structure(structure, filename)
        
halo_df_path = csv_folder + 'mw_halo_velocities.csv'
disk_df_path = csv_folder + 'mw_disk_velocities.csv'
bulge_df_path = csv_folder + 'mw_bulge_velocities.csv'

if os.path.exists(halo_df_path) and os.path.exists(disk_df_path) and os.path.exists(bulge_df_path):
    widgets = ['Found galaxy velocities data, loading: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        halo_df = pd.read_csv(halo_df_path)
        disk_df = pd.read_csv(disk_df_path)
        bulge_df = pd.read_csv(bulge_df_path)
else:
    widgets = ['Computing galaxy velocities: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        halo_dict, disk_dict, bulge_dict = da.galaxy_structure_velocity(mw, n_disk, n_bulge)
        halo_df = pd.DataFrame(halo_dict)
        disk_df = pd.DataFrame(disk_dict)
        bulge_df = pd.DataFrame(bulge_dict)
    widgets = ['Saving logs: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        halo_df.to_csv('{}{}.{}'.format(csv_folder, 'mw_halo_velocities', 'csv'), index=False)
        disk_df.to_csv('{}{}.{}'.format(csv_folder, 'mw_disk_velocities', 'csv'), index=False)
        bulge_df.to_csv('{}{}.{}'.format(csv_folder, 'mw_bulge_velocities', 'csv'), index=False)
        
widgets = ['Plotting bulge velocities: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
    com_distance = bulge_df['com_distance'].values    #deprecated in pandas 0.24.1 which I cannot install on STRW computer
    vel_components = [bulge_df['angular_velocity'], bulge_df['radial_velocity'], 
                      bulge_df['tangential_velocity'], bulge_df['total_velocity']]
    plot_prefix = 'mw_bulge_'
    for component, filename in zip(vel_components, ['angular_velocity', 'radial_velocity', 'tangential_velocity', 'total_velocity']):
        da.plot_velocity_component(com_distance, component, plot_prefix + filename)

vel_factor = 1/6.5
        
widgets = ['Plotting galaxy rotation curve (1/6.5 of v): ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
    com_distance = [#halo_df['com_distance'].values, 
                    disk_df['com_distance'].values * vel_factor, 
                    bulge_df['com_distance'].values * vel_factor]
    total_vel = [#halo_df['total_velocity'].values, 
                 disk_df['total_velocity'].values * vel_factor, 
                 bulge_df['total_velocity'].values * vel_factor]
    
    da.galaxy_rotation_curve(com_distance, total_vel, 'mw_rotation_curve', labels=['disk', 'bulge'])

###### running simulation ######

if not NOSIMULATION:
    import galaxies as gal
    print('Saving plots in: {}'.format(gal.__SCRIPT_PATH__), flush=True)
    
    t_end_int = int(np.round(t_end.value_in(units.Myr), decimals=0))
    t_step_int = int(np.round(t_step.value_in(units.Myr), decimals=0))
    print('Simulating mw (t = {} Myr, step = {} Myr) ...'.format(t_end_int, t_step_int), flush=True)
    gal.simulate_single_galaxy(mw, converter, n_halo, n_bulge, n_disk, t_end, interval=t_step, plot=True)