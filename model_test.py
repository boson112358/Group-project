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

"""
#creates model folder
glxy_folder = SCRIPT_PATH + '/data/model/'
if not os.path.exists(glxy_folder):
    os.makedirs(glxy_folder)

#creates plots folder
plot_folder = SCRIPT_PATH + '/data/test_plots/'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)
    
#creates csv folder
csv_folder = SCRIPT_PATH + '/data/test_csv/'
if not os.path.exists(csv_folder):
    os.makedirs(csv_folder)
    
#creates plots folder
sim_plot_folder = SCRIPT_PATH + '/data/sim_plots/'
if not os.path.exists(sim_plot_folder):
    os.makedirs(sim_plot_folder)
"""
    
#ignore warnings (AmuseWarning at end of simulation)
import warnings
#warnings.filterwarnings('ignore')

#defines parser for terminal usage
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--generation', 
                    help='Generate a new galaxy model', 
                    action='store_true')
parser.add_argument('--mway', 
                    help='Using Milky Way data', 
                    action='store_true')
parser.add_argument('--andromeda', 
                    help='Using Andromeda data', 
                    action='store_true')
parser.add_argument('--simulation', 
                    help='Run galaxy simulation', 
                    action='store_true')
parser.add_argument('--analysis', 
                    help='Analyse galaxy model', 
                    action='store_true')
parser.add_argument('--correction', 
                    help='Correct galaxy velocity and mass', 
                    action='store_true')
args = parser.parse_args()

GENERATION = args.generation
MWAY = args.mway
ANDROMEDA = args.andromeda
SIMULATION = args.simulation
ANALYSIS = args.analysis
CORRECTION = args.correction

from amuse.lab import units, Particles, nbody_system
from amuse.ext.galactics_model import new_galactics_model

#importing data analysis functions
import modules.data_analysis as da
import modules.galaxies as gal
import modules.simulations as sim
import modules.progressbar as pbar
import modules.progressbar.widgets as pbwg


###### initial conditions ######

#simulation parameters
t_end = 1000 | units.Myr
t_step = 0.1 | units.Myr

#GalacICs output dir
output_directory = '/data1/brentegani/'

#M31 displacement
rotation = np.array([[0.7703,  0.3244,  0.5490],
                     [-0.6321, 0.5017,  0.5905],
                     [-0.0839, -0.8019, 0.5915]])

traslation = [-379.2, 612.7, 283.1] | units.kpc
radial_velocity = 117 * np.array([0.4898, -0.7914, 0.3657]) | units.kms
transverse_velocity = 50 * np.array([0.5236, 0.6024, 0.6024]) | units.kms

#galaxy parameters
scale_mass_galaxy = 1e12 | units.MSun
scale_radius_galaxy = 80 | units.kpc
n_bulge = 10000
n_disk = 20000
n_halo = 40000

mw_parameters = {'name': 'mw',
                 #halo parameters
                 'n_halo': n_halo,
                 'halo_scale_radius': 12.96 | units.kpc,
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

m31_parameters = {'name': 'm31_not_displaced',
                  #halo parameters
                  'n_halo': n_halo,
                  'halo_scale_radius': 12.94 | units.kpc,
                  #disk parameters
                  'disk_number_of_particles' : n_disk,
                  'disk_mass' : 33.40 * 2.33 * 10e9 | units.MSun,
                  'disk_scale_length' : 5.577 | units.kpc,
                  'disk_outer_radius' : 30 | units.kpc, 
                  'disk_scale_height_sech2' : 0.3 | units.kpc,
                  'disk_central_radial_velocity_dispersion': 0.7,
                  #bulge parameters
                  'bulge_scale_radius' : 1.826 | units.kpc,
                  'bulge_number_of_particles' : n_bulge,
                  #unused parameters (unclear effect)
                  "halo_streaming_fraction": 0.,
                  "bulge_streaming_fraction": 0.,
                  #unused parameters (they cause problems)
                  'disk_scale_length_of_sigR2': 5.577 | units.kpc}


###### selecting galaxy ######

if MWAY:
    print('Using Milky Way data', flush=True)
    glxy_param = mw_parameters
elif ANDROMEDA:
    print('Using Andromeda data', flush=True)
    glxy_param = m31_parameters
else:
    raise ValueError('Choose a galaxy')

if SIMULATION:
    print('Simuation test run', flush=True)
else:
    print('Skipping galaxy simulation', flush=True)
    
if ANALYSIS:
    print('Analysis of galaxy model', flush=True)
else:
    print('Skipping model analysis', flush=True)


###### loading galaxy ######

glxy_name = glxy_param['name']

converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)

if GENERATION:
    glxy, glxy_path = gal.make_galaxy(glxy_param['n_halo'], converter, glxy_param['name'], test=True,
                                      #output dir
                                      output_directory = '/data1/brentegani/',
                                      #halo parameters
                                      halo_scale_radius = glxy_param['halo_scale_radius'],
                                      #disk parameters
                                      disk_number_of_particles = glxy_param['disk_number_of_particles'],
                                      disk_mass = glxy_param['disk_mass'],
                                      disk_scale_length = glxy_param['disk_scale_length'],
                                      disk_outer_radius = glxy_param['disk_outer_radius'],
                                      disk_scale_height_sech2 = glxy_param['disk_scale_height_sech2'],
                                      disk_central_radial_velocity_dispersion=glxy_param['disk_central_radial_velocity_dispersion'],
                                      #bulge paramaters
                                      bulge_scale_radius = glxy_param['bulge_scale_radius'],
                                      bulge_number_of_particles = glxy_param['bulge_number_of_particles'])

else:
    glxy, glxy_path = gal.load_galaxy_data(glxy_param['name'], test=True)


###### model analysis ######

if ANALYSIS:
    glxy_halo, glxy_disk, glxy_bulge = da.galaxy_structures(glxy, n_disk, n_bulge)
    print('Halo particles: {}; expected: {}'.format(len(glxy_halo), n_halo)) 
    print('Disk particles: {}; expected: {}'.format(len(glxy_disk), n_disk))
    print('Bulge particles: {}; expected: {}'.format(len(glxy_bulge), n_bulge))

    widgets = ['Plotting galaxy structures: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        for structure, filename in zip([glxy_halo, glxy_disk, glxy_bulge], [glxy_name + '_halo', 
                                                                            glxy_name + '_disk', 
                                                                            glxy_name + '_bulge']):
            da.plot_galaxy_structure(structure, glxy_path, filename)
            
    _, csv_folder = da.create_model_analysis_dirs(glxy_path)

    halo_df_path = csv_folder + glxy_name + '_halo_velocities.csv'
    disk_df_path = csv_folder + glxy_name + '_disk_velocities.csv'
    bulge_df_path = csv_folder + glxy_name + '_bulge_velocities.csv'

    if os.path.exists(halo_df_path) and os.path.exists(disk_df_path) and os.path.exists(bulge_df_path):
        widgets = ['Found galaxy velocities data, loading: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
        with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
            halo_df = pd.read_csv(halo_df_path)
            disk_df = pd.read_csv(disk_df_path)
            bulge_df = pd.read_csv(bulge_df_path)
    else:
        widgets = ['Computing galaxy velocities: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
        with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
            halo_dict, disk_dict, bulge_dict = da.galaxy_structure_velocity(glxy, n_disk, n_bulge)
            halo_df = pd.DataFrame(halo_dict)
            disk_df = pd.DataFrame(disk_dict)
            bulge_df = pd.DataFrame(bulge_dict)
        widgets = ['Saving logs: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
        with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
            halo_df.to_csv('{}{}.{}'.format(csv_folder, glxy_name + '_halo_velocities', 'csv'), index=False)
            disk_df.to_csv('{}{}.{}'.format(csv_folder, glxy_name + '_disk_velocities', 'csv'), index=False)
            bulge_df.to_csv('{}{}.{}'.format(csv_folder, glxy_name + '_bulge_velocities', 'csv'), index=False)
            
            
###### correcting velocities and mass ######

corr = ''

if CORRECTION:
    print('Correcting velocities and mass: ...', flush=True)
    vel_factor = 1/100
    mass_factor = 1/1000

    start_vel_ex = glxy.velocity[800]
    glxy.velocity = glxy.velocity * vel_factor
    print('v_ini = {}\nv_fin = {}'.format(start_vel_ex, glxy.velocity[800]), flush=True)
    
    if ANALYSIS:
        start_vel_ex2 = disk_df['angular_velocity'][800]
        for key in ['angular_velocity', 'radial_velocity', 'tangential_velocity', 'total_velocity']:
            halo_df[key] *= vel_factor
            disk_df[key] *= vel_factor
            bulge_df[key] *= vel_factor
        print('v_ini = {}\nv_fin = {}'.format(start_vel_ex2, disk_df['angular_velocity'][800]), flush=True)

    start_mass = da.galaxy_total_mass(glxy)
    glxy.mass = glxy.mass * mass_factor
    final_mass = da.galaxy_total_mass(glxy)

    corr = '_corr'
        

###### plotting velocity components ######

if ANALYSIS:
    widgets = ['Plotting halo velocities: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        com_distance = halo_df['com_distance'].values    #deprecated in pandas 0.24.1 which I cannot install on STRW computer
        vel_components = [halo_df['angular_velocity'], halo_df['radial_velocity'], 
                          halo_df['tangential_velocity'], halo_df['total_velocity']]
        plot_prefix = glxy_name + corr + '_halo_'
        for component, filename in zip(vel_components, ['angular_velocity', 'radial_velocity', 
                                                        'tangential_velocity', 'total_velocity']):
            da.plot_velocity_component(com_distance, component, glxy_path, plot_prefix + filename)
            
    widgets = ['Plotting disk velocities: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        com_distance = disk_df['com_distance'].values    #deprecated in pandas 0.24.1 which I cannot install on STRW computer
        vel_components = [disk_df['angular_velocity'], disk_df['radial_velocity'], 
                          disk_df['tangential_velocity'], disk_df['total_velocity']]
        plot_prefix = glxy_name + corr + '_disk_'
        for component, filename in zip(vel_components, ['angular_velocity', 'radial_velocity', 
                                                        'tangential_velocity', 'total_velocity']):
            da.plot_velocity_component(com_distance, component, glxy_path, plot_prefix + filename)

    widgets = ['Plotting bulge velocities: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        com_distance = bulge_df['com_distance'].values    #deprecated in pandas 0.24.1 which I cannot install on STRW computer
        vel_components = [bulge_df['angular_velocity'], bulge_df['radial_velocity'], 
                          bulge_df['tangential_velocity'], bulge_df['total_velocity']]
        plot_prefix = glxy_name + corr + '_bulge_'
        for component, filename in zip(vel_components, ['angular_velocity', 'radial_velocity', 
                                                        'tangential_velocity', 'total_velocity']):
            da.plot_velocity_component(com_distance, component, glxy_path, plot_prefix + filename)

    widgets = ['Plotting galaxy rotation curve: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        com_distance = [halo_df['com_distance'].values, 
                        disk_df['com_distance'].values, 
                        bulge_df['com_distance'].values]
        com_distance_full = np.concatenate((halo_df['com_distance'].values, 
                                            disk_df['com_distance'].values, 
                                            bulge_df['com_distance'].values), axis=0)
        
        total_vel = [halo_df['total_velocity'].values, 
                     disk_df['total_velocity'].values, 
                     bulge_df['total_velocity'].values]
        total_vel_full = np.concatenate((halo_df['total_velocity'].values, 
                                         disk_df['total_velocity'].values, 
                                         bulge_df['total_velocity'].values), axis=0)
        
        da.galaxy_rotation_curve(com_distance_full, total_vel_full, glxy_path, 
                                 glxy_name + corr +'_rotation_curve', labels=glxy_param['name'])
        

###### running simulation ######

if SIMULATION:
    print('Saving plots in: {}'.format(sim.create_single_gal_dir(glxy_path)), flush=True)
    
    t_end_int = int(np.round(t_end.value_in(units.Myr), decimals=0))
    dec_num = 0
    if t_step.value_in(units.Myr) < 1:
        dec_num = 2
    t_step_int = np.round(t_step.value_in(units.Myr), decimals=dec_num)
    
    print('Simulating mw (t = {} Myr, step = {} Myr) ...'.format(t_end_int, t_step_int), flush=True)
    sim.simulate_single_galaxy(glxy, converter, n_halo, n_bulge, n_disk, t_end, glxy_path, 
                               interval=t_step, plot=True, plot_freq=1000)
