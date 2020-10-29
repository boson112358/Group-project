###### importing modules ######

import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '' and os.name != 'nt':
    print('No display found. Using non-interactive Agg backend for matplotlib')
    mpl.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

import inspect
import sys

#finds script path
#filename = inspect.getframeinfo(inspect.currentframe()).filename
SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
#add progressbar path to possible paths to import modules
sys.path.insert(1, SCRIPT_PATH + '/python-progressbar/')

#creates plots folder
plot_folder = SCRIPT_PATH + '/plots/testrun/'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)
    
#ignore warnings (AmuseWarning at end of simulation)

import warnings
warnings.filterwarnings('ignore')

#progressbar import
import progressbar as pbar
import progressbar.widgets as pbwg

from amuse.lab import *

    
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

print('Simuation test run', flush=True)

converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)

mw_data_path = SCRIPT_PATH + '/galaxies/data/testrun/{}_test'.format(mw_parameters['name'])

widgets = ['Found galaxies data, loading: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
    mw = read_set_from_file(mw_data_path, "hdf5")
    
mw_halo = mw[n_disk+n_bulge:]
mw_bulge = mw[n_disk:n_disk+n_bulge]
mw_disk = mw[:n_disk]

print('Correcting particles velocities ...', flush=True)
mw[n_disk:n_disk+n_bulge].velocity = mw[n_disk:n_disk+n_bulge].velocity/10000
    
import galaxies as gal
print('Saving plots in: {}'.format(gal.__SCRIPT_PATH__), flush=True)

###### running simulation ######

t_end_int = int(np.round(t_end.value_in(units.Myr), decimals=0))
t_step_int = int(np.round(t_step.value_in(units.Myr), decimals=0))
print('Simulating mw (t = {} Myr, step = {} Myr) ...'.format(t_end_int, t_step_int), flush=True)
gal.simulate_single_galaxy(mw, converter, n_halo, n_bulge, n_disk, t_end, interval=t_step, plot=True)