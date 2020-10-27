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
filename = inspect.getframeinfo(inspect.currentframe()).filename
SCRIPT_PATH = os.path.dirname(os.path.abspath(filename))
#add progressbar path to possible paths to import modules
sys.path.insert(1, '/home/brentegani/Desktop/sma-group-e-project/' + '/python-progressbar/')

#creates plots folder
mw_folder = SCRIPT_PATH + '/plots_test/'
if not os.path.exists(mw_folder):
    os.makedirs(mw_folder)

#progressbar import
import progressbar as pbar
import progressbar.widgets as pbwg

from amuse.lab import units, nbody_system

###### plot function ######

def plot_single_galaxy(halo, disk, bulge, title, filename):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)

    ax.scatter(halo.x.value_in(units.kpc), halo.y.value_in(units.kpc),
                   c='tab:blue', alpha=1, s=1, lw=0, label='halo')
    ax.scatter(disk.x.value_in(units.kpc), disk.y.value_in(units.kpc),
                   c='tab:orange', alpha=1, s=1, lw=0, label='disk')
    ax.scatter(bulge.x.value_in(units.kpc), bulge.y.value_in(units.kpc),
                   c='tab:green', alpha=1, s=1, lw=0, label='bulge')
    
    plt.legend(loc='upper right')
    
    savepath = SCRIPT_PATH + '/plots_test/'
    
    plt.savefig(savepath + filename)
    
###### initial conditions ######

#simulation parameters
scale_mass_galaxy = 1.0e12 | units.MSun
scale_radius_galaxy = 100 | units.kpc
t_end = 1000 | units.Myr
n_bulge = 10000
n_disk = 20000
n_halo = 40000

#M31 displacement
rotation = np.array([[0.7703,  0.3244,  0.5490],
                     [-0.6321, 0.5017,  0.5905],
                     [-0.0839, -0.8019, 0.5915]])

traslation = [-379.2, 612.7, 283.1] | units.kpc
radial_velocity = 117 * np.array([0.4898, -0.7914, 0.3657]) | units.kms
transverse_velocity = 50 * np.array([0.5236, 0.6024, 0.6024]) | units.kms

mw_parameters = {'name': 'mw',
                 'n_halo': n_halo,
                 #'halo_outer_radius' : 244.48999 | units.kpc,
                 'halo_scale_length': 12.96 | units.kpc,
                 'disk_number_of_particles' : n_disk,
                 'disk_mass' : 19.66 * 2.33 * 10e9 | units.MSun,
                 'disk_scale_length' : 2.806 | units.kpc,
                 'disk_outer_radius' : 30 | units.kpc,
                 'disk_truncation_width': 1.,
                 'disk_scale_height_sech2' : 0.409 | units.kpc,
                 'disk_central_radial_velocity_dispersion': 0.7, # * 100 | units.kms,
                 'disk_scale_length_of_sigR2': 2.806 | units.kpc,
                 'bulge_number_of_particles' : n_bulge,
                 'bulge_scale_radius' : 0.788 | units.kpc}

###### building galaxy ######

from amuse.ext.galactics_model import new_galactics_model

converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)

galaxy_dict = mw_parameters

#n_halo = galaxy_dict['n_halo']
    
widgets = ['Building {} galaxy: '.format(galaxy_dict['name']), pbwg.AnimatedMarker(), ' ',
           pbwg.Timer(), pbwg.EndMsg()]
with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
    galaxy = new_galactics_model(n_halo,
                                 converter,
                                 disk_number_of_particles = galaxy_dict['disk_number_of_particles'],
                                 disk_mass = galaxy_dict['disk_mass'],
                                 disk_scale_length = galaxy_dict['disk_scale_length'],
                                 disk_outer_radius = galaxy_dict['disk_outer_radius'],
                                 #disk_truncation_width = galaxy_dict['disk_truncation_width'],
                                 disk_scale_height_sech2 = galaxy_dict['disk_scale_height_sech2'],
                                 disk_central_radial_velocity_dispersion = galaxy_dict['disk_central_radial_velocity_dispersion'],
                                 disk_scale_length_of_sigR2 = galaxy_dict['disk_scale_length_of_sigR2'],
                                 bulge_number_of_particles = galaxy_dict['bulge_number_of_particles'],
                                 bulge_scale_radius = galaxy_dict['bulge_scale_radius'],
                                )