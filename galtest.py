import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '' and os.name != 'nt':
    print('No display found. Using non-interactive Agg backend for matplotlib')
    mpl.use('Agg')


import numpy as np
import matplotlib.pyplot as plt
import random
import math
from scipy.stats import kde

import inspect
import sys

#finds script path
filename = inspect.getframeinfo(inspect.currentframe()).filename
SCRIPT_PATH = os.path.dirname(os.path.abspath(filename))
#add progressbar path to possible paths to import modules
sys.path.insert(1, SCRIPT_PATH + '/python-progressbar/')
sys.path.insert(1, SCRIPT_PATH + '/galaxies/')

from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model

#progressbar import
import progressbar as pbar
import progressbar.widgets as pbwg

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
    
    savepath = SCRIPT_PATH + '/plots/single_galaxy_test/'
    
    plt.savefig(savepath + filename)
    

def density_histogram(galaxy, title, filename):
    x = galaxy.x.value_in(units.kpc)
    y = galaxy.y.value_in(units.kpc)
    nbins = 5000
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    x_label = "X [kpc]"
    y_label = "Y [kpc]"

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-300, 300)
    plt.ylim(-300, 300)
   
    ax.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Greens_r)
    #plt.colorbar()
    #cbar.set_label('Density')
    plt.tight_layout()
    
    savepath = SCRIPT_PATH + '/plots/single_galaxy_test/'
    
    plt.savefig(savepath + filename)
    

def density_histogram2(galaxy, title, filename):
    x = galaxy.x.value_in(units.kpc)
    y = galaxy.y.value_in(units.kpc)
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    xy_range = [[-20, 20], [-20, 20]]
    
    counts,ybins,xbins,image = plt.hist2d(x, y, bins=200, range=xy_range,
                                          #normed=True
                                          )
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    cont = ax.contourf(counts, extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()], cmap='RdGy')
    
    cbar = plt.colorbar(cont)
    cbar.set_label('Density')
    plt.tight_layout()
    
    savepath = SCRIPT_PATH + '/plots/single_galaxy_test/'
    
    plt.savefig(savepath + filename)
    

def make_galaxy(converter, galaxy_dict, script_path):
    
    n_halo = galaxy_dict['n_halo']
    
    widgets = ['Building {} galaxy: '.format(galaxy_dict['name']), pbwg.AnimatedMarker(), ' ',
               pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy = new_galactics_model(n_halo,
                                     converter,
                                     #halo_scale_length = galaxy_dict['halo_scale_length'],
                                     disk_number_of_particles = galaxy_dict['disk_number_of_particles'],
                                     disk_mass = galaxy_dict['disk_mass'],
                                     disk_scale_length = galaxy_dict['disk_scale_length'],
                                     disk_outer_radius = galaxy_dict['disk_outer_radius'],
                                     disk_truncation_width = galaxy_dict['disk_truncation_width'],
                                     disk_scale_height_sech2 = galaxy_dict['disk_scale_height_sech2'],
                                     disk_central_radial_velocity_dispersion = galaxy_dict['disk_central_radial_velocity_dispersion'],
                                     disk_scale_length_of_sigR2 = galaxy_dict['disk_scale_length_of_sigR2'],
                                     bulge_number_of_particles = galaxy_dict['bulge_number_of_particles'],
                                     #bulge_scale_radius = galaxy_dict['bulge_scale_radius']
                                    )
    
    galaxy_data_path = script_path + '/galaxies/data/testrun/{}_testrun'.format(galaxy_dict['name'])
    
    widgets = ['Saving {} galaxy data: '.format(galaxy_dict['name']), 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        write_set_to_file(galaxy, galaxy_data_path, 'hdf5')

    return galaxy
    

def simulate_single_galaxy(galaxy1, converter, n_halo, n_bulge, n_disk, t_end, plot=False):
    
    dynamics_code = Gadget2(converter, number_of_workers=4)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = dynamics_code.dm_particles.add_particles(galaxy1)
    halo1 = set1[n_disk+n_bulge:]
    bulge1 = set1[n_disk:n_disk+n_bulge]
    disk1 = set1[:n_disk]
    
    dynamics_code.particles.move_to_center()
    #dynamics_code.parameters.timestep = 0.5 | units.Myr
    
    if plot == True:
        plot_number = 0
        plot_single_galaxy(halo1, disk1, bulge1,
             "MW\nt = 0 Myr", 
             'mw_testrun_' + str(plot_number).zfill(4))
    
    current_iter = 0
    interval = 0.5 | units.Myr
    t_end_in_Myr = t_end.as_quantity_in(units.Gyr)
    total_iter = int(t_end_in_Myr/interval) + 1
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    while dynamics_code.model_time < t_end:
        
        current_iter +=1
        
        dynamics_code.evolve_model(dynamics_code.model_time + interval)
        
        if plot == True:
            if current_iter in [10*i for i in range(1, total_iter)]:
                plot_number += 1
                plot_single_galaxy(halo1, disk1, bulge1, 
                             "MW\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                                                                  decimals=0))),
                             'mw_testrun_' + str(plot_number).zfill(4))
        
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        plot_number += 1
        plot_single_galaxy(halo1, disk1, bulge1,
                     "MW\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),                              
                     'mw_testrun_' + str(plot_number).zfill(4))
        
    dynamics_code.stop()
    
#creates plots folder
mw_folder = SCRIPT_PATH + '/plots/single_galaxy_test/'
if not os.path.exists(mw_folder):
    os.makedirs(mw_folder)

#simulation parameters
scale_mass_galaxy = 10e12 | units.MSun
scale_radius_galaxy = 100 | units.kpc
t_end = 4000 | units.Myr
n_bulge = 20000
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
                 'disk_central_radial_velocity_dispersion': 0.7 * 100 | units.kms,
                 'disk_scale_length_of_sigR2': 2.806 | units.kpc,
                 'bulge_number_of_particles' : n_bulge,
                 'bulge_scale_radius' : 0.788 | units.kpc}

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


#reads galaxy data
galaxy_data = SCRIPT_PATH + '/galaxies/data/mw_full' 
converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)
#m31_not_displaced = read_set_from_file(galaxy_data, "hdf5")

#from galaxies import galaxies as gal

#m31 = gal.displace_galaxy(m31_not_displaced, rotation, traslation, radial_velocity, transverse_velocity)
#mw = make_galaxy(converter, mw_parameters, SCRIPT_PATH)

mw = read_set_from_file(galaxy_data, "hdf5")

print('Simulating m31 (t = {} Myr) ...'.format(int(np.round(t_end.value_in(units.Myr), decimals=0))), flush=True)
simulate_single_galaxy(mw, converter, n_halo, n_bulge, n_disk, t_end, plot=True)