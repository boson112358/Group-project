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

from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model

#progressbar import
import progressbar as pbar
import progressbar.widgets as pbwg

def make_plot_mw(galaxy, disk, title, filename):
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

    ax.scatter(galaxy.x.value_in(units.kpc), galaxy.y.value_in(units.kpc),
                   c='tab:blue', alpha=1, s=1, lw=0)
    ax.scatter(disk.x.value_in(units.kpc), disk.y.value_in(units.kpc),
                   c='tab:orange', alpha=1, s=1, lw=0)
    
    savepath = SCRIPT_PATH + '/plots/mw_new_model/'
    
    plt.savefig(savepath + filename)

def make_mw_test(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk):
    converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    
    widgets = ['Building galaxy 1: ', pbwg.AnimatedMarker(), ' ',
               pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy1 = new_galactics_model(n_halo,
                                      converter,
                                      halo_outer_radius = 244.48999 | units.kpc,
                                      disk_number_of_particles = n_disk,
                                      disk_mass = 19.66 * 2.33 * 10e9 | units.MSun,
                                      disk_scale_length = 2.806 | units.kpc,
                                      disk_outer_radius = 30 | units.kpc, 
                                      disk_scale_height_sech2 = 0.409 | units.kpc,
                                      bulge_scale_radius = 0.788 | units.kpc,
                                      bulge_number_of_particles = n_bulge)
                                      #central_radial_vel_dispersion = 0.7 * 100 | units.kms,
                                      #scale_length_of_sigR2 = 2.806 | units.kpc)
    
    
    widgets = ['Saving galaxies data: ', 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        write_set_to_file(galaxy1, SCRIPT_PATH + '/galaxies/data/MW_new_model', 'hdf5')

    return galaxy1, converter

def simulate_mw(galaxy1, converter, n_halo, t_end, plot=False):
    #converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    dynamics_code = Gadget2(converter, number_of_workers=4)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = dynamics_code.dm_particles.add_particles(galaxy1)
    disk1 = set1[:n_halo]
    
    dynamics_code.particles.move_to_center()
    #dynamics_code.parameters.timestep = 0.5 | units.Myr
    
    if plot == True:
        plot_number = 0
        make_plot_mw(set1, disk1,
             "MW\nt = 0 Myr", 
             'mw_' + str(plot_number).zfill(4))
    
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
                make_plot_mw(set1, disk1, 
                             "MW\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                                                                  decimals=0))),
                             'mw_' + str(plot_number).zfill(4))
        
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        plot_number += 1
        make_plot_mw(set1, disk1, 
                     "MW\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),                              
                     'mw_' + str(plot_number).zfill(4))
        
    dynamics_code.stop()

M_galaxy = 1.0e12 | units.MSun
R_galaxy = 10 | units.kpc
n_bulge = 10000
n_disk = 10000
n_halo = 20000
t_end = 1100 | units.Myr

mw_folder = SCRIPT_PATH + '/plots/mw_new_model/'

if not os.path.exists(mw_folder):
    os.makedirs(mw_folder)

mw_data = SCRIPT_PATH + '/galaxies/data/MW_new_model' 

if not os.path.exists(mw_data):
    MW, converter = make_mw_test(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk)

if os.path.exists(mw_data):
    widgets = ['Loading galaxies data: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
        MW = read_set_from_file(mw_data, "hdf5")

simulate_mw(MW, converter, n_halo, t_end, plot=True)