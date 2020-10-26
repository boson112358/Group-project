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

def plot_single_galaxy(halo, disk, title, filename):
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

    ax.scatter(halo.x.value_in(units.kpc), halo.y.value_in(units.kpc),
                   c='tab:blue', alpha=1, s=1, lw=0)
    ax.scatter(disk.x.value_in(units.kpc), disk.y.value_in(units.kpc),
                   c='tab:orange', alpha=1, s=1, lw=0)
    
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

def simulate_single_galaxy(galaxy1, converter, n_halo, t_end, plot=False):
    
    dynamics_code = Gadget2(converter, number_of_workers=4)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = dynamics_code.dm_particles.add_particles(galaxy1)
    halo1 = set1[n_halo:]
    disk1 = set1[:n_halo]
    
    dynamics_code.particles.move_to_center()
    #dynamics_code.parameters.timestep = 0.5 | units.Myr
    
    if plot == True:
        plot_number = 0
        density_histogram2(set1,
             "M31\nt = 0 Myr", 
             'm31_' + str(plot_number).zfill(4))
        #plot_single_galaxy(halo1, disk1,
        #     "M31\nt = 0 Myr", 
        #     'm31_' + str(plot_number).zfill(4))
    
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
                density_histogram2(set1, 
                             "M31\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                                                                  decimals=0))),
                             'm31_' + str(plot_number).zfill(4))
                #plot_single_galaxy(halo1, disk1, 
                #             "M31\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                #                                                  decimals=0))),
                #             'm31_' + str(plot_number).zfill(4))
        
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        plot_number += 1
        density_histogram2(set1,
                     "M31\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),                              
                     'm31_' + str(plot_number).zfill(4))
        #make_plot_mw(halo1, disk1, 
        #             "M31\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),                              
        #             'm31_' + str(plot_number).zfill(4))
        
    dynamics_code.stop()
    
#creates plots folder
mw_folder = SCRIPT_PATH + '/plots/single_galaxy_test/'
if not os.path.exists(mw_folder):
    os.makedirs(mw_folder)

#simulation parameters
scale_mass_galaxy = 1.0e12 | units.MSun
scale_radius_galaxy = 10 | units.kpc
t_end = 400 | units.Myr
n_halo = 20000

#M31 displacement
rotation = np.array([[0.7703,  0.3244,  0.5490],
                     [-0.6321, 0.5017,  0.5905],
                     [-0.0839, -0.8019, 0.5915]])

traslation = [-379.2, 612.7, 283.1] | units.kpc
radial_velocity = 117 * np.array([0.4898, -0.7914, 0.3657]) | units.kms
transverse_velocity = 50 * np.array([0.5236, 0.6024, 0.6024]) | units.kms

#reads galaxy data
galaxy_data = SCRIPT_PATH + '/galaxies/data/m31_not_displaced_full' 
converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)
m31_not_displaced = read_set_from_file(galaxy_data, "hdf5")

from galaxies import galaxies as gal

m31 = gal.displace_galaxy(m31_not_displaced, rotation, traslation, radial_velocity, transverse_velocity)

print('Simulating m31 ...', flush=True)
simulate_single_galaxy(m31, converter, n_halo, t_end, plot=True)