#!/usr/bin/env python
# coding: utf-8

import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '' and os.name != 'nt':
    print('No display found. Using non-interactive Agg backend for matplotlib')
    mpl.use('Agg')


import numpy as np
import matplotlib.pyplot as plt

from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model

import inspect
import sys

#finds script path
filename = inspect.getframeinfo(inspect.currentframe()).filename
script_path = os.path.dirname(os.path.abspath(filename))
#add progressbar path to possible paths to import modules
sys.path.insert(1, script_path + '/python-progressbar/')


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


def make_plot(disk1, disk2, filename):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    
    fig = plt.figure()
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-300, 300)
    plt.ylim(-300, 300)

    plt.scatter(disk1.x.value_in(units.kpc), disk1.y.value_in(units.kpc),
                   c='tab:blue', alpha=1, s=1, lw=0)
    plt.scatter(disk2.x.value_in(units.kpc), disk2.y.value_in(units.kpc),
                   c='tab:orange', alpha=1, s=1, lw=0)
    
    savepath = script_path + '/plots/'
    
    plt.savefig(savepath + filename)

def make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk):
    converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    
    widgets = ['Building galaxy 1: ', pbwg.AnimatedMarker(), ' ',
               pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy1 = new_galactics_model(n_halo,
                                      converter,
                                      #do_scale = True,
                                      bulge_number_of_particles=n_bulge,
                                      disk_number_of_particles=n_disk)
    
    widgets = ['Building galaxy 2: ', pbwg.AnimatedMarker(),
               pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy2 = Particles(len(galaxy1))
        galaxy2.mass = galaxy1.mass
        galaxy2.position = galaxy1.position
        galaxy2.velocity = galaxy1.velocity
    
    widgets = ['Adjusting relative velocities and orientations: ', 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy1.rotate(0., np.pi/2, np.pi/4)
        galaxy1.position += [100.0, 100, 0] | units.kpc
        galaxy1.velocity += [-10.0, 0.0, -10.0] | units.km/units.s
        galaxy2.rotate(np.pi/4, np.pi/4, 0.0)
        galaxy2.position -= [100.0, 0, 0] | units.kpc
        galaxy2.velocity -= [0.0, 0.0, 0] | units.km/units.s

    return galaxy1, galaxy2, converter

def simulate_merger(galaxy1, galaxy2, converter, n_halo, t_end):
    converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    dynamics_code = Gadget2(converter, number_of_workers=4)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = dynamics_code.particles.add_particles(galaxy1)
    set2 = dynamics_code.particles.add_particles(galaxy2)
    
    dynamics_code.particles.move_to_center()
    
    disk1 = set1[:n_halo]
    disk2 = set2[:n_halo]
    
    if PLOT == True:
        make_plot(disk1, disk2, "Galaxy_merger_t0")
    
    current_iter = 0
    interval = 0.5 | units.Myr
    total_iter = int(t_end/interval) + 1
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    while dynamics_code.model_time < t_end:
        
        current_iter +=1
        
        dynamics_code.evolve_model(dynamics_code.model_time + interval)
                
        progress.update(current_iter)
        
    progress.finish()
    
    if PLOT == True:
        make_plot(disk1, disk2,
                  "Galaxy_merger_t" + str(t_end.value_in(units.Myr))+"Myr")
        
    dynamics_code.stop()


M_galaxy = 1.0e12 | units.MSun
R_galaxy = 10 | units.kpc
n_bulge = 10000
n_disk = 10000
n_halo = 20000
t_end = 200|units.Myr

galaxy1, galaxy2, converter = make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk)
simulate_merger(galaxy1, galaxy2, converter, n_halo, t_end)

