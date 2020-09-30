#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model

from progressbar import ProgressBar


def make_plot(disk1, disk2, filename):
    #print('Plotting {} ... '.format(filename), sep=' ', end='')
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
    
    plt.savefig(filename)
    #print('Done\n')

def make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk):
    converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    
    print('Building galaxy 1 ... ', sep=' ', end='')
    galaxy1 = new_galactics_model(n_halo,
                                  converter,
                                  #do_scale = True,
                                  bulge_number_of_particles=n_bulge,
                                  disk_number_of_particles=n_disk)
    print('Done\n')
    
    print('Building galaxy 2 ... ', sep=' ', end='')
    galaxy2 = Particles(len(galaxy1))
    galaxy2.mass = galaxy1.mass
    galaxy2.position = galaxy1.position
    galaxy2.velocity = galaxy1.velocity
    print('Done\n')
    
    print('Adjusting relative velocities and orientations ... ', sep=' ', end='')
    galaxy1.rotate(0., np.pi/2, np.pi/4)
    galaxy1.position += [100.0, 100, 0] | units.kpc
    #galaxy1.velocity += [-3000.0, 0.0, -3000.0] | units.km/units.s
    galaxy1.velocity += [-10.0, 0.0, -10.0] | units.km/units.s

    galaxy2.rotate(np.pi/4, np.pi/4, 0.0)
    galaxy2.position -= [100.0, 0, 0] | units.kpc
    galaxy2.velocity -= [0.0, 0.0, 0] | units.km/units.s
    print('Done\n')

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
    
    make_plot(disk1, disk2, "Galaxy_merger_t0")
    
    current_iter = 0
    interval = 0.5 | units.Myr
    total_iter = t_end/interval
    
    progress = ProgressBar(total_days, comp_line='Done', 
                           global_time_measure=True, iteration_time_measure=True)
    progress.show(current_day, prefix='Step {} of {}:'.format(current_iter, total_iter), suffix=' ... ')
    
    while dynamics_code.model_time < t_end:
        
        current_iter +=1
        
        progress.start_iteration_measure()
        
        dynamics_code.evolve_model(dynamics_code.model_time + interval)
                
        progress.show(current_day, prefix='Step {} of {}:'.format(current_iter, total_iter), suffix=' ... ')
        
    dynamics.stop()
    
    make_plot(disk1, disk2,
              "Galaxy_merger_t" + str(t_end.value_in(units.Myr))+"Myr")


M_galaxy = 1.0e12 | units.MSun
R_galaxy = 10 | units.kpc
n_bulge = 10000
n_disk = 10000
n_halo = 20000
t_end = 200|units.Myr

galaxy1, galaxy2, converter = make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk)
simulate_merger(galaxy1, galaxy2, converter, n_halo, t_end)


# Error in the terminal:
# 
# There are not enough slots available in the system to satisfy the 4
# slots that were requested by the application:
# 
#   /home/al/anaconda3/envs/amuse-rp/bin/python
# 
# Either request fewer slots for your application, or make more slots
# available for use.
# 
# A "slot" is the Open MPI term for an allocatable unit where we can
# launch a process.  The number of slots available are defined by the
# environment in which Open MPI processes are run:
# 
#   1. Hostfile, via "slots=N" clauses (N defaults to number of
#      processor cores if not provided)
#   2. The --host command line parameter, via a ":N" suffix on the
#      hostname (N defaults to 1 if not provided)
#   3. Resource manager (e.g., SLURM, PBS/Torque, LSF, etc.)
#   4. If none of a hostfile, the --host command line parameter, or an
#      RM is present, Open MPI defaults to the number of processor cores
# 
# In all the above cases, if you want Open MPI to default to the number
# of hardware threads instead of the number of processor cores, use the
# --use-hwthread-cpus option.
# 
# Alternatively, you can use the --oversubscribe option to ignore the
# number of available slots when deciding the number of processes to
# launch.



