import numpy as np
import matplotlib.pyplot as plt
import random
import math
import sys

from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model
from amuse.couple import bridge

#progressbar import

import progressbar as pbar
import progressbar.widgets as pbwg


def make_plot(disk1, disk2, script_path, filename):
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
    
    
def make_plot_testdisk(disk1, disk2, test_disk, script_path, filename):
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
    plt.scatter(test_disk.x.value_in(units.kpc), test_disk.y.value_in(units.kpc),
                   c='tab:green', alpha=1, s=1, lw=0)
    
    savepath = script_path + '/plots/'
    
    plt.savefig(savepath + filename)
    
def make_plot_galstars(disk, stars, script_path, filename):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    
    fig = plt.figure()
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-300, 300)
    plt.ylim(-300, 300)

    plt.scatter(disk.x.value_in(units.kpc), disk.y.value_in(units.kpc),
                   c='tab:blue', alpha=1, s=1, lw=0)
    plt.scatter(stars.x.value_in(units.kpc), stars.y.value_in(units.kpc),
                   c='tab:orange', alpha=1, s=1, lw=0)
    
    savepath = script_path + '/plots/'
    
    plt.savefig(savepath + filename)


def nstars(numparticles, radius, radius_var, vel):  #I kinda messed up, I shouldn't include Sag so we can remove that part and rename
    particles = Particles(numparticles)
    i=0
    while i < numparticles:
        #here we generate our N bodies, we can change code in here to reflect different things
        #We can distribut them randomly, have them near eachother, all at the same radius etc.
        #All particles in x-y plane, this should still be adjusted accordingly
        
        #For now I will keep the mass and radius the same and simulate stars at our radius +- a bit
        angle = random.random()*2*math.pi
        particles[i].mass = 5 | units.MSun 
        particles[i].radius = 1.5 | units.RSun   
        adjusted_rad = radius*(1 + radius_var*(2*(random.random()-0.5)))  #this gives a random radius between 1-radius_var and 1+radius_var times the radius
        particles[i].position = (adjusted_rad*math.cos(angle),adjusted_rad*math.sin(angle) , 0 ) | units.kpc
        particles[i].velocity = [math.sin(angle)*vel,-math.cos(angle)*vel,0] | (units.m/units.s)
        
        i += 1
    
    return particles


def make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk, script_path, test=False):
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
    
    if not test:
        galaxy1_name = script_path + '/galaxies/data/M31_full'
        galaxy2_name = script_path + '/galaxies/data/MW_full'
    if test:
        galaxy1_name = script_path + '/galaxies/data/M31_test'
        galaxy2_name = script_path + '/galaxies/data/MW_test'
    
    widgets = ['Saving galaxies data: ', 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        write_set_to_file(galaxy1, galaxy1_name, 'hdf5')
        write_set_to_file(galaxy2, galaxy2_name, 'hdf5')

    return galaxy1, galaxy2, converter


def test_particles(galaxy, n_halo, n_bulge, n_disk):
    _dsk = Particles(len(galaxy))
    _dsk.mass = galaxy.mass
    _dsk.position = galaxy.position
    _dsk.velocity = galaxy.velocity
    disk = _dsk[n_bulge:n_halo]
    return disk


def simulate_merger(galaxy1, galaxy2, converter, n_halo, t_end, script_path, plot=False):
    converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    dynamics_code = Fi(converter, number_of_workers=1)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = dynamics_code.particles.add_particles(galaxy1)
    set2 = dynamics_code.particles.add_particles(galaxy2)
    
    dynamics_code.particles.move_to_center()
    dynamics_code.timestep = 0.5 | units.Myr
    
    disk1 = set1[:n_halo]
    disk2 = set2[:n_halo]
    
    if plot == True:
        make_plot(disk1, disk2, script_path, "Galaxy_merger_t0")
    
    current_iter = 0
    interval = 0.5 | units.Myr
    total_iter = int(t_end/interval) + 10
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    while dynamics_code.model_time < t_end:
        
        current_iter +=1
        
        if current_iter >= progress.maxval:
            progress.maxval = current_iter + 1
        
        dynamics_code.evolve_model(dynamics_code.model_time + interval)
                
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        make_plot(disk1, disk2, script_path,
                  "Galaxy_merger_t" + str(t_end.value_in(units.Myr))+"Myr")
        
    dynamics_code.stop()
    

def simulate_merger_with_particles(galaxy1, galaxy2, converter, n_halo, n_bulge, n_disk, t_end, script_path, plot=False):
    converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    dynamics_code = Fi(converter, number_of_workers=1)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    set1 = dynamics_code.particles.add_particles(galaxy1)
    set2 = dynamics_code.particles.add_particles(galaxy2)
    
    dynamics_code.particles.move_to_center()
    
    test_disk = test_particles(set2, n_halo, n_bulge, n_disk)
    
    test_code = Fi(converter, number_of_workers=1)
    test_code.parameters.epsilon_squared = (100 | units.parsec)**2
    set3 = test_code.particles.add_particles(test_disk)
    
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(test_code, (dynamics_code,) )
    gravity.timestep = 0.5 | units.Myr
    
    disk1 = set1[:n_halo]
    disk2 = set2[:n_halo]
    
    if plot == True:
        make_plot_testdisk(disk1, disk2, set3, script_path, "test_merger_t0")
    
    current_iter = 0
    interval = 0.5 | units.Myr
    total_iter = int(t_end/interval) + 10
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    while gravity.model_time < t_end:
        
        current_iter +=1
        
        gravity.evolve_model(gravity.model_time + interval)
                
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        make_plot_testdisk(disk1, disk2, set3, script_path,
                  "test_merger_t" + str(t_end.value_in(units.Myr))+"Myr")
        
    gravity.stop()
    

def mw_and_stars(galaxy1, stars, galaxy_converter, n_halo, t_end, script_path, plot=False):
    galaxy_converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    galaxy_dynamics_code = Fi(galaxy_converter, number_of_workers=4)
    galaxy_dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = galaxy_dynamics_code.particles.add_particles(galaxy1)
    
    galaxy_dynamics_code.particles.move_to_center()
    
    disk1 = set1[:n_halo]
    
    star_converter=nbody_system.nbody_to_si(stars.mass.sum(), 
                                   stars.position.length())

    star_dynamics_code = BHTree(star_converter)
    star_dynamics_code.particles.add_particles(stars)
    ch_g2l = star_dynamics_code.particles.new_channel_to(stars)
    #star_dynamics_code.timestep = 1|units.Myr
    
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(star_dynamics_code, (galaxy_dynamics_code,) )
    gravity.add_system(galaxy_dynamics_code, (star_dynamics_code,) )
    gravity.timestep = 0.5 | units.Myr
    
    if plot == True:
        make_plot_galstars(disk1, stars, script_path, "galstars_t0")
        
    x = [] | units.kpc
    y = [] | units.kpc
    
    current_iter = 0
    interval = 0.5
    total_iter = t_end/interval + 1
    
    times = np.arange(0., t_end, interval) | units.Myr
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    for time in times:
        
        current_iter +=1
        
        gravity.evolve_model(time)
        
        ch_g2l.copy()
        x.append(stars.x)
        y.append(stars.y)
                
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        make_plot_galstars(disk1, stars, script_path,
                  "galstars_t" + str(t_end.value_in(units.Myr))+"Myr")
        
    galaxy_dynamics_code.stop()
    star_dynamics_code.stop()


def merger_and_igm(galaxy1, galaxy2, converter, sph_code, n_halo, t_end, script_path, plot=False):
    converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    dynamics_code = Fi(converter, number_of_workers=4)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = dynamics_code.particles.add_particles(galaxy1)
    set2 = dynamics_code.particles.add_particles(galaxy2)
    
    dynamics_code.particles.move_to_center()
    
    disk1 = set1[:n_halo]
    disk2 = set2[:n_halo]
    
    if plot == True:
        make_plot(disk1, disk2, script_path, "sph_merger_t0")
    
    current_iter = 0
    interval = 0.5 | units.Myr
    total_iter = int(t_end/interval) + 1
    
    gravity_sph = bridge.Bridge(use_threading=False)
    gravity_sph.add_system(dynamics_code, (sph_code,) )
    gravity_sph.add_system(sph_code, (dynamics_code,) )
    gravity_sph.timestep = 0.5 | units.Myr
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    while dynamics_code.model_time < t_end:
        
        current_iter +=1
        
        gravity_sph.evolve_model(gravity_sph.model_time + interval)
                
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        make_plot(disk1, disk2, script_path,
                  "sph_merger_t" + str(t_end.value_in(units.Myr))+"Myr")
        
    gravity_sph.stop()