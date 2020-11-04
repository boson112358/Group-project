from modules.common import __MERGER_DIR__, __PLOT_MERGER_DIR__, __CONTOUR_MERGER_DIR__, __ZOOM_MERGER_DIR__
from modules.progressbar import progressbar as pbar
from modules.progressbar import widgets as pbwg

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.colors import LogNorm

#from amuse.lab import units, Particles
#from amuse.community import Fi, Gadget2
from amuse.lab import *
from amuse.couple import bridge


###### create single galaxy plot output ######

def create_single_gal_dir(glxy_path):
    simul_dir = glxy_path + 'sim_plots/'
    if not os.path.exists(simul_dir):
        os.makedirs(simul_dir)
            
    return simul_dir


###### plot functions ######

def plot_merger(mw_halo, mw_disk, mw_bulge, m31_halo, m31_disk, m31_bulge, title, filename):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-500, 500)
    plt.ylim(-500, 500)

    #plotting mw_halo and mw_disk
    #ax.scatter(mw_halo.x.value_in(units.kpc), mw_halo.y.value_in(units.kpc),
               #c='tab:blue', alpha=0.6, s=1, lw=0)
    ax.scatter(mw_disk.x.value_in(units.kpc), mw_disk.y.value_in(units.kpc),
               c='cyan', alpha=1, s=1, lw=0, label='mw')
    ax.scatter(mw_bulge.x.value_in(units.kpc), mw_bulge.y.value_in(units.kpc),
               c='blue', alpha=1, s=1, lw=0, label='mw')
    
    #plotting m31_halo and m31_disk
    #ax.scatter(m31_halo.x.value_in(units.kpc), m31_halo.y.value_in(units.kpc),
               #c='tab:orange', alpha=0.6, s=1, lw=0)
    ax.scatter(m31_disk.x.value_in(units.kpc), m31_disk.y.value_in(units.kpc),
               c='orange', alpha=1, s=1, lw=0, label='m31')
    ax.scatter(m31_bulge.x.value_in(units.kpc), m31_bulge.y.value_in(units.kpc),
               c='red', alpha=1, s=1, lw=0, label='m31')
    
    plt.legend(loc='upper right')
    
    savepath = __SCRIPT_PATH__ + '/data/merger_plots/'
    
    plt.savefig(savepath + filename)
    
    
def plot_zoomed_merger(mw_halo, mw_disk, m31_halo, m31_disk, title, script_path, filename):
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

    #plotting mw_halo and mw_disk
    #ax.scatter(mw_halo.x.value_in(units.kpc), mw_halo.y.value_in(units.kpc),
               #c='tab:blue', alpha=0.6, s=1, lw=0)
    ax.scatter(mw_disk.x.value_in(units.kpc), mw_disk.y.value_in(units.kpc),
               c='tab:blue', alpha=1, s=1, lw=0, label='mw')
    
    #plotting m31_halo and m31_disk
    #ax.scatter(m31_halo.x.value_in(units.kpc), m31_halo.y.value_in(units.kpc),
               #c='tab:orange', alpha=0.6, s=1, lw=0)
    ax.scatter(m31_disk.x.value_in(units.kpc), m31_disk.y.value_in(units.kpc),
               c='tab:orange', alpha=1, s=1, lw=0, label='m31')
    
    plt.legend(loc='upper_right')
    
    savepath = script_path + '/plots/zoomed-merger-plots/'
    
    plt.savefig(savepath + filename)
    
    
def contour_merger(galaxy1, galaxy2, title, filename):
    _x1, _x2 = np.array(galaxy1.x.value_in(units.kpc)), np.array(galaxy2.x.value_in(units.kpc))
    x = np.concatenate((_x1, _x2))
    y = np.concatenate((galaxy1.y.value_in(units.kpc), galaxy2.y.value_in(units.kpc)))
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    xy_range = [[-760, 760], [-760, 760]]
    
    counts,ybins,xbins,image = plt.hist2d(x, y, bins=400, range=xy_range)
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    cont = ax.contourf(counts, norm=LogNorm(), levels=[1e-1, 1, 10, 100, 1e3, 1e4, 1e5],
                       extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()], cmap='RdGy')
    
    cbar = plt.colorbar(cont)
    cbar.set_label('Density')
    plt.tight_layout()
    
    savepath = __SCRIPT_PATH__ + '/data/merger_contour/'
    
    plt.savefig(savepath + filename)
    
    
def make_plot_testdisk(disk1, disk2, test_disk, title, script_path, filename):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    
    fig = plt.figure()
    
    plt.title(title)
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

    
def make_plot_galstars(disk, stars, title, script_path, filename):
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

    ax.scatter(disk.x.value_in(units.kpc), disk.y.value_in(units.kpc),
                   c='tab:blue', alpha=1, s=1, lw=0)
    ax.scatter(stars.x.value_in(units.kpc), stars.y.value_in(units.kpc),
                   c='tab:orange', alpha=1, marker='.', lw=1)
    
    savepath = script_path + '/plots/solar-system-plots/'
    
    plt.savefig(savepath + filename)
    
    
def plot_single_galaxy(halo, disk, bulge, title, glxy_path, filename):
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
    
    out_plot = create_single_gal_dir(glxy_path)
  
    plt.savefig(out_plot + filename)
    
    
###### plot condition function ######

def check_last_plot_time(current_time, last_plot_time, plot_interval, unit=units.Myr):
    if (current_time - last_plot_time).value_in(unit) >= plot_interval.value_in(unit):
        return True
    else:
        return False
    
    
###### merger function ######
    
def simulate_merger(galaxy1, galaxy2, n_halo, n_disk, n_bulge, t_end, script_path, interval=0.1|units.Myr, plot=False):
    converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    dynamics_code = Gadget2(converter, number_of_workers=4)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    if isinstance(dynamics_code, Gadget2) or isinstance(dynamics_code, Fi):
        #when using Gadget2 or Fi
        set1 = dynamics_code.dm_particles.add_particles(galaxy1)
        set2 = dynamics_code.dm_particles.add_particles(galaxy2)
    elif isinstance(dynamics_code, BHTree):
        #when using BHTree
        set1 = dynamics_code.particles.add_particles(galaxy1)
        set2 = dynamics_code.particles.add_particles(galaxy2)
    
    dynamics_code.particles.move_to_center()
    
    if isinstance(dynamics_code, Fi):
        dynamics_code.update_particle_set()
        
    dynamics_code.timestep = interval
    
    halo1 = set1[n_disk+n_bulge:]
    bulge1 = set1[n_disk:n_disk+n_bulge]
    disk1 = set1[:n_disk]
    
    halo2 = set2[n_disk+n_bulge:]
    bulge2 = set2[n_disk:n_disk+n_bulge]
    disk2 = set2[:n_disk]
    
    if plot == True:
        plot_number = 0
        last_plot_time = 0 | units.Myr
        contour_merger(set1, set2, 
                "MW M31 merger\nt = 0 Myr", 'mw_m31_cmerger_' + str(plot_number).zfill(4))
        plot_merger(halo1, disk1, bulge1, halo2, disk2, bulge2,
                    "MW M31 merger\nt = 0 Myr", 'mw_m31_merger_' + str(plot_number).zfill(4))
        #plot_zoomed_merger(halo1, disk1, halo2, disk2, 
        #                   "MW M31 merger\nt = 0 Myr", script_path, 'mw_m31_merger_' + str(plot_number).zfill(4))
    
    current_iter = 0
    t_end_in_Myr = t_end.as_quantity_in(units.Gyr)
    total_iter = int(t_end_in_Myr/interval) + 1
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    while dynamics_code.model_time < t_end:
        
        current_iter +=1
        
        dynamics_code.evolve_model(dynamics_code.model_time + interval)
        if isinstance(dynamics_code, Fi):
            dynamics_code.update_particle_set()
        
        if plot == True:
            if check_last_plot_time(dynamics_code.model_time, last_plot_time, t_end/100):
                plot_number += 1
                last_plot_time = dynamics_code.model_time
                contour_merger(set1, set2, 
                               "MW M31 merger\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                                                                               decimals=0))),
                               'mw_m31_cmerger_' + str(plot_number).zfill(4))
                plot_merger(halo1, disk1, bulge1, halo2, disk2, bulge2,
                            "MW M31 merger\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                                                                            decimals=0))),
                            'mw_m31_merger_' + str(plot_number).zfill(4))
                #plot_zoomed_merger(halo1, disk1, halo2, disk2, 
                #                   "MW M31 merger\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                #                                                                   decimals=0))),
                #                   script_path, 'mw_m31_merger_' + str(plot_number).zfill(4))
        
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        plot_number += 1
        contour_merger(set1, set2,
                       "MW M31 merger\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),
                       'mw_m31_cmerger_' + str(plot_number).zfill(4))
        plot_merger(halo1, disk1, bulge1, halo2, disk2, bulge2,
                    "MW M31 merger\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),
                    'mw_m31_merger_' + str(plot_number).zfill(4))
        #plot_zoomed_merger(halo1, disk1, halo2, disk2,
        #                   "MW M31 merger\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),
        #                   script_path, 'mw_m31_merger_' + str(plot_number).zfill(4))
        
    dynamics_code.stop()
    

def simulate_merger_with_particles(galaxy1, galaxy2, converter, n_halo, n_bulge, n_disk, t_end, script_path, plot=False):
    converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    dynamics_code = Fi(converter, redirection='none', number_of_workers=1)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    set1 = dynamics_code.dm_particles.add_particles(galaxy1)
    set2 = dynamics_code.dm_particles.add_particles(galaxy2)
    
    dynamics_code.particles.move_to_center()
    
    test_disk = test_particles(set2, n_halo, n_bulge, n_disk)
    
    test_code = Fi(converter, number_of_workers=1)
    test_code.parameters.epsilon_squared = (100 | units.parsec)**2
    set3 = test_code.dm_particles.add_particles(test_disk)
    
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(test_code, (dynamics_code,) )
    gravity.timestep = 0.5 | units.Myr
    
    disk1 = set1[:n_halo]
    disk2 = set2[:n_halo]
    
    if plot == True:
        make_plot_testdisk(disk1, disk2, set3, script_path, "test_merger_t0")
    
    current_iter = 0
    interval = 0.5 | units.Myr
    total_iter = int(t_end/interval) + 1
    
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
    
    
###### single galaxy simulation ######

def simulate_single_galaxy(galaxy1, converter, n_halo, n_bulge, n_disk, t_end, glxy_path,
                           solver=Gadget2, interval=0.5|units.Myr, plot=False, plot_freq=100):
    
    dynamics_code = solver(converter, number_of_workers=4)
    #dynamics_code = Fi(converter, redirection='none', number_of_workers=1)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = dynamics_code.dm_particles.add_particles(galaxy1)
    
    halo1 = set1[n_disk+n_bulge:]
    bulge1 = set1[n_disk:n_disk+n_bulge]
    disk1 = set1[:n_disk]
    
    dynamics_code.particles.move_to_center()
    
    if isinstance(dynamics_code, Gadget2):
        dynamics_code.timestep = interval
    
    if isinstance(dynamics_code, Fi):
        dynamics_code.parameters.timestep = interval
        dynamics_code.update_particle_set()
    
    if plot == True:
        plot_number = 0
        last_plot_time = 0 | units.Myr
        plot_single_galaxy(halo1, disk1, bulge1,
             "TEST\nt = 0 Myr", 
             glxy_path, 'mw_testrun_' + str(plot_number).zfill(4))
    
    current_iter = 0
    t_end_in_Myr = t_end.as_quantity_in(units.Gyr)
    total_iter = int(t_end_in_Myr/interval) + 10
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    while dynamics_code.model_time < t_end:
        
        current_iter +=1
        
        dynamics_code.evolve_model(dynamics_code.model_time + interval)
        
        if isinstance(dynamics_code, Fi):
            dynamics_code.update_particle_set()
        
        if plot == True:
            if check_last_plot_time(dynamics_code.model_time, last_plot_time, t_end/plot_freq):
                plot_number += 1
                last_plot_time = dynamics_code.model_time
                plot_single_galaxy(halo1, disk1, bulge1, 
                             "TEST\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                                                                  decimals=0))),
                             glxy_path, 'mw_testrun_' + str(plot_number).zfill(4))
        
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        plot_number += 1
        plot_single_galaxy(halo1, disk1, bulge1,
                     "TEST\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),                              
                     glxy_path, 'mw_testrun_' + str(plot_number).zfill(4))
        
    dynamics_code.stop()

def mw_and_stars(galaxy1, stars, galaxy_converter, star_solver, n_halo, t_end, script_path, plot=False):
    leapfrog = True
    if isinstance(star_solver, (BHTree, ph4)):
        leapfrog = False
    
    galaxy_converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    galaxy_dynamics_code = Fi(galaxy_converter, redirection='none', number_of_workers=1)
    galaxy_dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = galaxy_dynamics_code.particles.add_particles(galaxy1)
    
    galaxy_dynamics_code.particles.move_to_center()
    
    disk1 = set1[:n_halo]
    
    if leapfrog:
        galaxy_dynamics_code.parameters.timestep = 0.5 | units.Myr
        solver = galaxy_dynamics_code
    else:
        star_converter=nbody_system.nbody_to_si(stars.mass.sum(), 
                                                stars.position.length())
        star_dynamics_code = star_solver(star_converter)
        star_dynamics_code.particles.add_particles(stars)
        ch_g2l = star_dynamics_code.particles.new_channel_to(stars)
    
        gravity = bridge.Bridge(use_threading=False)
        gravity.add_system(star_dynamics_code, (galaxy_dynamics_code,) )
        gravity.timestep = 0.5 | units.Myr
        solver = gravity
    
    if plot == True:
        plot_number = 0
        make_plot_galstars(disk1, stars, 
                           "MW and Solar System\nt = 0 Myr", 
                           script_path, 'galstars_' + str(plot_number).zfill(4))
        
    x = [] | units.kpc
    y = [] | units.kpc
    
    current_iter = 0
    interval = 0.5
    total_iter = int(t_end/(interval | units.Myr)) + 1
    t_end_scalar = t_end.value_in(units.Myr)
    #print(t_end_scalar)
    
    times = np.arange(0., t_end_scalar, interval) | units.Myr
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    for time in times:
        
        solver.evolve_model(time)
        
        if leapfrog:
            star_solver(current_iter, stars, galaxy_dynamics_code)
        else:
            ch_g2l.copy()
        
        x.append(stars.x)
        y.append(stars.y)
        
        if plot == True:
            if current_iter in [10*i for i in range(1, total_iter)]:
                plot_number += 1
                make_plot_galstars(disk1, stars, 
                                   "MW and Solar System\nt = {} Myr".format(np.round(solver.model_time.value_in(units.Myr),
                                                                                     decimals=0)), 
                                   script_path, 'galstars_' + str(plot_number).zfill(4))
        
        current_iter +=1
        
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        plot_number += 1
        make_plot_galstars(disk1, stars, 
                           "MW and Solar System\nt = {} Myr".format(t_end.value_in(units.Myr)), 
                           script_path, 'galstars_' + str(plot_number).zfill(4))
    
    galaxy_dynamics_code.stop()
    if not leapfrog:
        star_dynamics_code.stop()