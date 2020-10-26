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


def plot_merger(mw_halo, mw_disk, m31_halo, m31_disk, title, script_path, filename):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-760, 760)
    plt.ylim(-760, 760)

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
    
    savepath = script_path + '/plots/merger_plots/'
    
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
    
    savepath = script_path + '/plots/zoomed_merger_plots/'
    
    plt.savefig(savepath + filename)
    
    
def contour_merger(galaxy1, galaxy2, title, script_path, filename):
    _x1, _x2 = np.array(galaxy1.x.value_in(units.kpc)), np.array(galaxy2.x.value_in(units.kpc))
    x = np.concatenate((_x1, _x2))
    y = np.concatenate((galaxy1.y.value_in(units.kpc), galaxy2.y.value_in(units.kpc)))
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    xy_range = [[-760, 760], [-760, 760]]
    
    counts,ybins,xbins,image = plt.hist2d(x, y, bins=400, range=xy_range,
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
    
    savepath = script_path + '/plots/merger_contour/'
    
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
    
    savepath = script_path + '/plots/solar_system_plots/'
    
    plt.savefig(savepath + filename)


def make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk, script_path, 
                  separation=780, 
                  M31_velocity=50, 
                  test=False):
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
    
    if not test:
        galaxy1_name = script_path + '/galaxies/data/M31_not_displaced_full'
        galaxy2_name = script_path + '/galaxies/data/MW_full'
    if test:
        galaxy1_name = script_path + '/galaxies/data/M31_not_displaced_test'
        galaxy2_name = script_path + '/galaxies/data/MW_test'
    
    widgets = ['Saving galaxies data: ', 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        write_set_to_file(galaxy1, galaxy1_name, 'hdf5')
        write_set_to_file(galaxy2, galaxy2_name, 'hdf5')

    return galaxy1, galaxy2, converter


def make_galaxy(converter, galaxy_dict, script_path, test=False):
    
    n_halo = galaxy_dict['n_halo']
    
    widgets = ['Building {} galaxy: '.format(galaxy_dict['name']), pbwg.AnimatedMarker(), ' ',
               pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy = new_galactics_model(n_halo,
                                     converter,
                                     halo_outer_radius = galaxy_dict['halo_outer_radius'],
                                     disk_number_of_particles = galaxy_dict['disk_number_of_particles'],
                                     disk_mass = galaxy_dict['disk_mass'],
                                     disk_scale_length = galaxy_dict['disk_scale_length'],
                                     disk_outer_radius = galaxy_dict['disk_outer_radius'],
                                     disk_scale_height_sech2 = galaxy_dict['disk_scale_height_sech2'],
                                     bulge_scale_radius = galaxy_dict['bulge_scale_radius'],
                                     bulge_number_of_particles = galaxy_dict['bulge_number_of_particles'])
                                     #central_radial_vel_dispersion = 0.7 * 100 | units.kms,
                                     #scale_length_of_sigR2 = 2.806 | units.kpc)
    
    if not test:
        galaxy_data_path = script_path + '/galaxies/data/{}_full'.format(galaxy_dict['name'])
    if test:
        galaxy_data_path = script_path + '/galaxies/data/{}_test'.format(galaxy_dict['name'])
    
    widgets = ['Saving {} galaxy data: '.format(galaxy_dict['name']), 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        write_set_to_file(galaxy, galaxy_data_path, 'hdf5')

    return galaxy


def displace_galaxy(galaxy, rotation_mat, translation_vector, radial_velocity, transverse_velocity):
    widgets = ['Adjusting relative velocities and orientations: ', 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        displaced_galaxy = Particles(len(galaxy))
        displaced_galaxy.mass = galaxy.mass
        displaced_galaxy.position = galaxy.position
        displaced_galaxy.velocity = galaxy.velocity
        
        for body in displaced_galaxy:
            body.position = ([body.x.value_in(units.kpc), 
                              body.y.value_in(units.kpc), 
                              body.z.value_in(units.kpc)] @ rotation_mat) | units.kpc
            body.velocity = ([body.vx.value_in(units.kms), 
                              body.vy.value_in(units.kms), 
                              body.vz.value_in(units.kms)] @ rotation_mat) | units.kms
        
        displaced_galaxy.position += translation_vector
        
        displaced_galaxy.velocity += radial_velocity
        displaced_galaxy.velocity += transverse_velocity
    
    return displaced_galaxy


def test_particles(galaxy, n_halo, n_bulge, n_disk):
    _dsk = Particles(len(galaxy))
    _dsk.mass = galaxy.mass
    _dsk.position = galaxy.position
    _dsk.velocity = galaxy.velocity
    disk = _dsk[n_bulge:n_halo]
    return disk


def simulate_merger(galaxy1, galaxy2, n_halo, t_end, script_path, plot=False):
    converter = nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
    
    dynamics_code = Gadget2(converter, number_of_workers=4)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = dynamics_code.dm_particles.add_particles(galaxy1)
    set2 = dynamics_code.dm_particles.add_particles(galaxy2)
    
    dynamics_code.particles.move_to_center()
    #dynamics_code.parameters.timestep = 0.5 | units.Myr
    
    halo1 = set1[n_halo:]
    disk1 = set1[:n_halo]
    
    halo2 = set2[n_halo:]
    disk2 = set2[:n_halo]
    
    if plot == True:
        plot_number = 0
        contour_merger(set1, set2, 
                "MW M31 merger\nt = 0 Myr", script_path, 'mw_m31_cmerger_' + str(plot_number).zfill(4))
        plot_merger(halo1, disk1, halo2, disk2, 
                    "MW M31 merger\nt = 0 Myr", script_path, 'mw_m31_merger_' + str(plot_number).zfill(4))
        plot_zoomed_merger(halo1, disk1, halo2, disk2, 
                           "MW M31 merger\nt = 0 Myr", script_path, 'mw_m31_merger_' + str(plot_number).zfill(4))
    
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
                contour_merger(set1, set2, 
                               "MW M31 merger\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                                                                               decimals=0))),
                               script_path, 'mw_m31_cmerger_' + str(plot_number).zfill(4))
                plot_merger(halo1, disk1, halo2, disk2, 
                            "MW M31 merger\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                                                                            decimals=0))),
                            script_path, 'mw_m31_merger_' + str(plot_number).zfill(4))
                plot_zoomed_merger(halo1, disk1, halo2, disk2, 
                                   "MW M31 merger\nt = {} Myr".format(int(np.round(dynamics_code.model_time.value_in(units.Myr), 
                                                                                   decimals=0))),
                                   script_path, 'mw_m31_merger_' + str(plot_number).zfill(4))
        
        progress.update(current_iter)
        
    progress.finish()
    
    if plot == True:
        plot_number += 1
        contour_merger(set1, set2,
                       "MW M31 merger\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),
                       script_path, 'mw_m31_cmerger_' + str(plot_number).zfill(4))
        plot_merger(halo1, disk1, halo2, disk2,
                    "MW M31 merger\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),
                    script_path, 'mw_m31_merger_' + str(plot_number).zfill(4))
        plot_zoomed_merger(halo1, disk1, halo2, disk2,
                           "MW M31 merger\nt = {} Myr".format(int(np.round(t_end.value_in(units.Myr), decimals=0))),
                           script_path, 'mw_m31_merger_' + str(plot_number).zfill(4))
        
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