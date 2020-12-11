###### importing modules ######

from modules.common import __MERGER_DIR__, __SOLAR_DIR__, galaxy_structures
from modules.common import __SCRIPT_PATH__, __ANIMATION_DIR__, __FRAME_DIR__
from modules.progressbar import progressbar as pbar
from modules.progressbar import widgets as pbwg

import sys
import os
import datetime
import imageio
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from amuse.lab import units, Particles, Fi, Gadget2, BHTree, ParticlesSuperset
from amuse.couple import bridge


###### create merger output dir ######

def create_merger_output_dir():
    parent = __MERGER_DIR__
    
    increasing = 1
    
    while True:
        current_merger = 'merger_' + str(datetime.date.today()) + '-' + str(increasing).zfill(4)
        out_dir = parent + current_merger + '/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)
            break
        else:
            increasing += 1
            
    return out_dir, current_merger

def create_merger_subdirs(current_out_dir):
    diskbulge_merger_dir = current_out_dir + '/merger_plots/'
    contour_merger_dir = current_out_dir + '/merger_contour/'
    zoom_merger_dir = current_out_dir + '/merger_zoom/'
    mw_zoom_dir = current_out_dir + '/merger_mwzoom/'
    
    histogram_dir = current_out_dir + '/merger_solar_histogram/'
    igm_dir = current_out_dir + '/merger_igm/'
    
    dirs = [diskbulge_merger_dir, contour_merger_dir, zoom_merger_dir, 
            mw_zoom_dir, histogram_dir, igm_dir]
    
    for direct in dirs:
        if not os.path.exists(direct):
            os.makedirs(direct, exist_ok=True)
            
    return diskbulge_merger_dir, contour_merger_dir, zoom_merger_dir, mw_zoom_dir, histogram_dir, igm_dir


###### create single galaxy plot output ######

def create_single_gal_dir(glxy_path):
    simul_dir = glxy_path + 'sim_plots/'
    if not os.path.exists(simul_dir):
        os.makedirs(simul_dir, exist_ok=True)
            
    return simul_dir


###### amination funnctions ######

class GalaxyAnimWriter(object):
    
    def __new__(self, out_name, out_format='.mp4', fps=20):
        return imageio.get_writer(__ANIMATION_DIR__ + out_name + out_format, fps=fps)
    
    def __init__(self, out_name, out_format='.mp4', fps=20):
        self.out_name = out_name
        self.out_format = out_format
        self.fps = fps
        
def read_frame(filename, parentdir=__FRAME_DIR__):
    return imageio.imread(parentdir + filename)


###### plot functions ######

def plot_diskbulge_merger(mw_halo, mw_disk, mw_bulge, m31_halo, m31_disk, m31_bulge, title, savepath, filename,
                          is_snapshot=True, is_frame=True):
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
    
    if is_snapshot:
        plt.savefig(savepath + filename)
    if is_frame:
        framename = '_diskbulge_temp_frame'
        plt.savefig(__FRAME_DIR__ + framename)
    
    
def plot_zoomed_merger(mw_disk, mw_bulge, m31_disk, m31_bulge, title, savepath, filename, particles=None,
                       is_snapshot=True, is_frame=True):
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
    ax.scatter(mw_bulge.x.value_in(units.kpc), mw_bulge.y.value_in(units.kpc),
               c='tab:blue', alpha=0.6, s=1, lw=0)
    ax.scatter(mw_disk.x.value_in(units.kpc), mw_disk.y.value_in(units.kpc),
               c='tab:blue', alpha=1, s=1, lw=0, label='mw')
    
    #plotting m31_halo and m31_disk
    ax.scatter(m31_bulge.x.value_in(units.kpc), m31_bulge.y.value_in(units.kpc),
               c='tab:orange', alpha=0.6, s=1, lw=0)
    ax.scatter(m31_disk.x.value_in(units.kpc), m31_disk.y.value_in(units.kpc),
               c='tab:orange', alpha=1, s=1, lw=0, label='m31')
    
    #plotting stars
    if particles != None:
        ax.scatter(particles.x.value_in(units.kpc), particles.y.value_in(units.kpc),
                   c='tab:pink', alpha=1, marker='.', lw=1, label='solar system')
    
    plt.legend(loc='upper right')
    
    if is_snapshot:
        plt.savefig(savepath + filename)
    if is_frame:
        framename = '_zoom_temp_frame'
        plt.savefig(__FRAME_DIR__ + framename)
    
    
def plot_contour_merger(galaxy1, galaxy2, title, savepath, filename,
                        is_snapshot=True, is_frame=True):
    _x1, _x2 = np.array(galaxy1.x.value_in(units.kpc)), np.array(galaxy2.x.value_in(units.kpc))
    x = np.concatenate((_x1, _x2))
    y = np.concatenate((galaxy1.y.value_in(units.kpc), galaxy2.y.value_in(units.kpc)))
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    xy_range = [[-500, 500], [-500, 500]]
    
    counts,xbins,ybins,image = plt.hist2d(y, x, bins=400, range=xy_range)
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    cont = ax.contourf(counts, norm=LogNorm(), levels=[1e-1, 1, 10, 100, 1e3, 1e4, 1e5],
                       extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()], cmap='cool')
    
    cbar = plt.colorbar(cont)
    cbar.set_label('Density')
    plt.tight_layout()
    
    if is_snapshot:
        plt.savefig(savepath + filename)
    if is_frame:
        framename = '_contour_temp_frame'
        plt.savefig(__FRAME_DIR__ + framename)
        
        
def plot_mw_zoom(mw_galaxy, mw_disk, mw_bulge, m31_disk, m31_bulge, title, savepath, filename,
                 particles=None, is_snapshot=True, is_frame=True):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    com = mw_galaxy.center_of_mass()
    xcom = com.x.value_in(units.kpc)
    ycom = com.y.value_in(units.kpc)
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-60 + xcom, xcom + 60)
    plt.ylim(-60 + ycom, ycom + 60)

    #plotting mw_halo and mw_disk
    ax.scatter(mw_disk.x.value_in(units.kpc), mw_disk.y.value_in(units.kpc),
               c='tab:blue', alpha=1, s=1, lw=0)
    ax.scatter(mw_bulge.x.value_in(units.kpc), mw_bulge.y.value_in(units.kpc),
               c='tab:blue', alpha=1, s=1, lw=0, label='mw')
    
    #plotting m31_halo and m31_disk
    ax.scatter(m31_disk.x.value_in(units.kpc), m31_disk.y.value_in(units.kpc),
               c='tab:orange', alpha=1, s=1, lw=0)
    ax.scatter(m31_bulge.x.value_in(units.kpc), m31_bulge.y.value_in(units.kpc),
               c='tab:orange', alpha=1, s=1, lw=0, label='m31')
    
    #plotting stars
    if particles != None:
        ax.scatter(particles.x.value_in(units.kpc), particles.y.value_in(units.kpc),
                   c='tab:pink', alpha=1, marker='.', lw=1, label='solar system')
    
    plt.legend(loc='upper right')
    
    if is_snapshot:
        plt.savefig(savepath + filename)
    if is_frame:
        framename = '_mwzoom_temp_frame'
        plt.savefig(__FRAME_DIR__ + framename)
        
        
def solar_plot(particles, is_snapshot=True):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    plt.title('sol_test')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    ax.scatter(particles.x.value_in(units.kpc), particles.y.value_in(units.kpc),
                   c='tab:pink', alpha=1, marker='.', lw=5, label='solar system')
    plt.legend(loc='upper right')
    if is_snapshot:
        plt.savefig(__MERGER_DIR__ + 'sol_test')
    
def make_plot_galstars(disk, bulge, stars, title, filename):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-40, 40)
    plt.ylim(-40, 40)

    ax.scatter(disk.x.value_in(units.kpc), disk.y.value_in(units.kpc),
                   c='tab:blue', alpha=1, s=1, lw=0)
    ax.scatter(bulge.x.value_in(units.kpc), bulge.y.value_in(units.kpc),
                   c='tab:blue', alpha=1, s=1, lw=0)
    ax.scatter(stars.x.value_in(units.kpc), stars.y.value_in(units.kpc),
                   c='tab:orange', alpha=1, marker='.', lw=1)
    
    savepath = __SOLAR_DIR__
    
    plt.savefig(savepath + filename)
    
    plt.close()
    
    
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
    
    
###### plot wrapper ######

def plt_anim_wrapper(galaxy1, galaxy2, n_disk, n_bulge, 
                     last_snap_number, current_time, last_plot_time, common_title, 
                     diskbulge_merger_dir, contour_merger_dir, zoom_merger_dir, mw_zoom_dir,
                     filename_prefix, t_end, snap_freq, particles=None,
                     snapshot=True, animation=True, 
                     diskbulge=True, contour=True, zoom=True, mwzoom=True,
                     diskbulge_writer=None, contour_writer=None, zoom_writer=None):
    
    is_snapshot = False
    
    if snapshot:
        if check_last_plot_time(current_time, last_plot_time, t_end/snap_freq):
            is_snapshot = True
            current_snap_number = last_snap_number + 1
            
    if not is_snapshot:
        current_snap_number = -1
    
    halo1, disk1, bulge1 = galaxy_structures(galaxy1, n_disk, n_bulge)
    halo2, disk2, bulge2 = galaxy_structures(galaxy2, n_disk, n_bulge)
    
    full_title = common_title + '\nt = {} Myr'.format(int(np.round(current_time.value_in(units.Myr), decimals=0))) 
    
    if diskbulge:
        try:
            plot_diskbulge_merger(halo1, disk1, bulge1, halo2, disk2, bulge2,
                                  full_title,
                                  diskbulge_merger_dir,
                                  filename_prefix + '_diskbulge_merger_' + str(current_snap_number).zfill(4),
                                  is_snapshot=is_snapshot, is_frame=animation)
        except:
            pass
        
        if animation:
            img = read_frame('_diskbulge_temp_frame.png')
            diskbulge_writer.append_data(img)
    
    if contour:
        try:
            plot_contour_merger(galaxy1, galaxy2, 
                                full_title,
                                contour_merger_dir,
                                filename_prefix + '_contour_merger_' + str(current_snap_number).zfill(4),
                                is_snapshot=is_snapshot, is_frame=animation)
        except:
            pass
        
        if animation:
            img = read_frame('_contour_temp_frame.png')
            contour_writer.append_data(img)
        
    if zoom:
        try:
            plot_zoomed_merger(disk1, bulge1, disk2, bulge2, 
                               full_title,
                               zoom_merger_dir, 
                               filename_prefix + '_zoomed_merger_' + str(current_snap_number).zfill(4),
                               particles=particles,
                               is_snapshot=is_snapshot, is_frame=animation)
        except:
            pass
        
        if animation:
            img = read_frame('_zoom_temp_frame.png')
            zoom_writer.append_data(img)
            
            
     
    if mwzoom:
        try:
            plot_mw_zoom(galaxy1, disk1, bulge1, disk2, bulge2, 
                         full_title, 
                         mw_zoom_dir, 
                         filename_prefix + '_mwzoomed_merger_' + str(current_snap_number).zfill(4),
                         particles=particles, 
                         is_snapshot=is_snapshot, is_frame=animation)
        except:
            pass
        
    #solar_plot(particles)
            
    plt.close('all')
    
    if is_snapshot:
        return current_time, current_snap_number
    else:
        return last_plot_time, last_snap_number
    
    
###### igm plot function ######

def plot_igm(sph_code, Ngrid, Lg, output_dir, filename):
    x, y, z, vx, vy, vz = setup_grid(Ngrid, Lg)
    rho, rhovx, rhovy, rhovz, rhoe = sph_code.get_hydro_state_at_point(x, y, z, vx, vy, vz)
    rho = rhoe.reshape((Ngrid + 1, Ngrid + 1))
    max_dens = np.log10(rho.value_in(units.m**-1 * units.s**-2 * units.kg)).max()
    min_dens = np.log10(rho.value_in(units.m**-1 * units.s**-2 * units.kg)).min()
    #print("Done running")
    # plt.scatter(sph_code.gas_particles.x.value_in(units.kpc),
    #             sph_code.gas_particles.y.value_in(units.kpc),
    #             c = 'r')
    # plt.scatter(sph_code.dm_particles.x.value_in(units.kpc),
    #             sph_code.dm_particles.y.value_in(units.kpc),
    #             c = 'b')
    cax = plt.imshow(np.log10(rho.value_in(units.m**-1 * units.s**-2 * units.kg)),
                     extent=[-Lg, Lg, -Lg, Lg], vmin=min_dens, vmax=max_dens,
                     origin='lower', cmap="hot")
    plt.colorbar(cax, ticks=[1.e-8, 0.5*max_dens,max_dens], orientation='vertical', fraction=0.045)
    plt.savefig(output_dir + filename)
    plt.close()

def setup_grid(N, L):
    x,y=np.indices( ( N+1,N+1 ))
    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=0.*x 
    vx=0.*x
    vy=0.*x
    vz=0.*x
    x=units.kpc(x)
    y=units.kpc(y)
    z=units.kpc(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)
    return x, y, z, vx, vy, vz

    
###### plot condition function ######

def check_last_plot_time(current_time, last_plot_time, plot_interval, unit=units.Myr):
    if current_time.value_in(unit) == 0:
        return True
    if (current_time - last_plot_time).value_in(unit) >= plot_interval.value_in(unit):
        return True
    else:
        return False
    
    
###### animation closing function ######

def close_animation(current_time, last_anim_time, anim_interval=(500 | units.Myr), unit=units.Myr):
     if (current_time - last_anim_time).value_in(unit) >= anim_interval.value_in(unit):
        return True
     else:
        return False
    
    
###### galaxy separation ######

def separation(glxy1, glxy2):
    cm1 = glxy1.center_of_mass()
    cm2 = glxy2.center_of_mass()
    sep = np.sqrt((cm2.x - cm1.x)**2 + (cm2.y - cm1.y)**2 + (cm2.z + cm1.z)**2)
    
    return sep

def sep_logs(sep_list, sep, current_time):
    sep_list[0].append(current_time)
    sep_list[1].append(sep)
    
    return sep_list


###### solar system position ######

def solar_logs(solar_pos_list, mw_com, particles, current_time):
    solar_pos_list[0].append(current_time)
    solar_pos_list[1].append(mw_com.x.value_in(units.kpc))
    solar_pos_list[2].append(mw_com.y.value_in(units.kpc))
    solar_pos_list[3].append(mw_com.z.value_in(units.kpc))
    for i in range(0, len(particles)):
        solar_pos_list[4+3*i].append(particles[i-1].x.value_in(units.kpc))
        solar_pos_list[4+3*i+1].append(particles[i-1].y.value_in(units.kpc))
        solar_pos_list[4+3*i+2].append(particles[i-1].z.value_in(units.kpc))
        
    return solar_pos_list


###### update merger dictionary ######

def update_merger_dict(mrgr_dict, current_time, glxy1, glxy2, solar=None):
    
    glxy1_com = glxy1.center_of_mass()
    glxy2_com = glxy2.center_of_mass()
    one_and_two = ParticlesSuperset([glxy1, glxy2])
    one_and_two_com = one_and_two.center_of_mass()
    
    key_list = ['time (Myr)', 'mw_com_x (kpc)', 'mw_com_y (kpc)', 'mw_com_z (kpc)', 
                'm31_com_x (kpc)', 'm31_com_y (kpc)', 'm31_com_z (kpc)', 
                'mw_m31_com_x (kpc)', 'mw_m31_com_y (kpc)', 'mw_m31_com_z (kpc)']
    
    value_list = [current_time.value_in(units.Myr), 
                  glxy1_com.x.value_in(units.kpc), 
                  glxy1_com.y.value_in(units.kpc), 
                  glxy1_com.z.value_in(units.kpc),
                  glxy2_com.x.value_in(units.kpc), 
                  glxy2_com.y.value_in(units.kpc), 
                  glxy2_com.x.value_in(units.kpc),
                  one_and_two_com.x.value_in(units.kpc), 
                  one_and_two_com.y.value_in(units.kpc), 
                  one_and_two_com.z.value_in(units.kpc)]
    
    #check if the key exists, if not creates that key with a list
    for key, value in zip(key_list, value_list):
        if key in mrgr_dict:
            mrgr_dict[key].append(value)
        else:
            mrgr_dict[key] = [value]
    
    if not solar == None:
        solar_com = solar.center_of_mass()
        
        solar_com_keys = ['solar_com_x (kpc)', 'solar_com_y (kpc)', 'solar_com_z (kpc)']                
        solar_com_values = [solar_com.x.value_in(units.kpc), 
                            solar_com.y.value_in(units.kpc), 
                            solar_com.z.value_in(units.kpc)]
        
        for key, value in zip(solar_com_keys, solar_com_values):
            if key in mrgr_dict:
                mrgr_dict[key].append(value)
            else:
                mrgr_dict[key] = [value]
        
        solar_position_keys = [['x_tr_' + str(index).zfill(5) + ' (kpc)', 
                                'y_tr_' + str(index).zfill(5) + ' (kpc)', 
                                'z_tr_' + str(index).zfill(5) + ' (kpc)'] for index in range(len(solar))]
        solar_position_values = [[tracker.x.value_in(units.kpc), 
                                  tracker.y.value_in(units.kpc), 
                                  tracker.z.value_in(units.kpc)] for tracker in solar]
        
        for tracker_key_list, tracker_value_list in zip(solar_position_keys, solar_position_values):
            for key, value in zip(tracker_key_list, tracker_value_list):
                if key in mrgr_dict:
                    mrgr_dict[key].append(value)
                else:
                    mrgr_dict[key] = [value]
                         
    return mrgr_dict

def pad_dict_list(dict_list, padel):
    lmax = 0
    for lname in dict_list.keys():
        lmax = max(lmax, len(dict_list[lname]))
    for lname in dict_list.keys():
        ll = len(dict_list[lname])
        if  ll < lmax:
            dict_list[lname] += [padel] * (lmax - ll)
    return dict_list

    
###### merger function ######

def simulate_merger(galaxy1, galaxy2, n_halo, n_disk, n_bulge, t_end, converter, solver, initial_conditions_dict,
                    interval=0.5|units.Myr, animation=False, snapshot=False, snap_freq=1000,
                    sol_system=None, igm_gas_particles=None, igm_dm_particles=None, box_grid=0):

    #sets up the gravity solver
    dynamics_code = solver(converter, number_of_workers=4)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    dynamics_code.timestep = interval
    #if igm_gas_particles != None and igm_dm_particles != None:
        #dynamics_code.parameters.self_gravity_flag = True
    
    #computes coordinate difference between mw center of mass and solar system and moves it
    if sol_system != None:
        deltapos = galaxy1.center_of_mass() - sol_system.position
        sol_system.position = galaxy1.center_of_mass() - deltapos
    
    #uploads the particles into the gravity solver
    set1 = dynamics_code.dm_particles.add_particles(galaxy1)
    set2 = dynamics_code.dm_particles.add_particles(galaxy2)
    if sol_system != None:
        set5 = dynamics_code.dm_particles.add_particles(sol_system)
    else:
        set5 = None

    #moves system to center of mass 
    dynamics_code.particles.move_to_center()
    
    #uploads igm particles
    if igm_gas_particles != None and igm_dm_particles != None:
        set3 = dynamics_code.dm_particles.add_particles(igm_gas_particles)
        set4 = dynamics_code.dm_particles.add_particles(igm_dm_particles)
    
    #creates channels
    mw_channel = dynamics_code.dm_particles.new_channel_to(set1)
    m31_channel = dynamics_code.dm_particles.new_channel_to(set2)
    if igm_gas_particles != None and igm_dm_particles != None:
        igm_gas_channel = dynamics_code.dm_particles.new_channel_to(set3)
        igm_dm_channel = dynamics_code.dm_particles.new_channel_to(set4)
    if sol_system != None:
        solar_channel = dynamics_code.dm_particles.new_channel_to(set5)
    
    #checks if the solver is Fi
    if isinstance(dynamics_code, Fi):
        dynamics_code.parameters.timestep = interval
        dynamics_code.update_particle_set()
        
    #creates output dirs
    out_dir, current_merger = create_merger_output_dir()
    diskbulge_merger_dir, contour_merger_dir, zoom_merger_dir, mw_zoom_dir, histogram_dir, igm_dir = create_merger_subdirs(out_dir)

    #initializes animations
    if animation:
        animation_number = 1
        last_anim_time = 0 | units.Myr
        db_writer = GalaxyAnimWriter(current_merger + '_diskbulge_anim' + str(animation_number).zfill(2))
        cr_writer = GalaxyAnimWriter(current_merger + '_contour_anim'+ str(animation_number).zfill(2))
        zm_writer = GalaxyAnimWriter(current_merger + '_zoom_anim'+ str(animation_number).zfill(2))
    else:
        db_writer = None
        cr_writer = None
        zm_writer = None

    #initializes snapshots
    start_zoom_plot = False
    last_snap_number = 0
    last_plot_time = 0 | units.Myr
    current_time = 0 | units.Myr
    _a, _b = plt_anim_wrapper(set1, set2, n_disk, n_bulge,
                              last_snap_number, current_time, last_plot_time, 'MW M31 merger',
                              diskbulge_merger_dir, contour_merger_dir, zoom_merger_dir, mw_zoom_dir,
                              'mw_m31', t_end, 100000, particles=set5,
                              animation=animation, snapshot=snapshot,
                              diskbulge_writer=db_writer, contour_writer=cr_writer, zoom_writer=zm_writer,
                              zoom=start_zoom_plot)
    
    #plots starting IGM
    #if igm_gas_particles != None and igm_dm_particles != None:
    #    plot_igm(set3|set4, 500, Lg, igm_dir, 'IGM_density_t0')

    #initializes separation logs and solar system positions
    merger_dict = initial_conditions_dict       
    merger_dict = update_merger_dict(merger_dict, current_time, set1, set2, solar=set5)

    #initializes progressbar
    current_iter = 0
    t_end_in_Myr = t_end.as_quantity_in(units.Myr)
    total_iter = int(t_end_in_Myr/interval) + 1
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '),
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    while dynamics_code.model_time < t_end:
        
        #evolves model
        dynamics_code.evolve_model(dynamics_code.model_time + interval)
        if isinstance(dynamics_code, Fi):
            dynamics_code.update_particle_set()
        
        #copies channels
        mw_channel.copy()
        m31_channel.copy()
        if sol_system != None:
            solar_channel.copy()
        if igm_gas_particles != None and igm_dm_particles != None:
            igm_gas_channel.copy()
            igm_dm_channel.copy()
        
        #current time of the solver
        current_time = dynamics_code.model_time
        
        #makes plots
        if last_plot_time.value_in(units.Myr) > 3000:
            start_zoom_plot = True
        last_plot_time, last_snap_number = plt_anim_wrapper(set1, set2, n_disk, n_bulge,
                                                            last_snap_number, current_time, last_plot_time,
                                                            'MW M31 merger',
                                                            diskbulge_merger_dir, contour_merger_dir,
                                                            zoom_merger_dir, mw_zoom_dir,
                                                            'mw_m31', t_end, snap_freq, particles=set5,
                                                            animation=animation, snapshot=snapshot,
                                                            diskbulge_writer=db_writer, contour_writer=cr_writer,
                                                            zoom_writer=zm_writer, zoom=start_zoom_plot)
        
        #updates separation logs and solar system positions
        merger_dict = update_merger_dict(merger_dict, current_time, set1, set2, solar=set5)
        
        #manages animations
        if animation:
            if close_animation(dynamics_code.model_time, last_anim_time):
                db_writer.close()
                cr_writer.close()
                zm_writer.close()
                animation_number += 1
                last_anim_time = dynamics_code.model_time
                db_writer = GalaxyAnimWriter(current_merger + '_diskbulge_anim' + str(animation_number).zfill(2))
                cr_writer = GalaxyAnimWriter(current_merger + '_contour_anim'+ str(animation_number).zfill(2))
                zm_writer = GalaxyAnimWriter(current_merger + '_zoom_anim'+ str(animation_number).zfill(2))
        
        #updates progressbar
        current_iter +=1
        progress.update(current_iter)
    
    #closes progressbar
    progress.finish()
    
    #plots final IGM
    #if igm_gas_particles != None and igm_dm_particles != None:
    #    plot_igm(igm_gas_particles|igm_dm_particles, 500, Lg, igm_dir, 'IGM_density_tfinal')
    
    #kills solver
    dynamics_code.stop()
    
    #checks if every dict key list has the same numer of elements
    merger_dict = pad_dict_list(merger_dict, -1)
    
    #saves merger logs into dataframe
    df_merger = pd.DataFrame(merger_dict)
    df_merger.to_csv('{}merger_logs.csv'.format(out_dir), index=False)
                           
    #closes animations
    if animation:
        db_writer.close()
        cr_writer.close()
        zm_writer.close()
    

"""
def simulate_merger(galaxy1, galaxy2, n_halo, n_disk, n_bulge, t_end, converter, 
                    solver=Gadget2, interval=0.5|units.Myr, 
                    animation=False, snapshot=False, snap_freq=5000,
                    particles=None, particles_solver=None):
    
    #sets up the gravity solver
    if particles == None:
        dynamics_code = Gadget2(converter, number_of_workers=4)
    else:
        dynamics_code = Fi(converter, redirection='none', number_of_workers=1)
    
    #dynamics_code = solver(converter, number_of_workers=4)
    dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    if isinstance(dynamics_code, Gadget2) or isinstance(dynamics_code, Fi):
        #when using Gadget2 or Fi
        set1 = dynamics_code.dm_particles.add_particles(galaxy1)
        set2 = dynamics_code.dm_particles.add_particles(galaxy2)
    elif isinstance(dynamics_code, BHTree):
        #when using BHTree
        set1 = dynamics_code.particles.add_particles(galaxy1)
        set2 = dynamics_code.particles.add_particles(galaxy2)
    
    #computes coordinate difference between mw center of mass and solar system
    if particles != None:
        deltapos = set1.center_of_mass() - particles.position
    
    #moves system to center of mass and creates channels
    dynamics_code.particles.move_to_center()
    mw_channel = dynamics_code.particles.new_channel_to(set1)
    m31_channel = dynamics_code.particles.new_channel_to(set2)
    
    #moves the solar system particles
    if particles != None:
        particles.position = set1.center_of_mass() - deltapos 
        
    if isinstance(dynamics_code, Gadget2):
        dynamics_code.timestep = interval
    elif isinstance(dynamics_code, Fi):
        dynamics_code.parameters.timestep = interval
        dynamics_code.update_particle_set()
    
    #creates output dirs
    out_dir, current_merger = create_merger_output_dir()
    diskbulge_merger_dir, contour_merger_dir, zoom_merger_dir, mw_zoom_dir = create_merger_subdirs(out_dir)
    
    if animation:
        #initializes animations
        animation_number = 1
        last_anim_time = 0 | units.Myr
        db_writer = GalaxyAnimWriter(current_merger + '_diskbulge_anim' + str(animation_number).zfill(2))
        cr_writer = GalaxyAnimWriter(current_merger + '_contour_anim'+ str(animation_number).zfill(2))
        zm_writer = GalaxyAnimWriter(current_merger + '_zoom_anim'+ str(animation_number).zfill(2))
    else:
        db_writer = None
        cr_writer = None
        zm_writer = None
    
    #initializes snapshots
    start_zoom_plot = False
    last_snap_number = 0
    last_plot_time = 0 | units.Myr
    current_time = 0 | units.Myr
    _a, _b = plt_anim_wrapper(set1, set2, n_disk, n_bulge, 
                              last_snap_number, current_time, last_plot_time, 'MW M31 merger', 
                              diskbulge_merger_dir, contour_merger_dir, zoom_merger_dir, mw_zoom_dir,
                              'mw_m31', t_end, 100000, particles=particles,
                              animation=animation, snapshot=snapshot,
                              diskbulge_writer=db_writer, contour_writer=cr_writer, zoom_writer=zm_writer,
                              zoom=start_zoom_plot)
    
    #initializes separation
    sep_list = [[], []]
    sep_list = sep_logs(sep_list, separation(set1, set2).value_in(units.kpc), current_time.value_in(units.Myr))
    
    current_iter = 0
    t_end_in_Myr = t_end.as_quantity_in(units.Myr)
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
            
        if particles != None:
            particles_solver(current_iter, particles, dynamics_code)
            
        mw_channel.copy()
        m31_channel.copy()
        
        if last_plot_time.value_in(units.Myr) > 3000:
            start_zoom_plot = True
        last_plot_time, last_snap_number = plt_anim_wrapper(set1, set2, n_disk, n_bulge, 
                                                            last_snap_number, dynamics_code.model_time, last_plot_time, 
                                                            'MW M31 merger', 
                                                            diskbulge_merger_dir, contour_merger_dir, 
                                                            zoom_merger_dir, mw_zoom_dir,
                                                            'mw_m31', t_end, snap_freq, particles=particles,
                                                            animation=animation, snapshot=snapshot,
                                                            diskbulge_writer=db_writer, contour_writer=cr_writer,
                                                            zoom_writer=zm_writer, zoom=start_zoom_plot)
        sep_list = sep_logs(sep_list, separation(set1, set2).value_in(units.kpc), 
                            dynamics_code.model_time.value_in(units.Myr))
        
        if animation:
            if close_animation(dynamics_code.model_time, last_anim_time):
                db_writer.close()
                cr_writer.close()
                zm_writer.close()
                animation_number += 1
                last_anim_time = dynamics_code.model_time
                db_writer = GalaxyAnimWriter(current_merger + '_diskbulge_anim' + str(animation_number).zfill(2))
                cr_writer = GalaxyAnimWriter(current_merger + '_contour_anim'+ str(animation_number).zfill(2))
                zm_writer = GalaxyAnimWriter(current_merger + '_zoom_anim'+ str(animation_number).zfill(2))
        
        progress.update(current_iter)
        
    progress.finish()       
    dynamics_code.stop()
    
    #creates separation dictionary
    sep_dict = {}
    sep_dict.update({'time (Myr)': sep_list[0]})
    sep_dict.update({'sep (kpc)': sep_list[1]})
    df_sep = pd.DataFrame(sep_dict)
    df_sep.to_csv('{}separation.csv'.format(out_dir), index=False)
    
    if animation:
        db_writer.close()
        cr_writer.close()
        zm_writer.close()    
"""


###### single galaxy simulation ######

def simulate_single_galaxy(galaxy1, converter, n_halo, n_bulge, n_disk, t_end, glxy_path,
                           solver=Gadget2, interval=0.5|units.Myr, plot=False, plot_freq=100):
    
    if isinstance(solver, Gadget2):
        nw = 4
    elif isinstance(solver, Fi):
        nw = 1
        
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
    elif isinstance(dynamics_code, Fi):
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

def mw_and_stars(galaxy1, stars, converter, n_disk, n_bulge, t_end, star_solver,
                 interval=1|units.Myr, snapshot=False, snap_freq=300):
    
    leapfrog = True
    if isinstance(star_solver, (BHTree)):
        leapfrog = False
    
    galaxy_dynamics_code = Fi(converter, redirection='none', number_of_workers=1)
    galaxy_dynamics_code.parameters.epsilon_squared = (100 | units.parsec)**2
    
    set1 = galaxy_dynamics_code.dm_particles.add_particles(galaxy1)
    mw_channel = galaxy_dynamics_code.particles.new_channel_to(set1)
    galaxy_dynamics_code.particles.move_to_center()
    #mw_channel.update()
    
    halo1, disk1, bulge1 = galaxy_structures(set1, n_disk, n_bulge)
    
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
    
    if snapshot:
        plot_number = 0
        last_plot_time = 0 | units.Myr
        make_plot_galstars(disk1, bulge1, stars, 
                           "MW and Solar System\nt = 0 Myr", 
                           'galstars_' + str(plot_number).zfill(4))
        
    x = [] | units.kpc
    y = [] | units.kpc
    
    current_iter = 0
    total_iter = int(t_end/interval) + 1
    t_end_scalar = t_end.value_in(units.Myr)
    
    times = np.arange(0., t_end_scalar, interval.value_in(units.Myr)) | units.Myr
    
    widgets = ['Step ', pbwg.SimpleProgress(), ' ',
               pbwg.Bar(marker='=', tip='>', left='[', right=']', fill=' '), 
               pbwg.Percentage(), ' - ', pbwg.ETA('ETA'), pbwg.EndMsg()]
    progress = pbar.ProgressBar(widgets=widgets, maxval=total_iter, fd=sys.stdout).start()
    
    for time in times:
        
        solver.evolve_model(time)
        mw_channel.copy()
        
        if leapfrog:
            star_solver(current_iter, stars, galaxy_dynamics_code)
        else:
            ch_g2l.copy()
        
        x.append(stars.x)
        y.append(stars.y)
        
        if snapshot:
            if check_last_plot_time(solver.model_time, last_plot_time, t_end/snap_freq):
                halo1, disk1, bulge1 = galaxy_structures(set1, n_disk, n_bulge)
                plot_number += 1
                last_plot_time = solver.model_time
                make_plot_galstars(disk1, bulge1, stars, 
                                   "MW and Solar System\nt = {} Myr".format(np.round(solver.model_time.value_in(units.Myr),
                                                                                     decimals=0)), 
                                   'galstars_' + str(plot_number).zfill(4))
        
        current_iter +=1
        
        progress.update(current_iter)
        
    progress.finish()
    
    galaxy_dynamics_code.stop()
    if not leapfrog:
        star_dynamics_code.stop()
