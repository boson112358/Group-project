import os
from data_analysis import __SCRIPT_PATH__

#creates plot folder
PLOT_FOLDER = __SCRIPT_PATH__ + '/plots/model-analysis/'
if not os.path.exists(PLOT_FOLDER):
    os.makedirs(PLOT_FOLDER)

import numpy as np
import matplotlib.pyplot as plt
import random, math

from amuse.lab import units, Particle


###### utility functions ######

def galaxy_total_mass(galaxy):
    print('Galaxy total mass: ' + str(np.sum(galaxy.mass).in_(units.MSun)))   #Check total mass, should be around 1.5e12 MSun
    return(np.sum(galaxy.mass).in_(units.MSun))


def galaxy_structures(galaxy, n_disk, n_bulge, n_halo=None):
    halo = galaxy[n_disk+n_bulge:]
    disk = galaxy[:n_disk]
    bulge = galaxy[n_disk:n_disk+n_bulge]
    return halo, disk, bulge


def position_limit(structure, axes=[0, 1]):
    max_coordinate = 0
    for coordinate in axes:
        _temp = np.amax(structure.position[coordinate].value_in(units.kpc))
        if _temp > max_coordinate:
            max_coordinate = _temp
    return max_coordinate


def average_velocity_at_radius(radius, velocity, interval_length=1/25, max_radius=30):
    interval_limit = np.arange(0, max_radius + interval_length, step=interval_length)
    average_vel = []
    for i in range(1, len(interval_limit), 1):
        bin_l = interval_limit[i-1]
        bin_r = interval_limit[i]
        bin_velocity = []
        for index, value in enumerate(radius):
            if value > bin_l and value <= bin_r:
                bin_velocity.append(velocity[index])
        bin_average = np.mean(bin_velocity)
        average_vel.append(bin_average)
    
    intervals = np.delete(interval_limit, -1) + (1/2 * interval_length)  #the point is in the middle of the interval
    
    return intervals, average_vel

###### plot functions ######
    
def plot_galaxy_structure(structure, filename, title=None, label=None):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    plot_limit = position_limit(structure) + 10
    
    if title == None:
        title = filename
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-plot_limit, plot_limit)
    plt.ylim(-plot_limit, plot_limit)

    ax.scatter(structure.x.value_in(units.kpc), structure.y.value_in(units.kpc),
               c='tab:blue', alpha=1, s=1, lw=0, label=label)
    
    if label == None:
        pass
    else:
        plt.legend(loc='upper right')
  
    plt.savefig(PLOT_FOLDER + filename)
    
    
def plot_velocity_component(distance, component, filename, title=None, label=None):
    x_label = "Distance from galaxy center of mass [kpc]"
    y_label = "Velocity [km/s]"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    x_max = np.amax(distance) + 0.1
    y_max = np.amax(component) + 10
    y_min = np.amin(component) - 10
    
    if title == None:
        title = filename
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(-0.1, x_max)
    plt.ylim(y_min, y_max)

    ax.scatter(distance, component,
               c='tab:blue', alpha=1, s=1, lw=0, label=label)
    
    if label == None:
        pass
    else:
        plt.legend(loc='upper right')
  
    plt.savefig(PLOT_FOLDER + filename)
    
    
def structure_rotation_curve(ax, distance, total_velocity, label=None, color='tab:blue'):
    intervals, average_vel = average_velocity_at_radius(distance, total_velocity)

    ax.scatter(intervals, average_vel, 
               c=color, alpha=1, s=5, label=label)
    
    
def galaxy_rotation_curve(structure_distances, structure_velocities, filename, title=None, labels=None):
    x_label = "Distance from galaxy center of mass [kpc]"
    y_label = "Velocity [km/s]"
    
    if labels == None:
        labels = [None for i in range(len(structure_distances))]
    
    default_colors = ['tab:blue', 'tab:orange', 'tab:green']
    used_colors = [default_colors[i] for i in range(len(structure_distances))]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if labels == None:
        labels = [None for i in range(len(structure_distances))]
    
    for distance, velocity, c, label, in zip(structure_distances, structure_velocities, used_colors, labels):
        structure_rotation_curve(ax, distance, velocity, label=label, color=c)
    
    #x_max = np.amax(distance) + 0.1
    #y_max = np.amax(total_velocity) + 10
    #y_min = np.amin(total_velocity) - 10
    
    if title == None:
        title = filename
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    #plt.xlim(-0.1, x_max)
    #plt.ylim(y_min, y_max)
    
    if label == None:
        pass
    else:
        plt.legend(loc='upper right')
  
    plt.savefig(PLOT_FOLDER + filename)
    
    
###### velocity functions ######
    
def velocity_components(body, galaxy_com):
    unit_v = 1 | units.m**-1 * units.s
    unit_m = 1 | units.m**-1
    
    body_displaced = body.position - galaxy_com   #Displace since we are considering it with relation to the center of mass
            
    if body_displaced.z < 0 | units.m:
        theta = math.atan(math.sqrt((body_displaced.x**2 + body_displaced.y**2)*unit_m**2)/(unit_m*body_displaced.z)) + math.pi
    else:
        theta = math.atan(math.sqrt((body_displaced.x**2 + body_displaced.y**2)*unit_m**2)/(unit_m*body_displaced.z))

    if body_displaced.x < 0 | units.m:
        phi = math.atan(body_displaced.y / body_displaced.x) + math.pi
    else:
        phi = math.atan(body_displaced.y / body_displaced.x)

    runit = (math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta))
    tunit = (math.cos(theta)*math.cos(phi), math.cos(theta)*math.sin(phi), -math.sin(theta))
    punit = (-math.sin(phi), math.cos(phi), 0)

    z1 = np.array([[runit[0],tunit[0],punit[0]], 
                   [runit[1],tunit[1],punit[1]], 
                   [runit[2],tunit[2],punit[2]]])

    z2 = np.array([unit_v * body.velocity.x, 
                   unit_v * body.velocity.y, 
                   unit_v * body.velocity.z])

    polar_v = 1/1000 * np.linalg.solve(z1, z2)

    com_distance = math.sqrt((body_displaced.x**2 + body_displaced.y**2 + body_displaced.z**2)*unit_m**2) | units.m
    total_v = math.sqrt(unit_v**2*(body.velocity.x**2 + body.velocity.y**2 + body.velocity.z**2)) * 1/1000 | units.km/units.s
    
    radial_v = polar_v[0] | units.km/units.s        #x component
    tang_v = polar_v[1] | units.km/units.s          #y component
    angular_v = polar_v[2] | units.km/units.s       #z component

    return com_distance, total_v, radial_v, tang_v, angular_v    
    

def galaxy_structure_velocity(galaxy, n_disk, n_bulge):
    com = galaxy.center_of_mass()
    halo, disk, bulge = galaxy_structures(galaxy, n_disk, n_bulge)
    
    halo_dict = {}
    disk_dict = {}
    bulge_dict = {}
    
    for structure, dictionary in zip([halo, disk, bulge], [halo_dict, disk_dict, bulge_dict]):
        com_dist_list = []    #| units.kpc
        angular_v_list = []   #| units.km/units.s
        radial_v_list = []    #| units.km/units.s
        tang_v_list = []      #| units.km/units.s
        total_v_list = []     #| units.km/units.s
        
        for body in structure:
            com_dist, total_v, radial_v, tang_v, angular_v = velocity_components(body, com)
            
            com_dist_list.append(com_dist.value_in(units.kpc))
            angular_v_list.append(angular_v.value_in(units.km/units.s))
            radial_v_list.append(radial_v.value_in(units.km/units.s))
            tang_v_list.append(tang_v.value_in(units.km/units.s))
            total_v_list.append(total_v.value_in(units.km/units.s))
        
        dictionary.update({'com_distance': com_dist_list})
        dictionary.update({'angular_velocity': angular_v_list})
        dictionary.update({'radial_velocity': radial_v_list})
        dictionary.update({'tangential_velocity': tang_v_list})
        dictionary.update({'total_velocity': total_v_list})
        
    return halo_dict, disk_dict, bulge_dict