###### importing modules ######

from modules.common import __TEST_MODEL_DIR__, __FULL_MODEL_DIR__
from modules.progressbar import progressbar as pbar
from modules.progressbar import widgets as pbwg

import os
import sys

import datetime

from amuse.lab import units, Particles, write_set_to_file, read_set_from_file
from amuse.ext.galactics_model import new_galactics_model


###### make model output dir ######

def create_output_dir(glxy_name, test):
    if test == True:
        parent = __TEST_MODEL_DIR__
        tf = '_test_'
    elif test == False:
        parent = __FULL_MODEL_DIR__
        tf = '_full_'
    
    increasing = 1
    
    while True:
        current_model = glxy_name + tf + str(datetime.date.today()) + '-' + str(increasing).zfill(4)
        out_dir = parent + current_model + '/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            break
        else:
            increasing += 1
            
    return out_dir, current_model, tf


###### load galaxy data ######

def load_galaxy_data(glxy_name, test=False, loaded='last'):
    if test == True:
        parent = __TEST_MODEL_DIR__
        tf = '_test_'
    elif test == False:
        parent = __FULL_MODEL_DIR__
        tf = '_full_'
    
    dirnames_list = []
    mw_dirnames_list = []
    m31_dirnames_list = []
    
    for (dirpath, dirnames, filenames) in os.walk(parent):
        #print(list(filenames), flush=True)
        dirnames_list.extend(dirnames)
        break
    
    #print(os.walk(parent), flush=True)
    
    for dirname in dirnames_list:
        if dirname[0:2] == 'mw':
            mw_dirnames_list.append(dirname)
        elif dirname[0:3] == 'm31':
            m31_dirnames_list.append(dirname)
            
    if glxy_name == 'mw':
        current_dirlist = mw_dirnames_list
    elif glxy_name == 'm31_not_displaced':
        current_dirlist = m31_dirnames_list
        
    if loaded == 'last':
        try:
            glxy_dir = sorted(current_dirlist)[-1]
        except:
            try:
                glxy_dir = sorted(current_dirlist)[0]
            except:
                glxy_dir = ''
    else:
        for dirname in current_dirlist:
            if loaded == dirname:
                glxy_dir = dirname
                break
    
    
    glxy_data_name = glxy_dir
    glxy_data_dir = parent + glxy_dir + '/'
    glxy_data_path = glxy_data_dir + glxy_data_name

    if os.path.exists(glxy_data_path):
        widgets = ['Found galaxy data in {}, loading: '.format(glxy_data_dir), pbwg.AnimatedMarker(), pbwg.EndMsg()]
        with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
            glxy = read_set_from_file(glxy_data_path, "hdf5")
    else:
        raise ValueError('No galaxy model found')
            
    return glxy, glxy_data_dir


###### galaxy functions ######

def make_galaxy(n_halo, converter, glxy_name, test=False, **kwargs):
    widgets = ['Building {} galaxy: '.format(glxy_name), pbwg.AnimatedMarker(), ' ',
           pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy = new_galactics_model(n_halo,
                                     converter,
                                     do_scale=True,
                                     output_directory = '/data1/brentegani/',
                                     **kwargs)
        
    out_path, current_model, tf = create_output_dir(glxy_name, test=test)
    galaxy_data_path = out_path + current_model
    
    widgets = ['Saving {} galaxy data in {}: '.format(glxy_name, out_path), 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        write_set_to_file(galaxy, galaxy_data_path, 'hdf5')

    return galaxy, out_path


"""
def make_galaxy(converter, galaxy_dict, test=False):
    
    n_halo = galaxy_dict['n_halo']
    
    widgets = ['Building {} galaxy: '.format(galaxy_dict['name']), pbwg.AnimatedMarker(), ' ',
               pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy = new_galactics_model(n_halo,
                                     converter,
                                     do_scale=True,
                                     #halo parameters
                                     #halo_outer_radius = galaxy_dict['halo_outer_radius'],
                                     #halo_scale_length = galaxy_dict['halo_scale_length'],
                                     #disk parameters
                                     disk_number_of_particles = galaxy_dict['disk_number_of_particles'],
                                     disk_mass = galaxy_dict['disk_mass'],
                                     disk_scale_length = galaxy_dict['disk_scale_length'],
                                     disk_outer_radius = galaxy_dict['disk_outer_radius'],
                                     disk_scale_height_sech2 = galaxy_dict['disk_scale_height_sech2'],
                                     disk_central_radial_velocity_dispersion = galaxy_dict['disk_central_radial_velocity_dispersion'],
                                     #bulge paramaters
                                     bulge_scale_radius = galaxy_dict['bulge_scale_radius'],
                                     bulge_number_of_particles = galaxy_dict['bulge_number_of_particles'],
                                     #unused parameters
                                     #disk_scale_length_of_sigR2 = galaxy_dict['disk_scale_length_of_sigR2'],
                                     #halo_streaming_fraction = galaxy_dict["halo_streaming_fraction"],
                                     #bulge_streaming_fraction = galaxy_dict["bulge_streaming_fraction"],
                                     #output_directory = '/data1/brentegani/'
                                     )
    
    out_path, current_model, tf = create_output_dir(galaxy_dict['name'], test=test)
    galaxy_data_path = out_path + current_model
    
    widgets = ['Saving {} galaxy data: '.format(galaxy_dict['name']), 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        write_set_to_file(galaxy, galaxy_data_path, 'hdf5')

    return galaxy, galaxy_data_path
"""

###### galaxy rototranslation ######

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


###### test disk function ######

def test_particles(galaxy, n_halo, n_bulge, n_disk):
    _dsk = Particles(len(galaxy))
    _dsk.mass = galaxy.mass
    _dsk.position = galaxy.position
    _dsk.velocity = galaxy.velocity
    disk = _dsk[n_bulge:n_halo]
    return disk