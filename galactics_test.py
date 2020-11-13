###### importing modules ######

from modules.progressbar import progressbar as pbar
from modules.progressbar import widgets as pbwg

import sys

from amuse.lab import units, Particles, write_set_to_file, nbody_system
from amuse.ext.galactics_model import new_galactics_model


###### mw starting conditions ######

scale_mass_galaxy = 1e12 | units.MSun
scale_radius_galaxy = 80 | units.kpc

mw_param = {'name': 'mw',
            #halo parameters
            'n_halo': 40000,
            'halo_scale_radius': 12.96 | units.kpc,
            #disk parameters
            'disk_number_of_particles' : 20000,
            'disk_mass' : 19.66 * 2.33 * 10e9 | units.MSun,
            'disk_scale_length' : 2.806 | units.kpc,
            'disk_outer_radius' : 30 | units.kpc, 
            'disk_scale_height_sech2' : 0.409 | units.kpc,
            'disk_scale_length_of_sigR2': 2.806 | units.kpc,
            'disk_central_radial_velocity_dispersion': 0.7,
            #bulge parameters
            'bulge_scale_radius' : 0.788 | units.kpc,
            'bulge_number_of_particles' : 10000}


###### galaxy function ######

def make_galaxy_TEST(n_halo, converter, glxy_name, **kwargs):
    widgets = ['Building {} galaxy: '.format(glxy_name), pbwg.AnimatedMarker(), ' ',
           pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy = new_galactics_model(n_halo,
                                     converter,
                                     do_scale=True,
                                     verbose=True,
                                     **kwargs)
    
    widgets = ['Saving {} galaxy data: '.format(glxy_name), 
               pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        write_set_to_file(galaxy, glxy_name + 'TEST', 'hdf5')

    return galaxy


###### main ######

converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)

glxy = make_galaxy_TEST(mw_param['n_halo'], converter, mw_param['name'],
                        #halo parameters
                        halo_scale_radius = mw_param['halo_scale_radius'],
                        #disk parameters
                        disk_number_of_particles = mw_param['disk_number_of_particles'],
                        disk_mass = mw_param['disk_mass'],
                        disk_scale_length = mw_param['disk_scale_length'],
                        disk_outer_radius = mw_param['disk_outer_radius'],
                        disk_scale_height_sech2 = mw_param['disk_scale_height_sech2'],
                        disk_central_radial_velocity_dispersion = mw_param['disk_central_radial_velocity_dispersion'],
                        #bulge paramaters
                        bulge_scale_radius = mw_param['bulge_scale_radius'],
                        bulge_number_of_particles = mw_param['bulge_number_of_particles'])