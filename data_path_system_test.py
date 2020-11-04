from modules import galaxies as gal
from amuse.lab import *

"""
import sys
from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model

import modules.progressbar as pbar
import modules.progressbar.widgets as pbwg
"""

#galaxy parameters
scale_mass_galaxy = 1e12 | units.MSun
scale_radius_galaxy = 80 | units.kpc
n_bulge = 1000
n_disk = 2000
n_halo = 4000

mw_parameters = {'name': 'mw',
                 #halo parameters
                 'n_halo': n_halo,
                 'halo_scale_length': 12.96 | units.kpc,
                 #disk parameters
                 'disk_number_of_particles' : n_disk,
                 'disk_mass' : 19.66 * 2.33 * 10e9 | units.MSun,
                 'disk_scale_length' : 2.806 | units.kpc,
                 'disk_outer_radius' : 30 | units.kpc, 
                 'disk_scale_height_sech2' : 0.409 | units.kpc,
                 'disk_central_radial_velocity_dispersion': 0.7,
                 #bulge parameters
                 'bulge_scale_radius' : 0.788 | units.kpc,
                 'bulge_number_of_particles' : n_bulge,
                 #unused parameters (unclear effect)
                 "halo_streaming_fraction": 0.,
                 "bulge_streaming_fraction": 0.,
                 #unused parameters (they cause problems)
                 'disk_scale_length_of_sigR2': 2.806 | units.kpc}

m31_parameters = {'name': 'm31_not_displaced',
                  #halo parameters
                  'n_halo': n_halo,
                  'halo_scale_length': 12.94 | units.kpc,
                  #disk parameters
                  'disk_number_of_particles' : n_disk,
                  'disk_mass' : 33.40 * 2.33 * 10e9 | units.MSun,
                  'disk_scale_length' : 5.577 | units.kpc,
                  'disk_outer_radius' : 30 | units.kpc, 
                  'disk_scale_height_sech2' : 0.3 | units.kpc,
                  'disk_central_radial_velocity_dispersion': 0.7,
                  #bulge parameters
                  'bulge_scale_radius' : 1.826 | units.kpc,
                  'bulge_number_of_particles' : n_bulge,
                  #unused parameters (unclear effect)
                  "halo_streaming_fraction": 0.,
                  "bulge_streaming_fraction": 0.,
                  #unused parameters (they cause problems)
                  'disk_scale_length_of_sigR2': 5.577 | units.kpc}

converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)

glxy, glxy_path = gal.test_make_galaxy(n_halo, converter, mw_parameters['name'], test=True,
                            disk_number_of_particles = mw_parameters['disk_number_of_particles'],
                            bulge_number_of_particles = mw_parameters['bulge_number_of_particles'])
print(glxy_path, flush=True)



"""
def make_galaxy():
    widgets = ['Building {} galaxy: '.format('test'), pbwg.AnimatedMarker(), ' ',
           pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        galaxy1 = new_galactics_model(n_halo,
                                    converter,
                                    do_scale=True,
                                    bulge_number_of_particles=n_bulge,
                                    disk_number_of_particles=n_disk,
                                    output_directory = '/data1/brentegani/')
        
make_galaxy()


import sys
from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model

import progressbar as pbar
import progressbar.widgets as pbwg

n_halo  = 4000
n_bulge = 2000
n_disk  = 1000
M_galaxy = 1e12 | units.MSun
R_galaxy = 80  | units.kpc
converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
widgets = ['Building {} galaxy: '.format('test'), pbwg.AnimatedMarker(), ' ',
           pbwg.Timer(), pbwg.EndMsg()]
with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
    galaxy1 = new_galactics_model(n_halo,
                                converter,
                                do_scale=True,
                                bulge_number_of_particles=n_bulge,
                                disk_number_of_particles=n_disk,
                                output_directory = '/data1/brentegani/')
galaxy1.move_to_center()
"""