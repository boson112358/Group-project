###### importing modules ######

import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '' and os.name != 'nt':
    print('No display found. Using non-interactive Agg backend for matplotlib')
    mpl.use('Agg')

import numpy as np
import sys
import argparse

import modules.galaxies as gal
import modules.simulations as sim
import modules.data_analysis as da
import modules.solar_system as sol
import modules.igmedium as igm
from modules.g2extension import Gadget2Gravity
from modules.progressbar import progressbar as pbar
from modules.progressbar import widgets as pbwg

from amuse.lab import units, Particles, nbody_system, read_set_from_file, Fi, Gadget2


###### parser for terminal usage ######

parser = argparse.ArgumentParser()
parser.add_argument('--generation', 
                    help='Generate a new galaxy model', 
                    action='store_true')
parser.add_argument('--animation', 
                    help='Allow animation output', 
                    action='store_true')
parser.add_argument('--snapshot', 
                    help='Allow simulation snapshots plots output', 
                    action='store_true')
parser.add_argument('--solar', 
                    help='Add solar system', 
                    action='store_true')
parser.add_argument('--mwsolar', 
                    help='Simulates only the MW with the Solar System', 
                    action='store_true')
parser.add_argument('--igm', 
                    help='Add IGM', 
                    action='store_true')
parser.add_argument('--test', 
                    help='Use test initial conditions', 
                    action='store_true')
parser.add_argument('--nomerger', 
                    help='Stop after galaxy initialization', 
                    action='store_true')
parser.add_argument('--radvel', 
                    help='Value of the m31 radial velocity factor, default: 1', 
                    nargs=1)
parser.add_argument('--transvel', 
                    help='Value of the m31 transverse velocity factor, default: 1', 
                    nargs=1)
parser.add_argument('--timefinal', 
                    help='End time of the simulation, units are in Myr, default: 15000', 
                    nargs=1)
parser.add_argument('--timestep', 
                    help='Timestep of the simulation, units are in Myr, default: 5', 
                    nargs=1)
parser.add_argument('-f', 
                    help='Foo value, defined for jupyter notebook compatibility')
args = parser.parse_args()

GENERATION = args.generation
ANIMATION = args.animation
SNAPSHOT = args.snapshot
SOLAR = args.solar
MWSOLAR = args.mwsolar
IGM = args.igm
TEST = args.test
NOMERGER = args.nomerger

if ANIMATION:
    print('Animations turned on', flush=True)
if SNAPSHOT:
    print('Merger snapshot plots turned on', flush=True)

    
###### starting conditions ######
    
#Galaxies starting conditions
if not TEST:
    n_bulge = 20000
    n_disk = 20000
    n_halo = 40000
elif TEST:
    print('Using test initial conditions', flush=True)
    n_bulge = 10000
    n_disk = 20000
    n_halo = 40000

mw_parameters = {'name': 'mw',
                 #halo parameters
                 'n_halo': n_halo,
                 'halo_scale_radius': 12.96 | units.kpc,
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
                  'halo_scale_radius': 12.94 | units.kpc,
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

#simulation parameters
scale_mass_galaxy = 1e12 | units.MSun
scale_radius_galaxy = 80 | units.kpc

#parses the simulation times
if args.timefinal == None:
    t_end = 15000 | units.Myr
else:
    t_end = float(args.timefinal[0]) | units.Myr
    
if args.timestep == None:
    t_step = 5. | units.Myr
else:
    t_step = float(args.timestep[0]) | units.Myr

#Solar system starting conditions
n_stars = 1000                                       #How many particles we will add
solar_radial_distance = 8.2
solar_position = (-np.sqrt(0.5 * solar_radial_distance**2),
                  np.sqrt(0.5 * solar_radial_distance**2),
                  0) | units.kpc                   #If we displace MW, we can do it through this vector aswell
system_radius = 0.100 | units.kpc                  #neighbourhood in which we distribute our stars
solar_tang_velocity = 220 | (units.km/units.s)     #This is roughly the velocity of solarsystem around MW
mw_velocity_vector = (0, 0, 0) | units.kms

#Intergal medium starting conditions
N1 = n_bulge + n_disk + n_halo
N2 = n_bulge + n_disk + n_halo
box_side = 1000 | units.kpc
box_grid = box_side.value_in(units.kpc)
rho = 770 | units.MSun / (units.kpc)**3
u = 3.724e+9 | (units.m)**2 / (units.s)**2

#M31 displacement
rotation = np.array([[0.7703,  0.3244,  0.5490],
                     [-0.6321, 0.5017,  0.5905],
                     [-0.0839, -0.8019, 0.5915]])

traslation = [-379.2, 612.7, 283.1] | units.kpc

#parses the velocity factors
if args.radvel == None:
    m31_radvel_factor = 1.
else:
    m31_radvel_factor = float(args.radvel[0])
    
if args.transvel == None:
    m31_transvel_factor = 1.
else:
    m31_transvel_factor = float(args.transvel[0])
    
radial_velocity =  m31_radvel_factor * 117 * np.array([0.4898, -0.7914, 0.3657]) | units.kms
transverse_velocity = m31_transvel_factor * 42 * np.array([0.5236, 0.6024, 0.6024]) | units.kms

initial_conditions_dict = {#mw parameters
                           'mw_n_halo': [mw_parameters['n_halo']],
                           'mw_disk_number_of_particles' : [mw_parameters['disk_number_of_particles']],
                           'mw_disk_mass (MSun)' :[mw_parameters['disk_mass'].value_in(units.MSun)],
                           'mw_bulge_number_of_particles' : [mw_parameters['bulge_number_of_particles']],
                           #m31 parameters
                           'm31_n_halo': [m31_parameters['n_halo']],
                           'm31_disk_number_of_particles' : [m31_parameters['disk_number_of_particles']],
                           'm31_disk_mass (MSun)' : [m31_parameters['disk_mass'].value_in(units.MSun)],
                           'm31_bulge_number_of_particles' : [m31_parameters['bulge_number_of_particles']],
                           #time parameters
                           't final (Myr)': [t_end.value_in(units.Myr)],
                           't step (Myr)': [t_end.value_in(units.Myr)],
                           #m31 dynamic parameters
                           'm31_radvel_factor (* 117 km/s)': [m31_radvel_factor],
                           'm31_transvel_factor (* 42 km/s)': [m31_transvel_factor],
                           #solar tracker parameters
                           'add_solar': [SOLAR],
                           'n_stars': [n_stars],                                      
                           'solar_radial_distance (kpc)': [solar_radial_distance],
                           'solar_system_radius (kpc)': [system_radius.value_in(units.kpc)],
                           #igm parameters
                           'add_igm': [IGM]}

###### galaxy initialization ######

converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)

if GENERATION:
    print('Generating new galaxies', flush=True)
    mw, _ = gal.make_galaxy(mw_parameters['n_halo'], converter, mw_parameters['name'], test=TEST,
                            #output dir
                            output_directory = '/data1/brentegani/',
                            #halo parameters
                            #halo_scale_radius = glxy_param['halo_scale_radius'],
                            #disk parameters
                            disk_number_of_particles = mw_parameters['disk_number_of_particles'],
                            disk_mass = mw_parameters['disk_mass'],
                            #disk_scale_length = glxy_param['disk_scale_length'],
                            #disk_outer_radius = glxy_param['disk_outer_radius'],
                            #disk_scale_height_sech2 = glxy_param['disk_scale_height_sech2'],
                            #disk_central_radial_velocity_dispersion=glxy_param['disk_central_radial_velocity_dispersion'],
                            #bulge paramaters
                            #bulge_scale_radius = glxy_param['bulge_scale_radius'],
                            bulge_number_of_particles = mw_parameters['bulge_number_of_particles'])
    
    m31_not_displaced, _ = gal.make_galaxy(m31_parameters['n_halo'], converter, m31_parameters['name'], test=TEST,
                                           #output dir
                                           output_directory = '/data1/brentegani/',
                                           #halo parameters
                                           #halo_scale_radius = glxy_param['halo_scale_radius'],
                                           #disk parameters
                                           disk_number_of_particles = m31_parameters['disk_number_of_particles'],
                                           disk_mass = m31_parameters['disk_mass'],
                                           #disk_scale_length = glxy_param['disk_scale_length'],
                                           #disk_outer_radius = glxy_param['disk_outer_radius'],
                                           #disk_scale_height_sech2 = glxy_param['disk_scale_height_sech2'],
                                           #disk_central_radial_velocity_dispersion=m31_param['disk_central_radial_velocity_dispersion'],
                                           #bulge paramaters
                                           #bulge_scale_radius = glxy_param['bulge_scale_radius'],
                                           bulge_number_of_particles = m31_parameters['bulge_number_of_particles'])
else:
    #mw, _ = gal.load_galaxy_data('mw', test=TEST)
    #m31_not_displaced, _ = gal.load_galaxy_data('m31_not_displaced', test=TEST)
    mw_data_dir = 'used_models/mw_test_2020-11-10-0001/'
    m31_data_dir = 'used_models/m31_not_displaced_test_2020-11-10-0001/'
    
    widgets = ['Found galaxy data in {}, loading: '.format(mw_data_dir), pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        mw = read_set_from_file(mw_data_dir + 'mw_test_2020-11-10-0001', "hdf5")
            
    widgets = ['Found galaxy data in {}, loading: '.format(m31_data_dir), pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        m31_not_displaced = read_set_from_file(m31_data_dir + 'm31_not_displaced_test_2020-11-10-0001', "hdf5")
    
mw_mass = da.galaxy_total_mass(mw)
m31_mass = da.galaxy_total_mass(m31_not_displaced)

initial_conditions_dict.update({'mw_mass': [mw_mass.value_in(units.MSun)], 'm31_mass': [m31_mass.value_in(units.MSun)]})

if NOMERGER:
    print('Quitting after galaxy initialization')
    quit()


###### main ######

if not MWSOLAR:
    m31 = gal.displace_galaxy(m31_not_displaced, rotation, traslation, radial_velocity, transverse_velocity)

    t_end_int = int(np.round(t_end.value_in(units.Myr), decimals=0))
    t_step_int = int(np.round(t_step.value_in(units.Myr), decimals=0))

    if IGM:
        txt_line1 = 'Simulating merger with IGM\n'
        widgets = ['Building IGM: ', pbwg.AnimatedMarker(), ' ',
                   pbwg.Timer(), pbwg.EndMsg()]
        with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
            igm_gas, igm_dm = igm.setup_sph_code(N1, N2, box_side, rho, u)
    else:
        igm_gas = None
        igm_dm = None

    if SOLAR:
        txt_line1 = 'Simulating merger with Solar System (n = {})\n'.format(n_stars)
        stars = sol.make_solar_system(n_stars, solar_position, system_radius, mw_velocity_vector, solar_tang_velocity)
    else:
        txt_line1 = 'Simulating merger with no additional components:\n'
        stars = None

    if SOLAR and IGM:
        txt_line1 = 'Simulating merger with IGM and Solar System (n = {})\n'.format(n_stars)

    txt_line2 = 't = {} Myr, t step = {}\n'.format(t_end_int, t_step_int)
    txt_line3 = 'MW mass = {}, M31 mass = {}\n'.format(mw_mass, m31_mass)
    txt_line4 = 'm31 radial velocity factor = {} * 117\n'.format(m31_radvel_factor)
    txt_line5 = 'm31 transverse velocity factor = {} * 42'.format(m31_transvel_factor)
    print(txt_line1 + txt_line2 + txt_line3 + txt_line4 + txt_line5, flush=True)

    sim.simulate_merger(mw, m31, n_halo, n_disk, n_bulge, t_end, converter, Gadget2Gravity, initial_conditions_dict,
                        interval=t_step, 
                        animation=ANIMATION, snapshot=SNAPSHOT, snap_freq=1000,
                        sol_system=stars,  
                        igm_gas_particles=igm_gas, igm_dm_particles=igm_dm, box_grid=box_grid)
else:
    print('Simulating MW with solar system ...', flush=True)
    mw_velocity_vector = (0, 0, 0) | units.kms
    stars = sol.make_solar_system(n_stars, solar_position, system_radius, mw_velocity_vector, solar_tang_velocity)
    sim.mw_and_stars(mw, stars, converter, n_disk, n_bulge, t_end, sol.leapfrog_alg, snapshot=SNAPSHOT)
