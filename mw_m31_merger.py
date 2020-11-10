###### importing modules ######

import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '' and os.name != 'nt':
    print('No display found. Using non-interactive Agg backend for matplotlib')
    mpl.use('Agg')

import numpy as np
import sys

import modules.galaxies as gal
import modules.simulations as sim

from amuse.lab import units, Particles, nbody_system

"""
import inspect


#finds script path
filename = inspect.getframeinfo(inspect.currentframe()).filename
SCRIPT_PATH = os.path.dirname(os.path.abspath(filename))


plot_folders = [SCRIPT_PATH + '/plots', 
                SCRIPT_PATH + '/data/merger_plots',
                SCRIPT_PATH + '/plots/solar-system_plots',
                SCRIPT_PATH + '/data/merger_contour/',
                SCRIPT_PATH + '/plots/zoomed-merger-plots/']

for folder in plot_folders:
    if not os.path.exists(folder):
        os.makedirs(folder)
    
gal_folder = SCRIPT_PATH + '/data/'

if not os.path.exists(gal_folder):
    os.makedirs(gal_folder)


from amuse.lab import *
#from amuse.ext.galactics_model import new_galactics_model
#from amuse.couple import bridge
"""

#ignore warnings (AmuseWarning at end of simulation)
import warnings

warnings.filterwarnings('ignore')


#defines parser for terminal usage
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--generation', 
                    help='Generate a new galaxy model', 
                    action='store_true')
parser.add_argument('--plot', 
                    help='Allow output plots', 
                    action='store_true')
parser.add_argument('--solar', 
                    help='Add solar system', 
                    action='store_true')
parser.add_argument('--disk', 
                    help='Add test disk to MW', 
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
parser.add_argument('--correction', 
                    help='Correct galaxy velocity and mass', 
                    action='store_true')
parser.add_argument('-f', 
                    help='Foo value, defined for jupyter notebook compatibility')
args = parser.parse_args()

GENERATION = args.generation
PLOT = args.plot
SOLAR = args.solar
DISK = args.disk
IGM = args.igm
TEST = args.test
NOMERGER = args.nomerger
CORRECTION = args.correction

if PLOT:
    print('Plots turned on', flush=True)

    
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

#simulation parameters
scale_mass_galaxy = 1e12 | units.MSun
scale_radius_galaxy = 80 | units.kpc
t_end = 8500 | units.Myr

#Solar system starting conditions
n_stars = 10                                       #How many particles we will add
solar_radial_distance = 8.2
solar_position = (-np.sqrt(0.5 * solar_radial_distance**2),
                  np.sqrt(0.5 * solar_radial_distance**2),
                  0) | units.kpc                 #If we displace MW, we can do it through this vector aswell
system_radius = 0.100 | units.kpc                #neighbourhood in which we distribute our stars
solar_tang_velocity = 220 | (units.km/units.s)   #This is roughly the velocity of solarsystem around MW

#Intergal medium starting conditions
N1 = 5000
N2 = 1000
L = 10 | units.kpc
rho = 1000 | units.MSun / (units.kpc)**3
u = 1.6e+15 | (units.m)**2 / (units.s)**2

#M31 displacement
rotation = np.array([[0.7703,  0.3244,  0.5490],
                     [-0.6321, 0.5017,  0.5905],
                     [-0.0839, -0.8019, 0.5915]])

traslation = [-379.2, 612.7, 283.1] | units.kpc
radial_velocity = 117 * np.array([0.4898, -0.7914, 0.3657]) | units.kms
transverse_velocity = 50 * np.array([0.5236, 0.6024, 0.6024]) | units.kms


###### galaxy initialization ######

converter = nbody_system.nbody_to_si(scale_mass_galaxy, scale_radius_galaxy)

if GENERATION:
    print('Generating new galaxies', flush=True)
    mw, _ = gal.make_galaxy(converter, mw_parameters, test=TEST)
    m31_not_displaced, _ = gal.make_galaxy(converter, m31_parameters, test=TEST)
else:
    mw, _ = gal.load_galaxy_data('mw', test=TEST)
    m31_not_displaced, _ = gal.load_galaxy_data('m31_not_displaced', test=TEST)

"""
if not TEST:
    mw_data_path = SCRIPT_PATH + '/data/{}_full'.format(mw_parameters['name'])
    m31_not_displaced_data_path = SCRIPT_PATH + '/data/{}_full'.format(m31_parameters['name'])
    
elif TEST:
    mw_data_path = SCRIPT_PATH + '/data/{}_test'.format(mw_parameters['name'])
    m31_not_displaced_data_path = SCRIPT_PATH + '/data/{}_test'.format(m31_parameters['name']) 

if os.path.exists(mw_data_path) and os.path.exists(m31_not_displaced_data_path):
    widgets = ['Found galaxies data, loading: ', pbwg.AnimatedMarker(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        mw = read_set_from_file(mw_data_path, "hdf5")
        m31_not_displaced = read_set_from_file(m31_not_displaced_data_path, 'hdf5')
else:
    mw = gal.make_galaxy(converter, mw_parameters, SCRIPT_PATH, test=TEST)
    m31_not_displaced = gal.make_galaxy(converter, m31_parameters, SCRIPT_PATH, test=TEST)
"""

if NOMERGER:
    print('Quitting after galaxy initialization')
    quit()

if CORRECTION:
    print('Correcting velocities and mass: ...', flush=True)
    vel_factor = 1/6.5
    mass_factor = 1/1000

    mw.velocity = mw.velocity * vel_factor
    m31_not_displaced.velocity =  m31_not_displaced.velocity * vel_factor

    mw.mass = mw.mass * mass_factor
    m31_not_displaced.mass =  m31_not_displaced.mass * mass_factor


if all(value == False for value in [SOLAR, DISK, IGM]):
    m31 = gal.displace_galaxy(m31_not_displaced, rotation, traslation, radial_velocity, transverse_velocity)
    print('Simulating merger with no additional components (t = {} Myr) ...'.format(int(np.round(t_end.value_in(units.Myr), 
                                                                                                 decimals=0))), flush=True)
    sim.simulate_merger(mw, m31, n_halo, n_disk, n_bulge, t_end, converter, 
                        interval=0.5|units.Myr, plot=PLOT, plot_freq=1000)
    
if DISK:
    print('Simulating merger with disk test particles ...', flush=True)
    test_disk = gal.test_particles(MW, n_halo, n_bulge, n_disk)
    gal.simulate_merger_with_particles(M31, MW, converter, n_halo, n_bulge, n_disk, t_end, SCRIPT_PATH, plot=PLOT)

if SOLAR:
    print('Simulating MW with solar system ...', flush=True)
    MW.position += [100.0, 100.0, 0] | units.kpc
    t_end = 400 | units.Myr
    #MW.rotate(-np.pi/4, -np.pi/4, 0.0)
    stars = sol.make_solar_system(n_stars, solar_position, system_radius, MW_velocity_vector, solar_tang_velocity)
    gal.mw_and_stars(MW, stars, converter, sol.leapfrog_alg, n_halo, t_end, SCRIPT_PATH, plot=PLOT)


if IGM:
    print('Simulating merger with IGM ...', flush=True)
    widgets = ['Building IGM: ', pbwg.AnimatedMarker(), ' ',
               pbwg.Timer(), pbwg.EndMsg()]
    with pbar.ProgressBar(widgets=widgets, fd=sys.stdout) as progress:
        sph_code = igm.setup_sph_code(Gadget2, N1, N2, L, rho, u)

    gal.merger_and_igm(galaxy1, galaxy2, converter, sph_code, n_halo, t_end, SCRIPT_PATH, plot=PLOT)

