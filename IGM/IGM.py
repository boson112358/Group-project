#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy
from amuse.lab import *
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
#
from amuse.units import units
from amuse.lab import Particles
from amuse.community.ph4.interface import ph4


# In[2]:


def setup_sph_code(sph_code, N1, N2, L, rho, u):
    #rho -- local group density
    converter = ConvertBetweenGenericAndSiUnits(L, rho, constants.G)
    sph_code = sph_code(converter, mode = 'periodic')#, redirection = 'none')
    sph_code.parameters.periodic_box_size = 10.0 | units.parsec
    dm = Particles(N1)
    dm.mass = (rho * L**3) / N1
    numpy.random.seed(12345)
    dm.x = L * numpy.random.uniform(0.0, 1.0, N1)
    dm.y = L * numpy.random.uniform(0.0, 1.0, N1)
    dm.z = L * numpy.random.uniform(0.0, 1.0, N1)
    dm.vx = numpy.zeros(N1) | units.km / units.s
    dm.vy = numpy.zeros(N1) | units.km / units.s
    dm.vz = numpy.zeros(N1) | units.km / units.s
    gas = Particles(N2)
    gas.mass = 0.0001*(rho * L**3) / N2
    gas.x = L * numpy.random.uniform(0.0, 1.0, N2)
    gas.y = L * numpy.random.uniform(0.0, 1.0, N2)
    gas.z = L * numpy.random.uniform(0.0, 1.0, N2)
    gas.vx = numpy.zeros(N2) | units.km / units.s
    gas.vy = numpy.zeros(N2) | units.km / units.s
    gas.vz = numpy.zeros(N2) | units.km / units.s
    gas.u = u
    if isinstance(sph_code, Fi):
        sph_code.parameters.self_gravity_flag = False
        sph_code.parameters.timestep = 0.1 | generic_unit_system.time
        gas.h_smooth = L / N**(1/3.0)
        gas.position -= 0.5 * L
        
    sph_code.gas_particles.add_particles(gas)
    sph_code.dm_particles.add_particles(dm)
    sph_code.commit_particles()
    return sph_code


# In[3]:


def MW_and_M31():
    MW_M31 = Particles(2)
    MW = MW_M31[0]
    MW.mass = 1*10**12 | units.Msun
    MW.position = (0,0,0) | units.kpc
    MW.velosity = (0,0,0) | units.kms
    M31 = MW_M31[1]
    M31.mass = 1.6*10**12 | units.Msun
    M31.position = (780,0,0) | units.kpc
    M31.velosity = (0,0,0) | units.kms
    


# In[4]:


def merge_two_particles(gravity, particles_in_encounter):
    new_particle = Particles(1)
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.position = particles_in_encounter.center_of_mass()
    new_particle.velosity = particles_in_encounter.center_of_mass_velocity()
    new_particle.radius = particles_in_encounter.radius.sum()
    gravity.particles.add_particles(new_particle)
    gravity.particles.remove_particles(particles_in_encounter)
    


# In[5]:


def resolve_collision(collision_detection, gravity, bodies):
    if collision_detection.is_set():
        for ci in range(len(collision_detection.particles(0))): 
            encountering_particles = Particles(particles=[collision_detection.particles(0)[ci],
                                                          collision_detection.particles(1)[ci]])
            colliding_particles = encountering_particles.get_intersecting_subset_in(bodies)
            merge_two_particles(bodies, colliding_particles)
            bodies.synchronize_to(gravity.particles)
            


# In[6]:


gravity = ph4()
gravity.particles.add_particles()#add particles

stopping_condition = gravity.stopping_conditions.collision_detection
stopping_condition.enable()

end_time = 10.0 | units.Myr
model_time = 0 | units.Myr

while(model_time<end_time):
    dt = gravity.particles.time_step.min()
    model_time += dt
    gravity.evolve_model(model_time)
    resolve_collision(stopping_condition, gravity, bodies)

gravity.stop()

