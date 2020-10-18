#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units
from amuse.lab import Particles
from amuse.units import nbody_system
from amuse.community.ph4.interface import ph4


# In[2]:


def MW_and_M31():
    MW_M31 = Particles(2)
    MW = MW_M31[0]
    MW.mass = 1*10**12 | units.MSun
    MW.position = (0,0,0) | units.kpc
    MW.velocity = (10,10,0) | units.kms
    M31 = MW_M31[1]
    M31.mass = 1.6*10**12 | units.MSun
    M31.position = (780,0,0) | units.kpc
    M31.velocity = (0,0,0) | units.kms
    return MW_M31


# In[3]:


def merge_two_particles(gravity, particles_in_encounter):
    new_particle = Particles(1)
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.position = particles_in_encounter.center_of_mass()
    new_particle.velosity = particles_in_encounter.center_of_mass_velocity()
    new_particle.radius = particles_in_encounter.radius.sum()
    gravity.particles.add_particles(new_particle)
    gravity.particles.remove_particles(particles_in_encounter)


# In[4]:


def resolve_collision(collision_detection, gravity, bodies):
    if collision_detection.is_set():
        for ci in range(len(collision_detection.particles(0))): 
            encountering_particles = Particles(particles=[collision_detection.particles(0)[ci],
                                                          collision_detection.particles(1)[ci]])
            colliding_particles = encountering_particles.get_intersecting_subset_in(bodies)
            merge_two_particles(bodies, colliding_particles)
            bodies.synchronize_to(gravity.particles)


# In[5]:


galaxies = MW_and_M31()
converter=nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
gravity = ph4(converter)
gravity.particles.add_particles(galaxies)#add particles

stopping_condition = gravity.stopping_conditions.collision_detection
stopping_condition.enable()

end_time = 10.0 | units.Myr
model_time = 0 | units.Myr

while(model_time<end_time):
    dt = 1.0 | units.Myr
    model_time += dt
    gravity.evolve_model(model_time)
    resolve_collision(stopping_condition, gravity, galaxies)
    plt.scatter(galaxies.x.value_in(units.kpc),galaxies.y.value_in(units.kpc),c = 'r')

    

gravity.stop()

