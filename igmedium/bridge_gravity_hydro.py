#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from amuse.plot import plot
from amuse.units import units
from amuse.lab import Particles
from amuse.units import nbody_system
from amuse.community.ph4.interface import ph4
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi


# In[2]:


def MW_and_M31():
    MW_M31 = Particles(2)
    MW = MW_M31[0]
    MW.mass = 1*10**6 | units.MSun
    MW.position = (0,0,0) | units.kpc
    MW.velocity = (0.1,0.1,0) | units.kms
    M31 = MW_M31[1]
    M31.mass = 1.6*10**6 | units.MSun
    M31.position = (7.8,0,0) | units.kpc
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
ch_g2l = gravity.particles.new_channel_to(galaxies)

stopping_condition = gravity.stopping_conditions.collision_detection
stopping_condition.enable()

end_time = 6600.0 | units.Myr
model_time = 0 | units.Myr
x = [] | units.kpc
y = [] | units.kpc

while(model_time<end_time):
    dt = 10.0 | units.Myr
    model_time += dt
    gravity.evolve_model(model_time)
    ch_g2l.copy()
    x.append(galaxies.x)
    y.append(galaxies.y)
    resolve_collision(stopping_condition, gravity, galaxies)


a = 0
b = 4

#plt.scatter(galaxies.x.value_in(units.kpc),galaxies.y.value_in(units.kpc),c = 'r')
plot(x, y, lw=1)
plt.gca().set_aspect("equal", adjustable="box")
plt.ylim((a,b))
plt.show()

    

gravity.stop()


# In[6]:


# from amuse.lab import BHTree
# galaxies = MW_and_M31()
# converter=nbody_system.nbody_to_si(galaxies.mass.sum(), galaxies[1].position.length())
# gravity = BHTree(converter)
# gravity.particles.add_particles(galaxies)#add particles
# ch_g2l = gravity.particles.new_channel_to(galaxies)
# stopping_condition = gravity.stopping_conditions.collision_detection
# stopping_condition.enable()
# timestep = 10 #units Myr
# end_time = 8000.0 #units Myr
# gravity.timestep = timestep |units.Myr
# times = numpy.arange(0., end_time, timestep) | units.Myr
# x = [] | units.kpc
# y = [] | units.kpc
# for time in times:
#     gravity.evolve_model(time)
#     ch_g2l.copy()
#     x.append(galaxies.x)
#     y.append(galaxies.y)
# gravity.stop()
# plot(x, y, lw=1)
# pyplot.gca().set_aspect("equal", adjustable="box")
# pyplot.show()


# In[7]:


def setup_sph_code(sph_code, N1, N2, L, rho,u):
    #rho -- local group density
    converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.kpc)
    sph_code = sph_code(converter, redirection = 'none')#, mode = 'periodic')#, redirection = 'none')#change periodic
    #sph_code.parameters.periodic_box_size = 10.0 | units.kpc
    dm = Particles(N1)
    dm.mass = (rho * L**3) / N1
    np.random.seed(12345)
    dm.x = L * np.random.uniform(0.0, 1.0, N1)
    dm.y = L * np.random.uniform(0.0, 1.0, N1)
    dm.z = L * np.random.uniform(0.0, 1.0, N1)
    dm.vx = np.zeros(N1) | units.km / units.s
    dm.vy = np.zeros(N1) | units.km / units.s
    dm.vz = np.zeros(N1) | units.km / units.s
    gas = Particles(N2)
    gas.mass = (rho * L**3) / N2
    gas.x = L * np.random.uniform(0.0, 1.0, N2)
    gas.y = L * np.random.uniform(0.0, 1.0, N2)
    gas.z = L * np.random.uniform(0.0, 1.0, N2)
    gas.vx = np.zeros(N2) | units.km / units.s
    gas.vy = np.zeros(N2) | units.km / units.s
    gas.vz = np.zeros(N2) | units.km / units.s
    gas.u = u 
    if isinstance(sph_code, Fi):
        sph_code.parameters.self_gravity_flag = True
        sph_code.parameters.timestep = 5 | units.Myr
        gas.h_smooth = L / N2**(1/3.0)
        dm.h_smooth = L / N1**(1/3.0)
        #gas.position -= 0.5 * L
        
    sph_code.gas_particles.add_particles(gas)
    sph_code.dm_particles.add_particles(dm)
    sph_code.commit_particles()
    return sph_code


# In[8]:


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


# In[9]:


N1 = 500
N2 = 100
L = 10 | units.kpc
Lg = 40
rho = 1000 | units.MSun / (units.kpc)**3
u = 1.6e+15 | (units.m)**2 / (units.s)**2
sph_code = setup_sph_code(Fi, N1, N2, L, rho, u)


# In[ ]:





# In[10]:


sph_code.evolve_model(10.0|units.Myr)
print("Done running")
plt.scatter(sph_code.gas_particles.x.value_in(units.kpc),
            sph_code.gas_particles.y.value_in(units.kpc),
            c = 'r')
plt.scatter(sph_code.dm_particles.x.value_in(units.kpc),
            sph_code.dm_particles.y.value_in(units.kpc),
            c = 'b')
#pyplot.imshow()


# In[11]:


gravity1 = ph4(converter)
gravity1.particles.add_particles(galaxies)
channel = {"from_galaxies": galaxies.new_channel_to(gravity1.particles),
"to_galaxies": gravity1.particles.new_channel_to(galaxies)}

channel.update({"from_gas": sph_code.gas_particles.new_channel_to(sph_code.particles)})
channel.update({"to_gas": sph_code.particles.new_channel_to(sph_code.gas_particles)})

galaxies.add_particles(sph_code.gas_particles)

from amuse.couple import bridge
from amuse.ext.composition_methods import *
gravhydro = bridge.Bridge(use_threading=False) #, method=SPLIT_4TH_S_M4)
gravhydro.add_system(gravity1, (sph_code,))
gravhydro.add_system(sph_code, (gravity1,))
gravhydro.timestep = 5.0 | units.Myr


# In[ ]:


def gravity_hydro_bridge(gravity, hydro, gravhydro, galaxies,
                         t_end):

    model_time = 0 | units.Myr
    dt = 10|units.Myr  #1.0*Pinner
    while model_time < t_end:

        model_time += dt
        gravhydro.evolve_model(model_time)
        channel["to_galaxies"].copy()
        channel["to_gas"].copy()
        #channel["to_moon"].copy()
        
        
    gravity.stop()
    hydro.stop()

t_end = 20.0 | units.Myr
gravity_hydro_bridge(gravity1, sph_code, gravhydro, 
                     galaxies, t_end)


# In[ ]:


plt.scatter(sph_code.gas_particles.x.value_in(units.kpc),
            sph_code.gas_particles.y.value_in(units.kpc),
            c = 'r')
plt.scatter(galaxies.x.value_in(units.kpc),
            galaxies.y.value_in(units.kpc),
            c = 'g')


# In[ ]:




