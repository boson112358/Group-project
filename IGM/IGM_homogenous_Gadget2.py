#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from amuse.lab import *
#from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.units import units
from amuse.lab import Particles
from amuse.community.ph4.interface import ph4
from amuse.community.gadget2.interface import Gadget2


# In[2]:


def setup_sph_code(sph_code, N1, N2, L, rho,u):
    #rho -- local group density
    converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.kpc)
    sph_code = sph_code(converter,mode = 'periodic')#, redirection = 'none')
    sph_code.parameters.periodic_box_size = 10.0 | units.kpc
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
    gas.mass = 0.0001*(rho * L**3) / N2
    gas.x = L * np.random.uniform(0.0, 1.0, N2)
    gas.y = L * np.random.uniform(0.0, 1.0, N2)
    gas.z = L * np.random.uniform(0.0, 1.0, N2)
    gas.vx = np.zeros(N2) | units.km / units.s
    gas.vy = np.zeros(N2) | units.km / units.s
    gas.vz = np.zeros(N2) | units.km / units.s
    gas.u = u 
    if isinstance(sph_code, Fi):
        sph_code.parameters.self_gravity_flag = False
        sph_code.parameters.timestep = 0.1 | units.s
        gas.h_smooth = L / N2**(1/3.0)
        dm.h_smooth = L / N1**(1/3.0)
        #gas.position -= 0.5 * L
        
    sph_code.gas_particles.add_particles(gas)
    sph_code.dm_particles.add_particles(dm)
    sph_code.commit_particles()
    return sph_code


# In[3]:


N1 = 5000
N2 = 1000
L = 10 | units.kpc
rho = 1000 | units.MSun / (units.kpc)**3
u = 1.6e+15 | (units.m)**2 / (units.s)**2
sph_code = setup_sph_code(Gadget2, N1, N2, L, rho, u)


# In[4]:


import os

folder = os.getcwd()+ '/plot_IGM'

if not os.path.exists(folder):
    os.makedirs(folder)


# In[5]:

def evolve_sph():
    sph_code.evolve_model(1.0|units.s)
    print("Done running")
    plt.scatter(sph_code.gas_particles.x.value_in(units.kpc),
                sph_code.gas_particles.y.value_in(units.kpc),
                c = 'r')
    plt.scatter(sph_code.dm_particles.x.value_in(units.kpc),
                sph_code.dm_particles.y.value_in(units.kpc),
                c = 'b')
    plt.savefig('./plot_IGM/gas_dm_periodic.png')
    print("done plotting")


# In[ ]:




