###### importing modules ######

import numpy as np
import matplotlib.pyplot as plt
from amuse.lab import units, Particles, nbody_system, Fi

from modules.common import __SCRIPT_PATH__, __MERGER_DIR__


###### setting up igm ######

def setup_sph_code(N1, N2, L, rho, u):
    #rho -- local group density
    #converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.kpc)
    #sph_code = sph_code(converter)#, mode = 'periodic')#, redirection = 'none')#change periodic
    #sph_code.parameters.periodic_box_size = 10.0 | units.kpc
    dm = Particles(N1)
    dm.mass = 0.8*(rho * L**3) / N1
    #np.random.seed(12345)
    dm.x = L * np.random.uniform(-0.5, 0.5, N1)
    dm.y = L * np.random.uniform(-0.5, 0.5, N1)
    dm.z = L * np.random.uniform(-0.5, 0.5, N1)
    dm.vx = np.zeros(N1) | units.km / units.s
    dm.vy = np.zeros(N1) | units.km / units.s
    dm.vz = np.zeros(N1) | units.km / units.s
    gas = Particles(N2)
    gas.mass = 0.2*(rho * L**3) / N2
    gas.x = L * np.random.uniform(-0.5, 0.5, N2)
    gas.y = L * np.random.uniform(-0.5, 0.5, N2)
    gas.z = L * np.random.uniform(-0.5, 0.5, N2)
    gas.vx = np.zeros(N2) | units.km / units.s
    gas.vy = np.zeros(N2) | units.km / units.s
    gas.vz = np.zeros(N2) | units.km / units.s
    gas.u = u
    gas.h_smooth = L / N2**(1/3.0)
    dm.h_smooth = L / N1**(1/3.0)
        
    #sph_code.gas_particles.add_particles(gas)
    #sph_code.dm_particles.add_particles(dm)
    #sph_code.commit_particles()
    
    return gas, dm


###### setting up grid ######

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


'''
N1 = 10000
N2 = 130000
Ngrid = 500
L = 1000 | units.kpc
Lg = 1000
rho = 770 | units.MSun / (units.kpc)**3
u = 3.724e+9 | (units.m)**2 / (units.s)**2
sph_code = setup_sph_code(Fi, N1, N2, L, rho, u)
'''

###### igm plot function ######

def plot_igm(sph_code, Ngrid, Lg, output_dir, filename):
    x, y, z, vx, vy, vz = setup_grid(Ngrid,Lg)
    rho, rhovx, rhovy, rhovz, rhoe = sph_code.get_hydro_state_at_point(x, y, z, vx, vy, vz)
    rho = rhoe.reshape((Ngrid + 1, Ngrid + 1))
    max_dens = np.log10(rho.value_in(units.m**-1 * units.s**-2 * units.kg)).max()
    min_dens = np.log10(rho.value_in(units.m**-1 * units.s**-2 * units.kg)).min()
    #print("Done running")
    # plt.scatter(sph_code.gas_particles.x.value_in(units.kpc),
    #             sph_code.gas_particles.y.value_in(units.kpc),
    #             c = 'r')
    # plt.scatter(sph_code.dm_particles.x.value_in(units.kpc),
    #             sph_code.dm_particles.y.value_in(units.kpc),
    #             c = 'b')
    cax = plt.imshow(np.log10(rho.value_in(units.m**-1 * units.s**-2 * units.kg)),
                     extent=[-Lg, Lg, -Lg, Lg], vmin=min_dens, vmax=max_dens,
                     origin='lower', cmap="hot")
    plt.colorbar(cax, ticks=[1.e-8, 0.5*max_dens,max_dens], orientation='vertical', fraction=0.045)
    plt.savefig(output_dir + filename)
    #print(np.log10(rho.value_in(units.m**-1 * units.s**-2 * units.kg)).max()
    #      ,np.log10(rho.value_in(units.m**-1 * units.s**-2 * units.kg)).min())
    #print("done plotting")