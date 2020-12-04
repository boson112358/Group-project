#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[ ]:





# In[1]:


import numpy
import random,math
from amuse.plot import plot,scatter
from matplotlib import pyplot
from amuse.units import nbody_system,units
from amuse.lab import Particles
from amuse.community.fi.interface import Fi
from amuse.community.ph4.interface import ph4


# In[2]:
'''

N= 40   #These values and code are just for my test galaxy, will be replaced by the Fi solver
r = 8
v=230 

def MW_and_M31():
    MW_M31 = Particles(2+2*N)
    MW = MW_M31[0]
    MW.mass = 1*10**11 | units.MSun
    MW.position = (0,0,0) | units.kpc
    MW.velocity = (10,10,0) | units.kms
    #MW.u = 0.5*(math.sqrt(2)*10)**2*1000*1000 | units.m**2 * units.s**-2
    for i in range(N):
        MW_M31[1+i].mass = 10 | units.MSun
        angle = random.random()*2*math.pi
        MW_M31[1+i].position = (r*math.cos(angle),r*math.sin(angle),0) | units.kpc
        MW_M31[1+i].velocity = (10+v*math.sin(angle),10 - math.cos(angle)*v,0) | units.kms
        #Tried to add energy thing since Fi wanted it but didn't work
        #MW_M31[1+i].u = 0.5*(math.sqrt((10 + v*math.sin(angle))**2+10 - math.cos(angle)*v))**2*1000*1000 | units.m**2 * units.s**-2
    M31 = MW_M31[N+1]
    M31.mass = 1.6*10**11 | units.MSun
    M31.position = (780,0,0) | units.kpc
    M31.velocity = (0,0,0) | units.kms
    for i in range(N):
        MW_M31[N+2+i].mass = 10 | units.MSun
        angle = random.random()*2*math.pi
        MW_M31[N+2+i].position = (780 + r*math.cos(angle),r*math.sin(angle),0) | units.kpc
        MW_M31[N+2+i].velocity = (v*math.sin(angle),-math.cos(angle)*v,0) | units.kms
        #MW_M31[N+2+i].u = 0.5*230**2*1000*1000 | units.m**2 * units.s**-2
    return MW_M31


Galaxies=MW_and_M31()
N1 = 20
def MW():
    MW = Particles(1+N1)
    MW1 = MW[0]
    MW1.mass = 1*10**11 | units.MSun
    MW1.position = (0,0,0) | units.kpc
    MW1.velocity = (10,10,0) | units.kms
    #MW.u = 0.5*(math.sqrt(2)*10)**2*1000*1000 | units.m**2 * units.s**-2
    for i in range(N1):
        MW[1+i].mass = 10 | units.MSun
        angle = random.random()*2*math.pi
        MW[1+i].position = (r*math.cos(angle),r*math.sin(angle),0) | units.kpc
        MW[1+i].velocity = (10+v*math.sin(angle),10 - math.cos(angle)*v,0) | units.kms
    return MW
        
MW = MW()
mwc = MW.center_of_mass()


# In[3]:


#print(MW.center_of_mass())


# In[4]:


#Here we define the track particles with some values which can change their location and radius of distribution


Solarposition = (-6,6,0) | units.kpc  #If we displace MW, we can do it through this vector aswell
neighbourradius = 0.100 | units.kpc   #neighbourhood in which we distribute our stars
Ntracker=100                          #How many particles we will add
VelocityMilky = (10,10,0) | (units.km/units.s)  #Velocity at which the milkyway is traveling
vel = 220 | (units.km/units.s)   #This is roughly the velocity for our solarsystem around MW
#If we want to use different position we have to make sure to properly change the velocity.
#I also assumed the velocity is roughly constant throughout the neighbourradius which only for small r is the case
def Nstars():  
    particles = Particles(Ntracker)
    i=0
    while i < Ntracker:
        #here we generate our N bodies, we can change code in here to reflect different things
        #We can distribut them randomly, have them near eachother, all at the same radius etc.
        #All particles are in the x-y plane, this should still be adjusted accordingly to give some z randomization
        
        #For now I will keep the mass the same and distribute uniformly over a circle
        smallangle = random.random()*2*math.pi   #smallangle  is angle around solarpos, angle is angle around Sag A
        radius = neighbourradius * math.sqrt(random.random()) #sqrt so it's uniform and not clumped at r=0
        particles[i].mass = 5 | units.MSun 
        particles[i].radius = 1.5 | units.RSun   
        #I have no z displacement, I can still add this though
        particles[i].position = (radius*math.cos(smallangle) + Solarposition[0], radius*math.sin(smallangle) + Solarposition[1],Solarposition[2]) 
        #arctan only -pi/2 to pi/2, this codes gives us from -pi/2 to 3pi/2, we need global angle to fix vel vector
        if particles[i].x < 0 | units.kpc:
            angle = math.atan(particles[i].y / particles[i].x) + math.pi
        else:
            angle = math.atan(particles[i].y/particles[i].x)
        
        #print(360*angle/(2*math.pi)) you can check that the angles are correct this way
        
        particles[i].velocity = [math.sin(angle)*vel + VelocityMilky[0],-math.cos(angle)*vel+VelocityMilky[1],VelocityMilky[2]] 
        i += 1
    return particles
Tracker = Nstars()
#print(Tracker)


# In[5]:


#This ph4 will be replaced by the Fi later
converter=nbody_system.nbody_to_si(1.0e12|units.MSun, 100|units.kpc)
gravity = ph4(converter)
gravity.particles.add_particles(Galaxies)
channel1 = gravity.particles.new_channel_to(Galaxies)
channel2 = Tracker.new_channel_to(Tracker)
DeltaT = 0.1   #This is the timestep we will use in Myr
#If we want to have a non constant DeltaT we would need to change some of the code since leap-frog is not stable
#enough for varying timesteps. Can be done by finishing the leap-frog and starting a new one

gravity.timestep = DeltaT|units.Myr
times = numpy.arange(0., 50, DeltaT) | units.Myr


unity = 1 | units.m**-1 * units.s**2 #We do this to remove the units inside the vector
unity2 = 1 | units.m * units.s**-2  #Here we add the units on the outside
Accelout = [0,0,0] 
#Now that we have initialised our particle set we have the main gravitational solver with a certain timestep DeltaT
#This is what will later be the main solving part of Fi.
k=0
for time in times:
    #This is just the gravity solving, this will be replaced by the Fi solver in the main code
    gravity.evolve_model(time)
    channel1.copy()
    
    #print(Galaxies.position[3]) #I added this to quickly see if it works, I'll add some actual plots later
    #here we define a leapfrog integration method for the tracker particles through the get_gravity_at_point()
    if k ==0:       #First timestep is only half for leapfrog
        for i in range(Ntracker):  #We have to repeat this for all particles in our set of trackers
            Accelout = [0,0,0] 
            #We get the acceleration at x_i
            accel = gravity.get_gravity_at_point(3 | units.km, Tracker[i].position[0],Tracker[i].position[1],Tracker[i].position[2])
            for f in range(3):
                
                Accelout[f] = (accel[f]*unity)  #removing units
                
            Accelout = Accelout*unity2 #units outside vector
            
            #We calculate v_{i+1/2} over the x y and z (I tried adding directly but got mismatched units)
            #and we calculate x_{i+1} with the velocity v_{i+1/2}
            
            
                
            Tracker[i].velocity = Tracker[i].velocity + Accelout*gravity.timestep*0.5
                
            Tracker[i].position = Tracker[i].position + Tracker[i].velocity*gravity.timestep
              
            
            
    else:    #rest of the leapfrog, I might still have to change something for the last step of the frog
        
        for i in range(Ntracker):
            #We get the acceleration at x_i
            Accelout = [0,0,0]
            accel = gravity.get_gravity_at_point(3 | units.km, Tracker[i].position[0],Tracker[i].position[1],Tracker[i].position[2])
            for f in range(3):
                
                Accelout[f] = (accel[f]*unity)
                
            Accelout = Accelout*unity2
            #We calculate v_{i+1/2} over the x y and z (I tried adding directly but got mismatched units)
            #and we calculate x_{i+1} with the velocity v_{i+1/2}
            Tracker[i].velocity = Tracker[i].velocity + Accelout*gravity.timestep
              
            Tracker[i].position = Tracker[i].position + Tracker[i].velocity*gravity.timestep
            
        
    k += 1
    
    
    
    
    
gravity.stop()


# In[6]:


#for i in range(100):
    #scatter(Tracker[i].position[0],Tracker[i].position[1])
#for i in range(41):
    #scatter(Galaxies[i].x,Galaxies[i].y)
#Quick plot to see end result, need to make some actual data analysis functions and such which will come later

#The other cell's are mostly just me randomly testing stuff to make sure things work, they can be ignored


# In[7]:


#print(Tracker.position)
'''

# In[12]:

import os

folder = os.getcwd()+ '/Histograms'

if not os.path.exists(folder):
    os.makedirs(folder)



a = 8
b = 9
def solardistance(centermass, trakerpos,a,b):
    pyplot.figure()
    unity = 1| units.kpc**-2
    r = []
    for i in range(len(trakerpos)):
        d = centermass - trakerpos[i]
        r.append((math.sqrt((d[0]**2+d[1]**2+d[2]**2)*unity)))
    pyplot.hist(r,bins = 100,histtype='step',range=(a,b))
    pyplot.xlabel('r')
    pyplot.ylabel('number')
    #savepath = '/Histograms/'
    pyplot.savefig('./Histograms/histogram_0.png')
    #print(r)

#solardistance(mwc,Tracker.position,a,b)

    


# In[ ]:

