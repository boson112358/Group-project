from modules import *

import random
import math

from amuse.lab import Particles, units, Fi

#Here we define the track particles with some values which can change their location and radius of distribution

#If we want to use different position we have to make sure to properly change the velocity.
#I also assumed the velocity is roughly constant throughout the neighbourradius which only for small r is the case

def make_solar_system(N, solar_position, system_radius, galaxy_velocity, solar_tang_velocity):
    """
    here we generate our N bodies, we can change code in here to reflect different things
    We can distribut them randomly, have them near eachother, all at the same radius etc.
    All particles are in the x-y plane, this should still be adjusted accordingly to give some z randomization
    """
    
    particles = Particles(N)
    
    for body in particles:
        #For now I will keep the mass the same and distribute uniformly over a circle
        smallangle = random.random() * 2 * math.pi   #smallangle is angle around solarpos, angle is angle around Sag A
        radius = system_radius * math.sqrt(random.random()) #sqrt so it's uniform and not clumped at r=0
        body.mass = 5 | units.MSun 
        body.radius = 1.5 | units.RSun   
        #I have no z displacement, I can still add this though
        body.position = (radius * math.cos(smallangle) + solar_position[0], 
                         radius  *math.sin(smallangle) + solar_position[1], 
                         solar_position[2]) 
        #arctan only -pi/2 to pi/2, this codes gives us from -pi/2 to 3pi/2, we need global angle to fix vel vector
        if body.x < 0 | units.kpc:
            angle = math.atan(body.y / body.x) + math.pi
        else:
            angle = math.atan(body.y / body.x)
        
        #print(360*angle/(2*math.pi)) you can check that the angles are correct this way
        
        body.velocity = [-math.sin(angle) * solar_tang_velocity + galaxy_velocity[0], 
                         math.cos(angle) * solar_tang_velocity + galaxy_velocity[1], 
                         galaxy_velocity[2]] 
    return particles


def leapfrog_alg(current_iteration, particles, galaxy_solver):
    """
    here we define a leapfrog integration method for the tracker particles 
    through the get_gravity_at_point() method
    """
    
    unity = 1 | units.m**-1 * units.s**2 #We do this to remove the units inside the vector
    unity2 = 1 | units.m * units.s**-2  #Here we add the units on the outside
    
    for body in particles:     #We have to repeat this for all particles in our set of trackers
        accel_out = [0,0,0] 
        #We get the acceleration at x_i
        accel = galaxy_solver.get_gravity_at_point(3 | units.km, 
                                                   body.position[0],
                                                   body.position[1],
                                                   body.position[2])
        
        for i in range(len(accel_out)):
            accel_out[i] = (accel[i] * unity)  #removing units

        accel_out = accel_out * unity2 #units outside vector

        #We calculate v_{i+1/2} over the x y and z (I tried adding directly but got mismatched units)
        #and we calculate x_{i+1} with the velocity v_{i+1/2}

        if current_iteration==0:       #First timestep is only half for leapfrog
            tstep_factor = 0.5
        else:
            tstep_factor = 1

        body.velocity = body.velocity + accel_out * galaxy_solver.parameters.timestep * tstep_factor
        body.position = body.position + body.velocity * galaxy_solver.parameters.timestep * tstep_factor 