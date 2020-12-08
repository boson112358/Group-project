from amuse.community.gadget2.interface import Gadget2 
from amuse.couple.bridge import CalculateFieldForParticles

class Gadget2Gravity(Gadget2): 
    """ Extension of Gadget2 with gravity included """ 
    def init(self, 
             converter, 
             #mode='normal', 
             **kwargs): 
        Gadget2.__init__(self, 
                         converter, 
                         #mode=mode, 
                         **kwargs) 
        
    def get_gravity_at_point(self, radius, x, y, z): 
        field_code = CalculateFieldForParticles(particles=self.dm_particles) 
        return field_code.get_gravity_at_point(radius, x, y, z) 
    
    def get_potential_at_point(self, radius, x, y, z): 
        field_code = CalculateFieldForParticles(particles=self.dm_particles) 
        return field_code.get_potential_at_point(radius, x, y, z)