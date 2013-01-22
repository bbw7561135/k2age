## @package k2age
#
import logging

logger = logging.getLogger('star')

__all__ = [Star]

class Star(object):
    
    ## Constructor for mass track loader
    #
    #  @param mass Mass of the star.
    #  @param radius Observed radius of the star.
    #  @param metallicity Metallicity of the star.
    #  @param angular_velocity Angular velocity of the star (optional).
    def __init__(self, mass=None, radius=None, metallicity=None, angular_velocity=None):
        
        if None in [mass, radius, metallicity]:
            logger.error('Must declare the properties of the star.')
        
        # stellar properties
        self.mass = mass
        self.radius = radius
        self.average_density = (mass*1.989e33)/(radius*6.956e10)**3
        self.metallicity = metallicity
        self.angular_velocity = angular_velocity
    
    ## Load the appropriate mass track for the star.
    #
    #  @param self Star object.
    #  @param model_set Stellar evolution model object.
    def massTrack(self, model_set=None):
        
        if None in [model_set]:
            logger.error('No stellar evolution models were selected')
        
        mass_track = model_set.getMassTrack(mass=self.mass, metallicity=self.metallicity)
        
        return mass_track
            
