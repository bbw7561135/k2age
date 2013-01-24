#!/usr/bin/python
#
import sys
from k2age import *

def defaultBinary(mass_a=None, mass_b=None, radius_a=None, radius_b=None, \
                  metallicity=None, eccentricity=None, semi_major_axis=None):
    
    if None in [mass_a, mass_b, radius_a, radius_b, metallicity, \
                eccentricity, semi_major_axis]:
        sys.exit('ERROR: Not all input variables were defined.')
    
    # Create a model object
    models = DsepModel()
    
    # Create a primary and secondary star
    primary_star = Star(mass=mass_a, radius=radius_a, metallicity=metallicity)
    secondary_star = Star(mass=mass_b, radius=radius_b, metallicity=metallicity)
    
    # Interpolate model grid to get individual mass tracks
    primary_track = primary_star.massTrack(model_set=models)
    secondary_track = secondary_star.massTrack(model_set=models)

    # Create an instance of a binary system from the stars
    binary = Binary(primary=primary_star, secondary=secondary_star, \
                    eccentricity=eccentricity, semi_major_axis=semi_major_axis)
    
    # Convolve mass tracks to derive weighted k2
    binary_track = binary.convolveTracks(primary_track=primary_track, \
                                         secondary_track=secondary_track)
                                         
    # Print new track to a file
    binary.printToFile(ages=models.log_ages, k2_track=binary_track)
    

if __name__ == '__main__':
    """ Command line call mirrors a call to defaultBinary()
    
    Calling this routine from the command line will default to calling 
    the specified routine above, if the command line arguments were 
    specified correctly. In total, 7 command line arguments are necessary
    and must be defined in the correct order. That is:
    
    python defaultBinary *mass_a* *mass_b* *radius_a* *radius_b* *metallicity*
                         *eccentricity* *semi_major_axis*
    
    Input
    -----
    mass_a           ::  Mass of the primary star (solar units)
    mass_b           ::  Mass of the secondary star (solar units)
    radius_a         ::  Radius of the primary star (solar units)
    radius_b         ::  Radius of the secondary star (solar units)
    metallicity      ::  Metallicity of the system (dex, relative to solar)
    eccentricity     ::  Eccentricity of the binary orbit
    semi_major_axis  ::  Semi-major axis of the binary orbit (solar radii)
    
    """
    if len(sys.argv) != 8:
        sys.exit('ERROR: Not all input variables were specified.')
    
    sys.argv.pop(0)
    arg = [float(x) for x in sys.argv]
    
    defaultBinary(mass_a=arg[0], mass_b=arg[1], radius_a=arg[2], radius_b=arg[3], \
                  metallicity=arg[4], eccentricity=arg[5], semi_major_axis=arg[6])
