## Example of how to use this software to calculate an age from k2.
#
#  Example of how to use k2age to extract the age a detachted eclipsing
#  binary from the evolution of the interior structure constants, k_2.
#  The particular example below reproduces the data shown in Figure 3 of
#  Feiden & Dotter (2013). 
#
from k2age import *

def exampleBinary():

    # Create instance of theoretical model set
    models = DsepModel()
    
    # Create individual stars
    metallicity = 0.0
    primary = Star(mass=0.55, radius=0.62, metallicity=metallicity)
    secondary = Star(mass=0.25, radius=0.41, metallicity=metallicity)
    
    # Interpolate model grid to get individual mass tracks
    primary_track = primary.massTrack(model_set=models)
    secondary_track = secondary.massTrack(model_set=models)

    # Create an instance of a binary system from the stars
    binary = Binary(primary=primary, secondary=secondary, eccentricity=0.2, \
                    semi_major_axis=3.0)
    
    # Convolve mass tracks to derive weighted k2
    binary_track = binary.convolveTracks(primary_track=primary_track, \
                                         secondary_track=secondary_track)
                                         
    # Print new track to a file
    binary.printToFile(ages=models.log_ages, k2_track=binary_track)
    

