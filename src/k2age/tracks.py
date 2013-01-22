## @package k2age
#
import logging
import numpy as np
import scipy.interpolate as sint

logger = logging.getLogger('tracks')

__all__ = [DsepModel]

class DsepModel(object):
    
    ## Constructor for mass track generator.
    #  
    #  @param self Model mass track object.
    def __init__(self):
        
        # map masses to string for file name construction
        self.mass_id = {0.15: 'm0150', 0.17: 'm0170', 0.19: 'm0190', 
                        0.21: 'm0210', 0.23: 'm0230', 0.25: 'm0250', 
                        0.27: 'm0270', 0.29: 'm0290', 0.31: 'm0310', 
                        0.33: 'm0330', 0.35: 'm0350', 0.37: 'm0370',
                        0.39: 'm0390', 0.41: 'm0410', 0.43: 'm0430', 
                        0.45: 'm0450', 0.47: 'm0470', 0.49: 'm0490', 
                        0.51: 'm0510', 0.53: 'm0530', 0.55: 'm0550', 
                        0.57: 'm0570', 0.59: 'm0590', 0.61: 'm0610', 
                        0.63: 'm0630', 0.65: 'm0650', 0.67: 'm0670', 
                        0.69: 'm0690', 0.71: 'm0710', 0.73: 'm0730', 
                        0.75: 'm0750'}
        
        # map metallicities to string for file name construction
        self.feh_id = {-1.0: 'm100', -0.5: 'm50', -0.3: 'm30', -0.1: 'm10', 
                        0.0: 'p0', 0.1: 'p10',  0.2: 'p20', 0.3: 'p30'}
        
        # generate list of available masses and metallicities
        self.mass_list = [round((0.15 + i*0.02), 2) for i in range(31)]
        self.feh_list = [-1.0, -0.5, -0.3, -0.1, 0.0, 0.1, 0.2, 0.3]
        
        # generate list of valid ages
        self.log_ages = [6.0 + i*0.05 for i in range(81)]
        self.ages = [10.**x for x in self.log_ages]
        
        # directory location of model tracks
        self.directory = '../trk/'
        
    ## Get model mass track by interpolating in model grid.
    #
    #  @param self Model set object.
    #  @param mass Mass of the requested track.
    #  @param metallicity Metallicity of the requested track.
    def getMassTrack(self, mass=None, metallicity=None):
        
        if None in [mass, metallicity]:
            logger.error('Need to specify stellar properties')
        
        # find location of mass and metallicity in the grid
        mass_idx = np.searchsorted(self.mass_list, mass, side='left')
        feh_idx = np.searchsorted(self.feh_list, metallicity, side='left')
        
        # interpolate in mass at each metallicity
        #--- low metallicity
        filename = self.getFileName(mass=self.mass_list[mass_idx],
                               metallicity=self.feh_list[feh_idx])
        mass_lo_feh_lo = self.loadMassTrack(filename=filename)
        
        filename = self.getFileName(mass=self.mass_list[mass_idx + 1],
                               metallicity=self.feh_list[feh_idx])
        mass_hi_feh_lo = self.loadMassTrack(filename=filename)
                                     
        mass_feh_lo = self.massTrackInterpolate(mass_track_1=mass_lo_feh_lo,
                                           mass_track_2=mass_hi_feh_lo,
                                           x_1=self.mass_list[mass_idx],
                                           x_2=self.mass_list[mass_idx+1],
                                           x_new=mass,
                                           set_equal_grid=True)
                                           
        #--- high metallicity                               
        filename = self.getFileName(mass=self.mass_list[mass_idx],
                               metallicity=self.feh_list[feh_idx + 1])
        mass_lo_feh_hi = self.loadMassTrack(filename=filename)
        
        filename = self.getFileName(mass=self.mass_list[mass_idx + 1],
                               metallicity=self.feh_list[feh_idx + 1])
        mass_hi_feh_hi = self.loadMassTrack(filename=filename)
                                     
        mass_feh_hi = self.massTrackInterpolate(mass_track_1=mass_lo_feh_hi,
                                           mass_track_2=mass_hi_feh_hi,
                                           x_1=self.mass_list[mass_idx],
                                           x_2=self.mass_list[mass_idx+1],
                                           x_new=mass,
                                           set_equal_grid=True)
                                           
        # interpolate in metallicity
        new_track = self.massTrackInterpolate(mass_track_1=mass_feh_lo, 
                                         mass_track_2=mass_feh_hi,
                                         x_1=self.feh_list[feh_idx],
                                         x_2=self.feh_list[feh_idx+1],
                                         x_new=metallicity,
                                         set_equal_grid=False)          
        return new_track
        
    ## Load requested mass track into an array.
    #  
    #  @param self Model track object.
    #  @param filename File name of the mass track to be loaded.
    def loadMassTrack(self, filename=None):
        
        if None in [filename]:
            logger.error('Please specify filename for mass track.')
        
        return np.genfromtxt(filename, comments='#', usecols=(0,1,2,3,4,12))
    
    ## Create filename for a star or for the two nearest grid points.
    #
    #  @param self Model set object.
    #  @param mass Mass of the model.
    #  @param metallicity Metallicity of the model.
    def getFileName(self, mass=None, metallicity=None):
    
        if None in [mass, metallicity]:
            logger.error('Mass track properties not specified.')
        
        directory = self.directory + 'feh' + self.feh_id[metallicity]
        filename = self.mass_id[mass] + '_GS98_' + self.feh_id[metallicity]
        filename =  directory + '/' + filename + '_T60.iso'
        
        return filename
        
    ## Set mass track onto an evenly spaced grid.
    #
    #  @param self Model track set object.
    #  @param mass_track Mass track to be set on an even grid.
    def setEqualGrid(self, mass_track=None):
        
        if None in [mass_track]:
            logger.error('Mass track not specified for equal grid transformation.')
    
        # interpolate k2 onto equal grid
        interpolation_curve = sint.interp1d(mass_track[:,0], mass_track[:,5], kind='linear')
        new_mass_track = interpolation_curve(self.ages)
        
        return new_mass_track
    
    ## Interpolate between two mass tracks.
    #
    #  @param self Model track set object.
    #  @param mass_track_1 Array containing the first mass track.
    #  @param mass_track_2 Array containing the second mass track.
    #  @param x_1 Value of the dependent variable for mass_track_1.
    #  @param x_2 Value of the dependent variable for mass_track_2.
    #  @param x_new Value of the new dependent variable.
    #  @param set_equal_grid Flag to set track on an equally spaced grid.
    def massTrackInterpolate(self, mass_track_1=None, mass_track_2=None, \
                              x_1=None, x_2=None, x_new=None, set_equal_grid=True):
        
        if None in [mass_track_1, mass_track_2, x_1, x_2, x_new]:
            logger.error('Mass track interpolation failed: No mass tracks declared.')
        
        # set mass tracks on an equal age grid
        if set_equal_grid:
            new_mass_track_1 = self.setEqualGrid(mass_track_1)
            new_mass_track_2 = self.setEqualGrid(mass_track_2)
        else:
            new_mass_track_1 = mass_track_1
            new_mass_track_2 = mass_track_2
        
        # confirm mass tracks have the same length
        if len(new_mass_track_1) != len(new_mass_track_2):
            logger.warning('Mass tracks have different lengths!')
            new_mass_track_1 = self.setEqualGrid(new_mass_track_1)
            new_mass_track_2 = self.setEqualGrid(new_mass_track_2)
        else:
            pass
        
        # calculate the slope of the mass track at each time step
        slopes = [(new_mass_track_1[i] - new_mass_track_2[i])/(x_1 - x_2) \
                  for i in range(len(new_mass_track_1))]
        
        # compute the new mass track
        new_track = [(new_mass_track_1[i] + slopes[i]*(x_new - x_1)) \
                     for i in range(len(slopes))]
                     
        return new_track
        
