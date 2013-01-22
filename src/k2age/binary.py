## @package k2age
#
import logging
import numpy as np
import scipy.interpolate as sint

logger = logging.getLogger('binary')

__all__ = [Binary]

class Binary(object):
    
    ## Constructor for the binary system
    #
    #  @param self Binary object reference.
    #  @param primary Star object for primary.
    #  @param secondary Star object for secondary.
    #  @param eccentricity Eccentricity of the binary orbit.
    #  @param semi_major_axis Semi-major axis of the binary orbit (Rsun).
    #  @param v_anuglar_orbit Mean orbital angular velocity (optional).
    def __init__(self, primary=None, secondary=None, eccentricity=None, \
                 semi_major_axis=None, v_angular_orbit=None):
                 
        if None in [primary, secondary, eccentricity, semi_major_axis]:
            logger.error('Not all of the binary properties were specified.')
            
        if primary.metallicity != secondary.metallicity:
            logger.error('Your stars have different metallicities!')
        
        # Set defined binary properties.
        self.primary_mass = primary.mass
        self.secondary_mass = secondary.mass
        self.primary_radius = primary.radius
        self.secondary_radius = secondary.radius
        self.metallicity = primary.metallicity
        self.eccentricity = eccentricity
        self.semi_major_axis = semi_major_axis
        self.v_angular_A = primary.angular_velocity
        self.v_angular_B = secondary.angular_velocity
        self.v_angular_orbit = v_angular_orbit
        
        # Compute c_2,i coefficients.
        self.c2 = self.getC2Coefficients()
        
        
    ## Compute \f$c_{2, i}\f$ coefficients.
    #
    #  These coefficients weight the individual k2 values to derive a 
    #  weighted average of the two binary components. They are computed
    #  according to Equation (2) in Feiden & Dotter (2013),
    #  \f$c_{2,\, i} = \left[\left(\frac{\Omega_i}{\Omega_K}\right)^2 
    #  \left(1 + \frac{m_{3-i}}{m_i}\right) f(e) + \frac{15 m_{3-i}}{m_i}
    #   g(e) \right]\left(\frac{R_i}{A}\right)^5\f$.
    #
    #  @param self Binary system object.
    def getC2Coefficients(self):
        
        # Compute functions f(e) and g(e) 
        # Equations (3) and (4) in Feiden & Dotter (2013)
        f = (1. - self.eccentricity**2)**(-2)
        g = (8. + 12.*self.eccentricity**2 + self.eccentricity**4)*f**(5./2.)/8.
        
        # Determine ratio of rotational and orbital angular velocity.
        if None in [self.v_angular_A, self.v_angular_orbit]:
            # assume pseudo-synchronization
            omega_ratio_1 = (1. + self.eccentricity)/(1. - self.eccentricity)**3
        else:
            omega_ratio_1 = (v_angular_A/v_angular_orbit)**2
        
        if None in [self.v_angular_B, self.v_angular_orbit]:
            # assume pseudo-synchronization
            omega_ratio_2 = (1. + self.eccentricity)/(1. - self.eccentricity)**3
        else:
            omega_ratio_2 = (v_angular_B/v_angular_orbit)**2
        
        omega_ratio = [omega_ratio_1, omega_ratio_2]
        
        # Compute c_2,1
        c21 = (omega_ratio[0]*(1. + self.secondary_mass/self.primary_mass)*f + \
              15.*self.secondary_mass*g/self.primary_mass) * \
              (self.primary_radius/self.semi_major_axis)**5
        # Compute c_2,2
        c22 = (omega_ratio[1]*(1. + self.primary_mass/self.secondary_mass)*f + \
              15.*g*self.primary_mass/self.secondary_mass) * \
              (self.secondary_radius/self.semi_major_axis)**5
              
        return c21, c22
    
    ## Convolve two model mass tracks to get the k2-age relation.
    #
    #  The convolution of the two mass tracks takes the weighted average
    #  of the individual \f$k_{2,\, i}\f$ values. This is presented in 
    #  Equation (5) of Feiden & Dotter (2013).
    #  \f$\overline{k_{2}} = \frac{c_{2,\,1}k_{2,\,1} + 
    #      c_{2,\,2}k_{2,\,2}}{c_{2,\,1} + c_{2,\,2}}\f$.
    #
    #  @param self Binary object.
    #  @param primary_track Mass track of the primary star.
    #  @param secondary_track Mass track of the secondary star.
    def convolveTracks(self, primary_track=None, secondary_track=None):
        
        if None in [primary_track, secondary_track]:
            logger.error('Not all mass tracks were defined')
        
        track = [((self.c2[0]*primary_track[i] + self.c2[1]*secondary_track[i]) /
                 (self.c2[0] + self.c2[1])) for i in range(len(primary_track))]
                         
        return track
    
    ## Print binary track out to a file
    #
    #  @param self Binary object.
    #  @param ages Set of ages for convolved binary track.
    #  @param k2_track List of k2 values to be output.
    def printToFile(self, ages=None, k2_track=None):
        
        filename = './binary_star_track.out'
        fout = open(filename, 'w')
        
        for i in range(len(ages)):
            s = '{:4.2f} \t {:10.8f} \n'.format(ages[i], k2_track[i])
            fout.write(s)
        fout.close()
    
    ## Interpolate along binary track to find the age using \f$\overline{k_{2}}\f$.
    #
    #  @param self Binary object.
    #  @param binary_track Convolved binary track.
    #  @param age_list List of ages for convolved binary track.
    #  @param k2 Observed weighted apsidal motion constant.
    def k2Age(self, binary_track=None, age_list=None, k2=None):
        
        if None in [binary_track, age_list, k2]:
            logger.error('Not all input parameters were declared.')
        
        interp_curve = sint.interp1d(binary_track, age_list , kind='linear')
        
        return interp_curve(k2)
        
