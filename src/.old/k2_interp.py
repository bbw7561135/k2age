#!/usr/bin/python
"""
  Mass-metallicity interpolation routine for mass tracks presented in Feiden & 
  Dotter (2012). 

  author: Gregory A. Feiden :: <gregory.a.feiden.gr@dartmouth.edu>
    date: 28 July 2012
 
  update:

  dependencies: SciPy (interpolation routines)
                NumPy (file handling)

"""

__all__ = ['__name__', 'find_age', 'mass_feh_interp']
        
def mass_feh_interp(mass, feh):
    """ mass and metallicity interpolation routine """
    import os, sys, k2_subs

    print '\n\t Interpolating for: \n'
    print '\t\t M/Mo   = ', mass
    print '\t\t [Fe/H] = ', feh

    # confirm mass and feh are within acceptable bounds
    if (k2_subs.value_check(mass, feh) == 1):
        print '\n\t Exiting Program \n'
        sys.exit(1)
    else:
        pass

    # determine directory structure
    src = os.path.dirname(k2_subs.__file__)
    trk = src[:(len(src)-4)] + '/trk'
    out = trk + '/usr_int'
    os.chdir(src)

    # locate file keys for interpolation
    m_files = k2_subs.find_mass(mass)
    f_files = k2_subs.find_feh(feh)

    # interpolate in metallicity
    k2_subs.feh_interp(m_files, f_files, feh, trk)
    
    # setup file name
    filename = [('./tmp/' + x + '_feh' + str(feh) + '.tmp') for x in m_files]
    
    # read in track files to single array
    # structure of array: data[-file][-row][-column]
    data = k2_subs.read_mass_tracks(filename)
    
    # mass interpolation
    track_final = k2_subs.mass_interp(data, m_files, mass, feh)
    
    # final file name creation
    if (feh < 0):
        fid = 'm' + str(int(abs(feh)*100))
    else:
        fid = 'p' + str(int(feh*100))
        
    ff_name = '/m' + str(int(mass*1000)) + '_feh' + fid + '_k2.trk'
    
    # clean up
    os.rename(track_final, out+ff_name)
    os.remove(filename[0])
    os.remove(filename[1])
    
    print '\n \t Interpolation complete! \n'
    print '\t Output directed to', ff_name
    print
    
    return ff_name

#-------------------------------------------------------------------------------#

def find_age(k2, fname):
    from numpy import genfromtxt as gft
    from scipy.interpolate import interp1d
    
    data = genfromtxt(fname, comments='#')
    icurve = interp1d(data[:,1],data[:,0],kind='linear')
    age = icurve(k2)
    
    return age

#-------------------------------------------------------------------------------#

if __name__ == "__main__":
    import sys
    # ensure correct number of command line arguments are issued
    if (len(sys.argv) != 3):
        print "Whoops: Incorrect number of args passed to the script. \n"
        print "\t Proper call is: python k2_interp.py *mass* *feh* \n"
        print "\t Exiting Program. \n"
        sys.exit(2)
    else:
        mass_feh_interp(float(sys.argv[1]),float(sys.argv[2]))
