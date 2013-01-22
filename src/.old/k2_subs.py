
__all__ = ['equal_grid', 'feh_interp', 'feh_val_key', 'find_feh', 
           'find_mass', 'mass_interp', 'mass_val_key', 'read_mass_tracks',
           'value_check']

# Dictionaries required for file identification
#----------------------------------------------------------------------#

feh_id = {-1.0: 'm100', -0.5: 'm50', -0.3: 'm30', -0.1: 'm10', 
           0.0: 'p0', 0.1: 'p10',  0.2: 'p20', 0.3: 'p30'}

mass_id = {0.15: 'm0150', 0.17: 'm0170', 0.19: 'm0190', 0.21: 'm0210',
           0.23: 'm0230', 0.25: 'm0250', 0.27: 'm0270', 0.29: 'm0290', 
           0.31: 'm0310', 0.33: 'm0330', 0.35: 'm0350', 0.37: 'm0370',
           0.39: 'm0390', 0.41: 'm0410', 0.43: 'm0430', 0.45: 'm0450', 
           0.47: 'm0470', 0.49: 'm0490', 0.51: 'm0510', 0.53: 'm0530', 
           0.55: 'm0550', 0.57: 'm0570', 0.59: 'm0590', 0.61: 'm0610', 
           0.63: 'm0630', 0.65: 'm0650', 0.67: 'm0670', 0.69: 'm0690', 
           0.71: 'm0710', 0.73: 'm0730', 0.75: 'm0750'}


# Subroutines 
#----------------------------------------------------------------------#

def value_check(mass, feh):
    """ Error handling routine to ensure Mass and [Fe/H] are within the
        bounds required by the program """
    
    # check [Fe/H] is within file set
    f_list = feh_id.keys()
    f_list.sort()
    
    if (f_list[0] <= feh <= f_list[len(f_list) - 1]):
        err = 0
    else:
        print '\n'
        print 'Whoops: [Fe/H] value outside of acceptable bounds.\n'
        print '\t Must have: -1.0 <= [Fe/H] <= +0.3 \n \n'
        err = 1
        
    # check mass is within file set
    m_list = mass_id.keys()
    m_list.sort()

    if (m_list[0] <= mass <= m_list[len(m_list) - 1]):
        err = 0
    else:
        print 'Whoops: Mass outside of allowed bounds. \n'
        print '\t Allowed masses: 0.15 <= M/Mo <= 0.75  \n\n'
        err = 1
        
    return err
    
#----------------------------------------------------------------------#

def find_mass(mass):
    """ find nearest files appropriate for mass interpolation """
        
    # transfer dictionary keys to list
    m_list = mass_id.keys()
    
    m_list.insert(0, mass)  # insert test mass
    m_list.sort()           # re-sort mass list
    m = m_list.index(mass)  # locate index of test mass
        
    m_files = [mass_id[m_list[m - 1]], mass_id[m_list[m + 1]]]
        
    return m_files
    
#----------------------------------------------------------------------#

def find_feh(feh):
    """ locate metallicit(y/ies) to be used in interpolation """
    
    # transfer dictionary keys to list
    f_list = feh_id.keys()
    
    # locate surrounding metallicities
    f_list.insert(0, feh)  # insert test [Fe/H]
    f_list.sort()          # re-sort [Fe/H] list
    m = f_list.index(feh)  # locate index of test [Fe/H]
        
    f_files = [feh_id[f_list[m-1]], feh_id[f_list[m+1]]]
        
    return f_files
    
#----------------------------------------------------------------------#

def feh_interp(m_files, f_files, feh, trk):
    """ interpolate each set of masses in metallicity """
    
    for i in m_files:
        # equally space time grid
        age = [(1.e6 + h*1.e4) for h in range(int(9.9e5))]
        feh_data = [equal_grid(i, j, trk) for j in f_files]
        
        # convert dictionary values back to real values
        f_val = feh_val_key(f_files)
        
        # interpolate in metallicity
        m = (feh_data[1] - feh_data[0])/(f_val[1] - f_val[0])
        par = [(feh_data[0][j] + m[j]*(feh - f_val[0])) for j in range(len(age))]
        
        # create tmp file containing feh interpolated mass track
        file_new = './tmp/' + i + '_feh' + str(feh) + '.tmp'
        f_out = open(file_new, 'w')
        
        # output new list of k2 vs age
        for j in range(len(age)):
            s = str(age[j]) + '\t' + str(par[j]) + '\n'
            f_out.write(s)
        f_out.close()
    
    return

#----------------------------------------------------------------------#

def read_mass_tracks(filename):
    """ reads in mass tracks for masses in m_files """
    from numpy import genfromtxt
    
    # read in data to single array
    # structure of array: data[-file][-row][-column]
    data = [(genfromtxt(x, comments='#')) for x in filename]
    
    return data

#----------------------------------------------------------------------#

def mass_interp(data, m_files, mass, feh):
    """ interpolates DSEP mass tracks to yield a final mass track """
    
    # equally space time grid
    age = [(1.e6 + h*1.e4) for h in range(int(9.9e5))]
    
    # convert dictionary value to actual mass
    m_val = mass_val_key(m_files)
    
    # interpolate in mass
    m = [(data[1][j][1] - data[0][j][1])/(m_val[1] - m_val[0]) \
           for j in range(len(age))]
    par = [(data[0][j][1] + m[j]*(mass - m_val[0])) for j in range(len(age))]

    # create new tmp file
    file_new = './tmp/m' + str(mass*1000) + '_feh' + str(feh) + '.tmp'
    f_out = open(file_new, 'w')
    
    # write file header
    s = '#  Age \t k2 \n# (yr)\n'
    f_out.write(s)
    
    # output new list of k2 vs age
    for j in range(len(age)):
        s = str(age[j]) + '\t' + str(par[j]) + '\n'
        f_out.write(s)
    f_out.close()
    
    return file_new

#----------------------------------------------------------------------#

def equal_grid(m_file, f_file, trk):
    """ place all mass tracks on equal time grid """
    from scipy.interpolate import interp1d
    from numpy import genfromtxt as gft
    
    # setup file name and read in data
    file_name = trk + '/feh' + f_file + '/' + m_file + '_GS98_' + \
                f_file + '_T60.iso'
    data = gft(file_name, comments='#', skiprows=1, usecols=(0,12))
    
    # generate interpolation curves
    icurve = interp1d(data[:,0], data[:,1], kind='linear')
    
    # place data on equally spaced grid
    age_grid = [(1.e6 + i*1.e4) for i in range(int(9.9e5))]
    feh_data = icurve(age_grid) 
    
    return feh_data

#----------------------------------------------------------------------#

def feh_val_key(dict_val):
    """ converts feh dictionary value back to a real quantity """
    
    dict_key = [float(dict_val[i][1:])/100. for i in range(len(dict_val))]
    
    for i in range(len(dict_val)):
        if (dict_val[i][0] == 'm'):
            dict_key[i] = -1.*dict_key[i]
        else:
            pass
    
    return dict_key

#----------------------------------------------------------------------#

def mass_val_key(dict_val):
    """ converts mass dictionary value back to real quantity """
    
    dict_key = [float(dict_val[i][1:])/1000. for i in range(len(dict_val))]
    
    return dict_key

#----------------------------------------------------------------------#
