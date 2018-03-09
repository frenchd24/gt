#!/usr/bin/env python


'''

Test diameter chooser

'''

from utilities import *

def choose_diameter_measurement(options):
    # choose which measurement of diameter to use from the list 'options'
    # returns the list [key,(major,minor),dRatio,pa] for the chosen measurement
    
    order = ['K_s (2MASS "total")',\
    'K_s (LGA/2MASS "total")',\
    'K_s (2MASS isophotal)',\
    'K_s (LGA/2MASS "isophotal")',\
    'POSS1 103a-O',\
    'POSS1 103a-E',\
    'ESO-LV "Quick Blue" IIa-O',\
    'r (SDSS Isophotal)',\
    'RC3 D_0 (blue)',\
    'RC3 D_25, R_25 (blue)',\
    'r (SDSS Petrosian)',\
    'r (SDSS deVaucouleurs)',\
    'r (SDSS de Vaucouleurs)',\
    'RC3 A_e (Johnson B)',\
    'r (SDSS Exponential)',\
    'B (Johnson)',\
    'R (Kron-Cousins)',\
    'ESO-Uppsala "Quick Blue" IIa-O',\
    'ESO-LV IIIa-F',\
    ]
    
    
    
    # first prepare the list. If there are multiple surveys under the same name,
    # keep only the best
    d = {}
    for n in options:
        key = n[0]
        
        # if you encounter a duplicate survey, take the best pieces
        if d.has_key(key):
            other = d[key]
            
            otherMajor,otherMinor = other[1]
            otherPA = other[3]
            
            major,minor = eval(n[1])
            pa = n[3]
            
            # best case, just choose the larger measurement
            if isNumber(major) and isNumber(minor) and isNumber(pa):
                if isNumber(otherMajor) and isNumber(otherMinor) and isNumber(otherPA):
                    if major != minor:
                        if float(major) > float(otherMajor):
                            # if this one is better than the last one, replace it
                            d[key] = n
                        elif float(otherMajor) == float(otherMinor):
                            # if the current major axis is not the larger, but is 
                            # different from the minor axis value whereas the other
                            # survey's is not, replace it
                            d[key] = n
                    
                    # if they're both equal, choose the bigger one
                    elif float(otherMajor) == float(otherMinor):
                        if float(major) > float(otherMajor):
                            d[key] = n
                
                # if not isNumber(*) for all the things, choose the new one
                else:
                    d[key] = n
            
            # if the current one sucks, see how much the old one sucks
            elif isNumber(major) and isNumber(otherMajor):
                # both have major
                if isNumber(minor) and not isNumber(otherMinor):
                    # OTHER does not have minor
                    if isNumber(pa) and not isNumber(otherPA):
                        # if the other one just sucks, choose this one
                        d[key] = n
                    else:
                        # if OTHER has a PA and THIS doesn't, but otherwise
                        # THIS is better, just grab the PA from OTHER
                        new_n = n
                        new_n[3] = otherPA
                        d[key] = new_n
                        
                # if THIS one has a PA, but not a minor, take the PA
                elif isNumber(otherMinor) and not isNumber(minor):
                    if isNumber(pa) and not isNumber(otherPA):
                        new_n = other
                        new_n[3] = pa
                        d[key] = new_n
                        
            elif isNumber(major) and not isNumber(otherMajor):
                # this one has a major, the other does not
                
                if not isNumber(pa) and isNumber(otherMajor):
                    d[key] = n
                
                
        else:
            # if d does not yet have this key yet
            d[key] = n
            
            
            
    
    found = False
    found_no_pa = False
    found_no_minor = False
    found_only_major = False
    for i in order:
        for n in options:
            # n is input string, i is what i'm looking for
            if bfind(n,i) and not found:
                major = n[1][0]
                minor = n[1][1]
                pa = n[3]

                if isNumber(major) and isNumber(minor) and isNumber(pa):
                    # skip the weird cases where one measurement option has major = minor
                    if float(major) != float(minor):
                        best = n
                        found = True
                        break
                    
                elif isNumber(minor) and not isNumber(pa) and not found_no_pa:
                    found_no_pa = n
                
                elif not isNumber(minor) and isNumber(pa) and not found_no_minor:
                    found_no_minor = n
            
                elif not found_only_major:
                    found_only_major = n

    if found:
        to_return = best
        keys = [to_return[0],to_return[0]]
        
    elif found_no_pa:
        if found_no_minor:
            # grab the pa from the next best survey
            pa = found_no_minor[3]
            pa_key = found_no_minor[0]
            
            # grab the rest from the one with intact (major, minor), but no pa
            to_return = found_no_pa
            
            # replace the missing pa
            to_return[3] = pa
            
            # return a key saying what came from where
            keys = [to_return[0],pa_key]
            
        
    elif found_no_minor:
        to_return = found_no_minor
        keys = [found_no_minor[0],found_no_minor[0]]
    
    else:
        to_return = options[0]
        keys = 'No best fit'

        sys.stdout.write('\n No fit was found for the following: {0}'.format(options))
        sys.stdout.write('\n Exiting....')
        sys.exit()

    return to_return, keys
    
    
def main():
    
#     options = [['r (SDSS Isophotal)', (5.35, 4.91), 0.92, 152.0], ['r (SDSS Petrosian)', (1.18, None), None, None], ['r (SDSS deVaucouleurs)', (0.01, 0.01), 1.0, 163.0], ['r (SDSS Exponential)', (0.0, 0.0), 0.81, 163.0]]
#     options = [['r (SDSS Exponential)', (0.0, 0.0), 0.81, 163.0], ['r (SDSS Petrosian)', (1.18, None), None, None], ['r (SDSS deVaucouleurs)', (0.01, 0.01), 1.0, 163.0]]

#     options = [['POSS1 103a-O', (1380.0, 480.0), 0.35, 162.0], ['POSS1 103a-O', (270.0, 33.0), 0.12, None], ['POSS1 103a-O', (1350.0, 372.0), 0.28, None], ['RC3 D_25, R_25 (blue)', (1312.7, 434.5), 0.33, 162.0], ['RC3 A_e (Johnson B)', (476.6, None), None, None], ['RC3 D_0 (blue)', (1312.7, None), None, None], ['K_s (LGA/2MASS isophotal)', (146.4, 146.4), 1.0, 16.0], ['K_s (LGA/2MASS "total")', (250.2, 250.2), 1.0, 16.0]]
    options = [('FUV (GALEX) AB', 20.6129, '+/-0.217785', 'mag', 1950000000000000.0, 2.0599999999999999e-05, '+/-4.14E-06', 'Jy', '2012GASC..C...0000S', 'uncertainty', '1538.6     A', 'Broad-band measurement', 'Flux integrated from map', 'Kron flux in elliptical aperture', 'From new raw data'), ('FUV (GALEX) AB', 21.918099999999999, '+/-0.369618', 'mag', 1950000000000000.0, 6.2099999999999998e-06, '+/-2.11E-06', 'Jy', '2012GASC..C...0000S', 'uncertainty', '1538.6     A', 'Broad-band measurement', 'Flux in fixed aperture', 'Flux in 7.5 arcsec diameter aperture', 'From new raw data'), ('NUV (GALEX) AB', 19.1554, '+/-0.0683089', 'mag', 1290000000000000.0, 7.8999999999999996e-05, '+/-4.97E-06', 'Jy', '2012GASC..C...0000S', 'uncertainty', '2315.7     A', 'Broad-band measurement', 'Flux integrated from map', 'Kron flux in elliptical aperture', 'From new raw data'), ('NUV (GALEX) AB', 20.533100000000001, '+/-0.118782', 'mag', 1290000000000000.0, 2.2200000000000001e-05, '+/-2.43E-06', 'Jy', '2012GASC..C...0000S', 'uncertainty', '2315.7     A', 'Broad-band measurement', 'Flux in fixed aperture', 'Flux in 7.5 arcsec diameter aperture', 'From new raw data'), ('u (SDSS PSF) AB', 18.254000000000001, '+/-0.021', 'asinh mag', 836000000000000.0, 0.00018799999999999999, '+/-3.63E-06', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '3585       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('u (SDSS Model) AB', 16.298999999999999, '+/-0.007', 'asinh mag', 836000000000000.0, 0.00114, '+/-7.85E-06', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '3585       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('u (SDSS CModel) AB', 15.965, '', 'asinh mag', 836000000000000.0, 0.00149, '', 'Jy', '2005SDSS4.C...0000:', 'no uncertainty reported', '3585       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('u (SDSS Petrosian)AB', 16.199000000000002, '+/-0.099', 'asinh mag', 836000000000000.0, 0.00125, '+/-1.14E-04', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '3585       A', 'Broad-band measurement', 'Modelled datum', 'SDSS Petrosian Radius =    4.34".         \nSDSS flags: CHILD - object is part of a blended parent object; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent\'s children;', 'From new raw data'), ('u (SDSS PSF) AB', 18.254000000000001, '+/-0.021', 'asinh mag', 836000000000000.0, 0.00018799999999999999, '+/-3.64E-06', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '3585       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('u (SDSS CModel) AB', 15.965, '', 'asinh mag', 836000000000000.0, 0.0015499999999999999, '', 'Jy', '2007SDSS6.C...0000:', 'no uncertainty reported', '3585       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('u (SDSS Model) AB', 16.298999999999999, '+/-0.007', 'asinh mag', 836000000000000.0, 0.00114, '+/-7.34E-06', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '3585       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('m_p', 14.699999999999999, '+/-0.4', 'mag', 681000000000000.0, 0.00562, '+/-2.50E-03', 'Jy', '1966CGCG3.C...0000Z', 'rms noise', '4400       A', 'Broad-band measurement', 'Estimated by eye', '', 'From new raw data'), ('B (m_B)', 14.800000000000001, '+/-0.30', 'mag', 681000000000000.0, 0.0051200000000000004, '+/-1.63E-03', 'Jy', '1991RC3.9.C...0000d', 'rms uncertainty', '4400       A', 'Broad-band measurement', 'Multiple methods', '', 'Homogenized from new and previously published data; StandardJohnson UBVRI filters assumed'), ('g (SDSS PSF) AB', 16.597000000000001, '+/-0.019', 'asinh mag', 617000000000000.0, 0.000834, '+/-1.43E-05', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '4858       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('g (SDSS Model) AB', 14.426, '+/-0.002', 'asinh mag', 617000000000000.0, 0.0061599999999999997, '+/-1.19E-05', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '4858       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('g (SDSS CModel) AB', 14.372, '', 'asinh mag', 617000000000000.0, 0.0064700000000000001, '', 'Jy', '2005SDSS4.C...0000:', 'no uncertainty reported', '4858       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('g (SDSS Petrosian)AB', 14.423, '+/-0.088', 'asinh mag', 617000000000000.0, 0.0061799999999999997, '+/-5.03E-04', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '4858       A', 'Broad-band measurement', 'Modelled datum', 'SDSS Petrosian Radius =    3.69".         \nSDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent\'s children;', 'From new raw data'), ('g (SDSS PSF) AB', 16.597000000000001, '+/-0.019', 'asinh mag', 617000000000000.0, 0.000834, '+/-1.46E-05', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '4858       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('g (SDSS CModel) AB', 14.372, '', 'asinh mag', 617000000000000.0, 0.0064700000000000001, '', 'Jy', '2007SDSS6.C...0000:', 'no uncertainty reported', '4858       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('g (SDSS Model) AB', 14.426, '+/-0.002', 'asinh mag', 617000000000000.0, 0.0061599999999999997, '+/-1.13E-05', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '4858       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('r (SDSS Model) AB', 13.6, '+/-0.002', 'asinh mag', 477000000000000.0, 0.0132, '+/-2.36E-05', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '6290       A', 'Broad-band measurement', 'Modelled datum', 'SDSS Effective Radius =   2.37" x   0.72".;SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent\'s children; CANONICAL_BAND - this band was primary (usually r);', 'From new raw data'), ('r (SDSS PSF) AB', 15.593999999999999, '+/-0.014', 'asinh mag', 477000000000000.0, 0.0020999999999999999, '+/-2.70E-05', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '6290       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children; CANONICAL_BAND - this band was primary (usually r);", 'From new raw data'), ('r (SDSS CModel) AB', 13.632, '', 'asinh mag', 477000000000000.0, 0.012800000000000001, '', 'Jy', '2005SDSS4.C...0000:', 'no uncertainty reported', '6290       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children; CANONICAL_BAND - this band was primary (usually r);", 'From new raw data'), ('r (SDSS Petrosian)AB', 13.615, '+/-0.081', 'asinh mag', 477000000000000.0, 0.012999999999999999, '+/-9.73E-04', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '6290       A', 'Broad-band measurement', 'Modelled datum', 'SDSS Petrosian Radius =    3.42".         \nSDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent\'s children; CANONICAL_BAND - this band was primary (usually r);', 'From new raw data'), ('r (SDSS PSF) AB', 15.593999999999999, '+/-0.014', 'asinh mag', 477000000000000.0, 0.0020999999999999999, '+/-2.71E-05', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '6290       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children; CANONICAL_BAND - this band was primary (usually r);", 'From new raw data'), ('r (SDSS CModel) AB', 13.632, '', 'asinh mag', 477000000000000.0, 0.012800000000000001, '', 'Jy', '2007SDSS6.C...0000:', 'no uncertainty reported', '6290       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children; CANONICAL_BAND - this band was primary (usually r);", 'From new raw data'), ('r (SDSS Model) AB', 13.6, '+/-0.002', 'asinh mag', 477000000000000.0, 0.0132, '+/-2.43E-05', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '6290       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children; CANONICAL_BAND - this band was primary (usually r);", 'From new raw data'), ('i (SDSS PSF) AB', 15.307, '+/-0.014', 'asinh mag', 389000000000000.0, 0.0027399999999999998, '+/-3.43E-05', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '7706       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('i (SDSS Petrosian)AB', 13.218, '+/-0.085', 'asinh mag', 389000000000000.0, 0.0229, '+/-1.80E-03', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '7706       A', 'Broad-band measurement', 'Modelled datum', 'SDSS Petrosian Radius =    3.40".         \nSDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent\'s children;', 'From new raw data'), ('i (SDSS Model) AB', 13.196, '+/-0.002', 'asinh mag', 389000000000000.0, 0.019099999999999999, '+/-3.41E-05', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '7706       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('i (SDSS PSF) AB', 15.307, '+/-0.014', 'asinh mag', 389000000000000.0, 0.0027399999999999998, '+/-3.53E-05', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '7706       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('i (SDSS CModel) AB', 13.226000000000001, '', 'asinh mag', 389000000000000.0, 0.018599999999999998, '', 'Jy', '2007SDSS6.C...0000:', 'no uncertainty reported', '7706       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('i (SDSS CModel) AB', 13.226000000000001, '', 'asinh mag', 389000000000000.0, 0.018599999999999998, '', 'Jy', '2005SDSS4.C...0000:', 'no uncertainty reported', '7706       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('i (SDSS Model) AB', 13.196, '+/-0.002', 'asinh mag', 389000000000000.0, 0.019099999999999999, '+/-3.52E-05', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '7706       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('z (SDSS PSF) AB', 15.028, '+/-0.013', 'asinh mag', 325000000000000.0, 0.00347, '+/-4.31E-05', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '9222       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('z (SDSS CModel) AB', 12.882, '', 'asinh mag', 325000000000000.0, 0.025499999999999998, '', 'Jy', '2005SDSS4.C...0000:', 'no uncertainty reported', '9222       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('z (SDSS PSF) AB', 15.028, '+/-0.013', 'asinh mag', 325000000000000.0, 0.00347, '+/-4.16E-05', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '9222       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('z (SDSS Model) AB', 12.879, '+/-0.002', 'asinh mag', 325000000000000.0, 0.025100000000000001, '+/-5.52E-05', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '9222       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('z (SDSS CModel) AB', 12.882, '', 'asinh mag', 325000000000000.0, 0.025100000000000001, '', 'Jy', '2007SDSS6.C...0000:', 'no uncertainty reported', '9222       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('z (SDSS Petrosian)AB', 12.938000000000001, '+/-0.078', 'asinh mag', 325000000000000.0, 0.023800000000000002, '+/-1.71E-03', 'Jy', '2005SDSS4.C...0000:', 'based on count statistics only', '9222       A', 'Broad-band measurement', 'Modelled datum', 'SDSS Petrosian Radius =    3.16".         \nSDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent\'s children;', 'From new raw data'), ('z (SDSS Model) AB', 12.879, '+/-0.002', 'asinh mag', 325000000000000.0, 0.025100000000000001, '+/-4.63E-05', 'Jy', '2007SDSS6.C...0000:', 'based on count statistics only', '9222       A', 'Broad-band measurement', 'Modelled datum', "SDSS flags: CHILD - object is part of a blended parent object; COSMIC_RAY - contains a cosmic ray pixel; BAD_RADIAL - some low S/N radial points; INTERP - object contains interpolated-over pixels; BINNED1 - detected at >=5 sigma in original imaging frame; BINNED_CENTER - image was binned while centroiding; BRIGHTEST_GALAXY_CHILD - brightest child among one parent's children;", 'From new raw data'), ('J_total (2MASS)', 11.481, '+/-0.019', 'mag', 240000000000000.0, 0.0407, '+/-7.18E-04', 'Jy', '20032MASX.C.......:', '1 sigma uncert.', '1.25      microns', 'Broad-band measurement', 'Total flux', '', 'From new raw data'), ('J_14arcsec (2MASS)', 11.962, '+/-0.007', 'mag', 240000000000000.0, 0.026100000000000002, '+/-1.69E-04', 'Jy', '20032MASX.C.......:', '1 sigma uncert.', '1.25      microns', 'Broad-band measurement', 'Flux in fixed aperture', '14.0 x 14.0 arcsec aperture.', 'From new raw data'), ('J_K20 (2MASS)', 11.58, '+/-0.016', 'mag', 240000000000000.0, 0.037100000000000001, '+/-5.52E-04', 'Jy', '20032MASX.C.......:', '1 sigma uncert.', '1.25      microns', 'Broad-band measurement', 'Flux integrated from map', '48.8 x   15.6 arcsec integration area.', 'From new raw data'), ('H_K20 (2MASS)', 10.856, '+/-0.022', 'mag', 182000000000000.0, 0.0465, '+/-9.53E-04', 'Jy', '20032MASX.C.......:', '1 sigma uncert.', '1.65      microns', 'Broad-band measurement', 'Flux integrated from map', '48.8 x   15.6 arcsec integration area.', 'From new raw data'), ('H_total (2MASS)', 10.744999999999999, '+/-0.022', 'mag', 182000000000000.0, 0.0516, '+/-1.06E-03', 'Jy', '20032MASX.C.......:', '1 sigma uncert.', '1.65      microns', 'Broad-band measurement', 'Total flux', '', 'From new raw data'), ('H_14arcsec (2MASS)', 11.244, '+/-0.009', 'mag', 182000000000000.0, 0.032599999999999997, '+/-2.71E-04', 'Jy', '20032MASX.C.......:', '1 sigma uncert.', '1.65      microns', 'Broad-band measurement', 'Flux in fixed aperture', '14.0 x 14.0 arcsec aperture.', 'From new raw data'), ('K_s', 10.608000000000001, '+/-0.029', 'mag', 138000000000000.0, 0.038100000000000002, '+/-1.03E-03', 'Jy', '20032MASX.C.......:', '1 sigma uncert.', '2.17      cm', 'Broad-band measurement', 'Flux in fixed aperture', '48.8 x   15.6 arcsec integration area.', 'From new raw data'), ('K_s_14arcsec (2MASS)', 10.960000000000001, '+/-0.010', 'mag', 138000000000000.0, 0.0275, '+/-2.55E-04', 'Jy', '20032MASX.C.......:', '1 sigma uncert.', '2.17      microns', 'Broad-band measurement', 'Flux in fixed aperture', '14.0 x 14.0 arcsec aperture.', 'From new raw data'), ('K_s', 10.538, '+/-0.034', 'mag', 138000000000000.0, 0.040599999999999997, '+/-1.29E-03', 'Jy', '20032MASX.C.......:', '1 sigma uncert.', '2.17      cm', 'Broad-band measurement', 'Flux in fixed aperture', '', 'From new raw data'), ('1.4GHz', 2.5, '+/-0.14', 'milliJy', 1400000000.0, 0.0025000000000000001, '+/-1.35E-04', 'Jy', '1995ApJ...450..559B', 'rms noise', '1.4        GHz', 'Broad-band measurement; synthetic band', 'Peak flux', '', 'From new raw data; Corrected for contaminating sources'), ('1.4GHz', 2.5800000000000001, '', 'milliJy', 1400000000.0, 0.0025799999999999998, '', 'Jy', '1995ApJ...450..559B', 'no uncertainty reported', '1.4        GHz', 'Broad-band measurement; synthetic band', 'From fitting to map', '', 'From new raw data; Corrected for contaminating sources')]
    
    best = choose_diameter_measurement(options)
    
    print 'best: ',best
    
    print
    print 'Options were: '
    for i in options:
        print i
    
    print
    print 'done'
    
if __name__ == '__main__':
    main()
    