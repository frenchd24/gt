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

    options = [['POSS1 103a-O', (1380.0, 480.0), 0.35, 162.0], ['POSS1 103a-O', (270.0, 33.0), 0.12, None], ['POSS1 103a-O', (1350.0, 372.0), 0.28, None], ['RC3 D_25, R_25 (blue)', (1312.7, 434.5), 0.33, 162.0], ['RC3 A_e (Johnson B)', (476.6, None), None, None], ['RC3 D_0 (blue)', (1312.7, None), None, None], ['K_s (LGA/2MASS isophotal)', (146.4, 146.4), 1.0, 16.0], ['K_s (LGA/2MASS "total")', (250.2, 250.2), 1.0, 16.0]]
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
    