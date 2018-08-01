#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_countDiamSource.py, v1.1 07/07/18

Count how many galaxy ratios, PAs and diameters came from each survey (09/26/17)

v1.1: Update for thesis (07/07/18)

'''

import sys
import os
import csv
# import string
import warnings
import numpy
# import atpy
import getpass
from utilities import *
import math
from astropy.io import ascii
from astropy import table

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################


    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':
        inFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable13_filtered.csv'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2_group2.csv'
        
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5.dat'
#         outFilename2 = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5_altnames.dat'
#         outFilename3 = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5_altNamesLEFT.dat'


    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the files
    csv.field_size_limit(sys.maxsize)
    
    # basic info file
    inFile = open(inFilename,'rU')
    reader = csv.DictReader(inFile)
    
    d = {}
    
    d_general = {}
    
    count = 0
    RID_count = 0
    # loop through and do the work
    for line in reader:
        # update the counter
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
        
        MType = line['MType'].lower()
        
        if d_general.has_key(MType):
            prev = d_general[MType]
            prev +=1
            d_general[MType] = prev
        else:
            d_general[MType] = 1
            
        
        if MType[0] == 's' and MType[:2] != 's0':
            if d.has_key('s'):
                c = d['s']
                c+=1
                d['s'] = c
            else:
                d['s'] = 1
                
        elif bfind(MType, 'sb') or bfind(MType, 'sc') or bfind(MType, 'sd') or bfind(MType, 'sa'):
            if d.has_key('s'):
                c = d['s']
                c+=1
                d['s'] = c
            else:
                d['s'] = 1
                
        elif MType[0] == 'e':
            if d.has_key('e'):
                c = d['e']
                c+=1
                d['e'] = c
            else:
                d['e'] = 1

        elif MType[:2] == 's0':
            if d.has_key('s0'):
                c = d['s0']
                c+=1
                d['s0'] = c
            else:
                d['s0'] = 1

        elif MType[0] == 'i':
            if d.has_key('i'):
                c = d['i']
                c+=1
                d['i'] = c
            else:
                d['i'] = 1
                
        else:
            pass


            
    print
    print 'Total d_general: ',len(d_general)
    print 'Total d: ',len(d)
    print
    print
        
    dkeys = d.keys()
    dvals = d.values()
        
    d_generalkeys = d_general.keys()
    d_generalvals = d_general.values()
    
    for k,v in zip(dkeys,dvals):
        print '{0} = {1}'.format(k,v)
    
    print
    print
    print

    for k,v in zip(d_generalkeys,d_generalvals):
        print '{0} = {1}'.format(k,v)
    
    print
    print
    print
    
    # close the files
    inFile.close()

    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
