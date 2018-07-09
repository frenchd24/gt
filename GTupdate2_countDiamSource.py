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
    
    d = {'diameterKey':{},'ratioKey':{}, 'paKey':{}}
    
    count = 0
    RID_count = 0
    # loop through and do the work
    for line in reader:
        # update the counter
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
        
        RID_median = line['RID_median']
        
        if float(RID_median) >=0:
            RID_count +=1
        
        diameterKey = line['diameterKey']
        ratioKey = line['ratioKey']
        paKey = line['paKey']
        
        if d['diameterKey'].has_key(diameterKey):
            prev = d['diameterKey'][diameterKey]
            prev +=1
            d['diameterKey'][diameterKey] = prev
        else:
            d['diameterKey'][diameterKey] = 1
            
        if d['ratioKey'].has_key(ratioKey):
            prev = d['ratioKey'][ratioKey]
            prev +=1
            d['ratioKey'][ratioKey] = prev
        else:
            d['ratioKey'][ratioKey] = 1
            
        if d['paKey'].has_key(paKey):
            prev = d['paKey'][paKey]
            prev +=1
            prev = d['paKey'][paKey] = prev
        else:
            d['paKey'][paKey] = 1
            
            
    print
    print 'Total: ',count
    print 'RID_count: ',RID_count
    print 'Percentage: ',float(RID_count)/float(count)
    print
    print
        
    diameterKeys = d['diameterKey'].keys()
    diameterVals = d['diameterKey'].values()
    
    ratioKeys = d['ratioKey'].keys()
    ratioVals = d['ratioKey'].values()

    paKeys = d['paKey'].keys()
    paVals = d['paKey'].values()
    
    print 'For Diameter keys: '
    for k,v in zip(diameterKeys,diameterVals):
        print '{0} = {1}'.format(k,v)
    
    print
    print
    print
    
    print 'For Ratio keys: '
    for k,v in zip(ratioKeys,ratioVals):
        print '{0} = {1}'.format(k,v)
        
    print
    print
    print
    print
    
    print 'For PA keys: '
    for k,v in zip(paKeys,paVals):
        print '{0} = {1}'.format(k,v)
    
    
    # close the files
    inFile.close()
#     outFile.close()
#     outFile2.close()

    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
