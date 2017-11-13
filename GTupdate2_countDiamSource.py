#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_countDiamSource.py, v1.0 09/26/17

Count how many galaxy ratios, PAs and diameters came from each survey

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


        inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2_group2.csv'
        
        outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5.dat'
        outFilename2 = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5_altnames.dat'
        outFilename3 = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5_altNamesLEFT.dat'


    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the files
    csv.field_size_limit(sys.maxsize)
    
    # basic info file
    inFile = open(inFilename,'rU')
    reader = csv.DictReader(inFile)
    
    d = {'diameterKeys':{},'ratioKey':{}, 'paKey':{}}
    
    count = 0
    # loop through and do the work
    for line in reader:
        # update the counter
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
        
        diameterKey = line['diameterKey']
        ratioKey = line['ratioKey']
        paKey = line['paKey']
        
        if d['diameterKeys'].has_key(diameterKey):
            prev = d['diameterKeys'][diameterKey]
            prev +=1
        else:
            d['diameterKeys'][diameterKey] = 1
            
        if d['ratioKey'].has_key(ratioKey):
            prev = d['ratioKey'][ratioKey]
            prev +=1
        else:
            d['ratioKey'][ratioKey] = 1
            
        if d['paKey'].has_key(paKey):
            prev = d['paKey'][paKey]
            prev +=1
        else:
            d['paKey'][paKey] = 1
            
        
    diameterKeys = d['diameterKeys'].keys()
    diameterVals = d['diameterKeys'].values()
    
    ratioKeyKeys = d['ratioKeyKeys'].keys()
    ratioKeyVals = d['ratioKeyKeys'].values()

    paKeys = d['paKeys'].keys()
    paVals = d['paKeys'].values()
    
    print 'For Diameter keys: '
    for k,v in zip(diameterKeys,diameterVals):
        print '{0} = {1}'.format(k,v)
    
    print 
    print 'For Ratio keys: '
    for k,v in zip(ratioKeys,ratioVals):
        print '{0} = {1}'.format(k,v)
        
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
