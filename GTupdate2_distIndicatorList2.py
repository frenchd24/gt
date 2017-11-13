#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_distIndicatorList2.py, v1 10/17/2017

List all the distIndicators in NED-D_edit.csv

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


# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################


def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':
        filename = '/usr/data/moosejaw/frenchd/GT_update2/NED-D_edit.csv'
    
    elif user =='David':
        filename = '/Users/David/Research_Documents/GT_update2/NED-D_edit.csv'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the files
    csv.field_size_limit(sys.maxsize)
    
    # basic info file
    basicFile = open(filename,'rU')
    basicReader = csv.DictReader(basicFile)
  
    d = {}
      
    count = 0
    # loop through and do the work
    for line in basicReader:
        # update the counter
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
        
        # stuff in basicReader:
        method = line['Method']
            
        if d.has_key(method):
            n = d[method]
            n +=1
            d[method] = n
        else:
            d[method] = 1

    print 'Results - lumClass: ',
    for i in d:
        print i

    basicFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
