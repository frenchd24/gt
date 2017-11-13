#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: PA_ascii_read_test.py, v1 10/04/17

Test how to read NEDs ascii_bar output format as a replacement for XML.
    
'''

import sys
# import os
import csv
# import string
# import warnings
import numpy as np
# import atpy
import getpass
from utilities import *
import math
from astropy.io import ascii
from astropy import table
import urllib2


# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################

    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':
        inFilename = '/usr/data/moosejaw/frenchd/GT_update2/PA_ascii_test.txt'
        inFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/PA_ascii_test_out.csv'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    host = "http://ned.ipac.caltech.edu/cgi-bin/datasearch?objname={0}&search_type=Diameters&of=ascii_bar".format('NGC1574')
    url = urllib2.urlopen(host)
        
    # the header/column names
    fieldnames = (\
    'No.',\
    'Frequency targeted',\
    'Refcode',\
    'Major Axis',\
    'Major Axis Flag',\
    'Major Axis Unit',\
    'Minor Axis',\
    'Minor Axis Flag',\
    'Minor Axis Unit',\
    'Axis Ratio',\
    'Axis Ratio Flag',\
    'Major Axis Uncertainty',\
    'Ellipticity',\
    'Eccentricity',\
    'Position Angle',\
    'Equinox',\
    'Reference Level',\
    'NED Frequency',\
    'NED Major Axis',\
    'NED Major Axis Uncertainty',\
    'NED Axis Ratio',\
    'NED Ellipticity',\
    'NED Eccentricity',\
    'NED cos-1_axis_ratio',\
    'NED Position Angle',\
    'NED Minor Axis',\
    'Minor Axis Uncertainty',\
    'NED Minor Axis Uncertainty',\
    'Axis Ratio Uncertainty',\
    'NED Axis Ratio Uncertainty',\
    'Ellipticity Uncertainty',\
    'NED Ellipticity Uncertainty',\
    'Eccentricity Uncertainty',\
    'NED Eccentricity Uncertainty',\
    'Position Angle Uncertainty',\
    'NED Position Angle Uncertainty',\
    'Significance',\
    'Frequency',\
    'Frequency Unit',\
    'Frequency Mode',\
    'Detector Type',\
    'Fitting Technique',\
    'Features',\
    'Measured Quantity',\
    'Measurement Qualifiers',\
    'Targeted RA',\
    'Targeted DEC',\
    'Targeted Equinox',\
    'NED Qualifiers',\
    'NED Comment')
        
    # open the files
    csv.field_size_limit(sys.maxsize)
    
    # basic info file
#     inFile = open(inFilename,'rU')
#     reader = csv.DictReader(inFile,'|')
    
    inFile = url
    
    reader = ascii.read(inFile, data_start=9, delimiter='|')

    
#     t_altNames = table.Table(list(reader_altNames), names=fieldnames_altNames,dtype=(d['Name'],'U2038'))
# 
#     ascii.write(t_altNames, outFilename_altNames,format='fixed_width',overwrite=True,delimiter=' ',\
#     bookend=False,delimiter_pad=None,formats={'Name':'%-30s','altNames':'%-2038s'})

    print 'Here it is: '
    print
#     print 'NED Frequency vs Frequency Targeted'

    for i in reader:
        print 'Frequency Targeted: ',i[1]
        print "NED Major: ",i[18]
        print 'NED Minor: ',i[25]
        print 'NED Ratio: ',i[20]
        print 'NED PA: ',i[24]
        print
    
    # close the files
    inFile.close()

    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
