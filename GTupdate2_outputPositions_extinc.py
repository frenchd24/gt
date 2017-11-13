#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_outputPositions_extinc.py, v1 04/14/2017

Take all the ra and dec values from return_basic_full.csv
and output an ascii table with them, separated by only spaces (i.e., [ra dec] per line)

output: all_ra_dec.txt

Also did rejected_results_redo.csv (09/08/17)

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

def absoluteMag_noExtinc(m,d):
    # m is apparent magnitude, d is distance in Mpc
    M = float(m) - 5*math.log10((float(d)*10**6)/10)
    return M


def absoluteMag(m,d,e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)
    return M


def printOutInfo(line,t):
    # grab all the information in the table for an object and return it all if t == f, 
    # or just the basics if t == b
    
    preferredName = line['preferredName']
    oldName = line['oldName']
    redshift = line['redshift']
    degreesJ2000 = eval(line['degreesJ2000RA_Dec'])
    J2000 = eval(line['J2000RA_Dec'])
    gLongLat = eval(line['galacticLong_Lat'])
    radVel = line['radialVelocity (km/s)']
    vcorr = line['vcorr (km/s)']
    distvcorr = line['distvcorr (Mpc)']
    rid = eval(line['rIndependentDistMean_sd_min_max (Mpc)'])
    bestDist = line['Best Distance (Mpc)']
    angDiameter = line['angDiameters (arcsec)']
    linDiameter = eval(line['linDiameters (kpc)'])
    inclination = line['inclination (deg)']
    positionAngle = line['positionAngle (deg)']
    diameterKey = line['diameterKey']
    RC3flag = line['RC3flag']
    RC3type = line['RC3type']
    RC3inc = line['RC3inc (deg)']
    RC3pa = line['RC3pa (deg)']
    Groups = eval(line['Groups_Dist_std (Mpc)'])
    groupsInfo = line['groupsInfo']
    morphology = line['morphology']
    distanceIndicator = line['distanceIndicator']
    luminosityClass = line['luminosityClass']
    EBminusV = line['EBminusV']
    EBminusV_new = line['E(B-V)_new']
    B_median = line['B_median']
    B_sdss_median = line['B_sdss_median']
    B_median_Lstar = line['B_median_Lstar']
    B_max_Lstar = line['B_max_Lstar']
    B_sdss_median_Lstar = line['B_sdss_median_Lstar']
    B_sdss_max_Lstar = line['B_sdss_median_Lstar']
    Lstar = line['Lstar']
    photometry = eval(line['photometry'])
    altNames = eval(line['alternativeNames'])


    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full.csv'
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/all_ra_dec'
        basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo.csv'
        outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_ra_dec'
        
    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the files
    basicFile = open(basicFilename,'rU')
    basicReader = csv.DictReader(basicFile)
    
    
    outputFile = open(outputFilename+'{0}-{1}.txt'.format('0','20000'),'wt')
    
    # write the header
    outputFile.write('| ra       | dec      |\n')
    outputFile.write('| double   | double   |\n')
    
    # start the count on two because of the headers
    count = 2
    maxCount = 20000
    
    # loop through and do the work
    for bline in basicReader:
        
        if count % maxCount == 0:
            # close the old file
            outputFile.close()
            
            # open a new one
            outputFile = open(outputFilename+'{0}-{1}.txt'.format(count,count+maxCount),'wt')
    
            # write the header
            outputFile.write('| ra       | dec      |\n')
            outputFile.write('| double   | double   |\n')
            count +=2
        
        # grab the ra and dec
        ra,dec = eval(bline['degreesJ2000RA_Dec'])
        
        # write them to file
        outputFile.write('{0} {1}\n'.format(ra,dec))
    
        count +=1
        
    basicFile.close()
    outputFile.close()

    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
