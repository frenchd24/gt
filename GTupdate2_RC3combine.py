#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_RC3combine.py, v 1 05/24/2017

Add RC3 info to return_basic_full_extinc.csv

Output: return_basic_full_extinc_rc3.csv

ran on rejected_results_redo_extinc.csv (09/11/17)


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

def readlines_into_lists(file, header, delimiter):
    # take in a file. Reads each line and splits on spaces. returns a list of rows,
    # each of which is list of values
    
    # header is a list of characters that indicates header lines that should be skipped
    # example: header = ['\','|']
    
    outList = []
    lines = file.readlines()
    for l in lines:
        isHeader = False
        for c in header:
            # if a header character, 'c', is found at the start of the line, 'l'
            if str(l)[0] == c:
                isHeader = True
                print 'was header: ',str(l)
                
        if not isHeader:
            splitLine = l.split(delimiter)
            if len(splitLine)>4:
                outList.append(splitLine)
        
    return outList


    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc2.csv'
#         rc3Filename = '/usr/data/moosejaw/frenchd/GT_update2/rc3catalog_edit.txt'
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc_rc3.csv'

        basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_extinc.csv'
        rc3Filename = '/usr/data/moosejaw/frenchd/GT_update2/rc3catalog_edit.txt'
        outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_extinc_rc3.csv'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the files
    csv.field_size_limit(sys.maxsize)
    
    # basic info file
    basicFile = open(basicFilename,'rU')
    basicReader = csv.DictReader(basicFile)

    # new fieldnames for updated galaxy table
    fieldnames = ('preferredName',\
    'oldName',\
    'redshift',\
    'degreesJ2000RA_Dec',\
    'J2000RA_Dec',\
    'galacticLong_Lat',\
    'rIndependentDistMean_sd_min_max (Mpc)',\
    'morphology',\
    'distanceIndicator',\
    'luminosityClass',\
    'EBminusV',\
    'radialVelocity (km/s)',\
    'vcorr (km/s)',\
    'angDiameters (arcsec)',\
    'linDiameters (kpc)',\
    'distvcorr (Mpc)',\
    'inclination (deg)',\
    'RC3_type',\
    'RC3_d25 (arcsec)',\
    'RC3_r25',\
    'RC3_pa (deg)',\
    'alternativeNames')

    writerOutFile = open(outputFilename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    # rc3 file
    header = ['#']
    rc3File = open(rc3Filename,'rU')
    rc3Lines = readlines_into_lists(rc3File, header, ';')
    
    count = 0
    # loop through and do the work
    for bline in basicReader:
        count +=1
        # stuff in basicReader:
        bpreferredName = bline['preferredName']
        boldName = bline['oldName']
        z = bline['redshift']
        degreesJ2000RA_Dec = bline['degreesJ2000RA_Dec']
        J2000RA_Dec = bline['J2000RA_Dec']
        galacticLong_Lat = bline['galacticLong_Lat']
        RID = bline['rIndependentDistMean_sd_min_max (Mpc)']
        morph = bline['morphology']
        distanceIndicator = bline['distanceIndicator']
        lumClass = bline['luminosityClass']
        EBminusV = bline['EBminusV']
        helioVel = bline['radialVelocity (km/s)']
        vcorr = bline['vcorr (km/s)']
        angDiameters = bline['angDiameters (arcsec)']
        linDiameters = bline['linDiameters (kpc)']
        distvcorr = bline['distvcorr (Mpc)']
        inc = bline['inclination (deg)']
        altNames = eval(bline['alternativeNames'])
        
        rc3_diam = 'x'
        rc3_type = 'x'
        rc3_ratio = 'x'
        rc3_pa = 'x'
#         print 'alternativeNames: ',altNames
        
        for altName in altNames:
            if bfind(altName,'PGC'):
#                 print 'found one: ',altName
                for rline in rc3Lines:
                    # order in the RC3 file
                    # RA2000;DE2000;name;altname;PGC;type;T;lumcl;D25;R25;PA;BT;BT_code;BoT;V21;cz
        
                    rc3RA = rline[0]
                    rc3Dec = rline[1]
                    name = rline[2]
                    altname = rline[3]
                    pgcname = rline[4]
                    type = rline[5]
                    t = rline[6]
                    lumcl = rline[7]
                    d25 = rline[8]
                    r25 = rline[9]
                    pa = rline[10]
                    bt = rline[11]
                    bt_code = rline[12]
                    bot = rline[13]
                    v21 = rline[14]
                    cz = rline[15]
        
                    pgcname1,pgcname2 = pgcname.split('PGC')
                    pgcname2 = pgcname2.strip()

                    if len(pgcname2) == 1:
                        newPGCname = 'PGC00000'+pgcname2
        
                    elif len(pgcname2) == 2:
                        newPGCname = 'PGC0000'+pgcname2

                    elif len(pgcname2) == 3:
                        newPGCname = 'PGC000'+pgcname2
            
                    elif len(pgcname2) == 4:
                        newPGCname = 'PGC00'+pgcname2
            
                    elif len(pgcname2) == 5:
                        newPGCname = 'PGC0'+pgcname2
            
                    elif len(pgcname2) == 6:
                        newPGCname = 'PGC'+pgcname2
                    else:
                        newPGCname = pgcname
                    
#                     print newPGCname
        
                    if altName == newPGCname:
                        print 'found: ',newPGCname
                        if isNumber(d25):
                            rc3_diam = (10**(float(d25)-1))*60
                        else:
                            rc3_diam = 'x'
            
                        if isNumber(r25):
                            rc3_ratio = 1/(10**float(r25))
                        else:
                            rc3_ratio = 'x'
                            
                        rc3_pa = pa
                        rc3_type = type
                
                        break
                    
        # now write it to file
        objectInfoList = (\
        bline['preferredName'],\
        bline['oldName'],\
        bline['redshift'],\
        bline['degreesJ2000RA_Dec'],\
        bline['J2000RA_Dec'],\
        bline['galacticLong_Lat'],\
        bline['rIndependentDistMean_sd_min_max (Mpc)'],\
        bline['morphology'],\
        bline['distanceIndicator'],\
        bline['luminosityClass'],\
        bline['EBminusV'],\
        bline['radialVelocity (km/s)'],\
        bline['vcorr (km/s)'],\
        bline['angDiameters (arcsec)'],\
        bline['linDiameters (kpc)'],\
        bline['distvcorr (Mpc)'],\
        bline['inclination (deg)'],\
        rc3_type,\
        rc3_diam,\
        rc3_ratio,\
        rc3_pa,\
        altNames)
                      
        row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
        writer.writerow(row)

        # update counter
        percentComplete = round((float(count)/130759)*100,2)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        sys.stdout.flush()

    basicFile.close()
    writerOutFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
