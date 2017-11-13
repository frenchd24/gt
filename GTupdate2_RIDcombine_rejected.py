#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_RIDcombine_rejected.py, v 1 09/22/2017

Do the below on this file: rejected_results_redo_extinc_rc3.csv
Output: rejected_results_redo_extinc_rc3_rid.csv


Add redshift independent distances to info to return_basic_full_extinc_rc3.csv

Output: return_basic_full_extinc_rc3_rid.csv (09/20/2017)

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
        basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_extinc_rc3.csv'
        ridFilename = '/usr/data/moosejaw/frenchd/GT_update2/NED-D_edit.csv'
        outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_extinc_rc3_rid.csv'

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
    'RID',\
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
    
    # RID file
    ridFile = open(ridFilename,'rU')
    ridReader = csv.DictReader(ridFile)
    
    # RID file fieldnames:
    # Exclusion Code    D   G   Galaxy ID   m-M err D (Mpc) Method  REFCODE SN ID   redshift (z)    Hubble const.   Adopted LMC modulus Date (Yr. - 1980)   Notes
    
    d= {}
    for r in ridReader:
        name = r['Galaxy ID'].replace(' ','')
        dist = r['D (Mpc)']
        err = r['err']
        
#         print 'name: ',name
#         print 'dist: ',dist
#         print 'err:',err
        
        if not isNumber(err) or err == 0 or err == '0':
            err = 0.5

        if d.has_key(name):
            n = d[name]
            
            nd = n['dist']
            ne = n['err']
            
            nd.append(float(dist))
            ne.append(float(err))
                    
        else:
            d[name] = {'dist':[float(dist)],'err':[float(err)]}
    
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
        RID = eval(bline['rIndependentDistMean_sd_min_max (Mpc)'])
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
        
        meanDist = 'x'
        medianDist = 'x'
        stdDist = 'x'
        minDist = 'x'
        maxDist = 'x'
        
        for a in altNames:
            if d.has_key(a):
                n = d[a]
                dists = n['dist']
                errs = n['err']
                
                try:
                    meanDist = round(np.mean(dists),3)
                    medianDist = round(np.median(dists),3)
                    stdDist = round(np.std(dists),3)
                    minDist = min(dists)
                    maxDist = max(dists)
                    
#                     print
#                     print 'm-M = ',5*math.log10(meanDist*10**6 / 10)
                    
                    if len(dists) <2:
                        mM = 5*math.log10(meanDist*10**6 / 10)
                        stdDist = round(meanDist - (10**((mM - float(errs[0]))/5 + 1))/10**6,3)
                        
                except Exception,e:
                    print 'error: ',e
                    print 'dists = ',dists
                    print
                    
                break
        
#         if bpreferredName == 'LEOIVDWARFSPHEROIDAL':
#             print
#             print 'RID[0]: ',RID[0]
#             print 'isNumber(RID[0]) : ',isNumber(RID[0])
#             print 'isNumber(meanDist): ',isNumber(meanDist)
#             print
#             
        if isNumber(RID[0]) and not isNumber(meanDist):
            RID_mean = RID[0]
            RID_std = RID[1]
            RID_min = RID[2]
            RID_max = RID[3]
            
            meanDist = RID_mean
            medianDist = RID_mean
            stdDist = RID_std
            minDist = RID_min
            maxDist = RID_max
            
            print 'preferred Name : ',bpreferredName
            
        # now write it to file
        objectInfoList = (\
        bline['preferredName'],\
        bline['oldName'],\
        bline['redshift'],\
        bline['degreesJ2000RA_Dec'],\
        bline['J2000RA_Dec'],\
        bline['galacticLong_Lat'],\
        bline['rIndependentDistMean_sd_min_max (Mpc)'],\
        [meanDist,medianDist,stdDist,minDist,maxDist],\
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
