 #!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  add_extinction_column.py, v 2.0 04/24/17

Add extinction column to return_basic_full.csv based on new IRSA values

Returns: return_basic_full_extinc.csv

** Successfully made /usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc.csv
    on 04/24/17
**

ADOPTED FROM: add_extinction_column.py, v 1.0 02/06/14 = 
"
Adds a column in NewGalaxyTable2.csv for E(B-V) based on values from IRSA.

**Successfully made NewGalaxyTable3.csv on 02/13/15**
"

'''

import sys
import os
import csv
# import string
# import warnings
# import urllib
# import numpy
from pylab import *
# import atpy
import math
import getpass
import itertools
from utilities import *

# import scipy.optimize as optimization
# import pickle


# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.ticker import NullFormatter


def schechter(m,phi,mstar,alpha):
    # construct a Schechter luminosity function and return the associated density function
    
#     s = 0.4*log(10)*phi*(10**(0.4*(mstar-m)))**(alpha +1) * exp(-10**(0.4*(mstar-m)))
#     s = 0.4*log(10)*phi*(10**(0.4*(m-mstar)*(alpha +1))) * exp(-10**(0.4*(m-mstar)))
    s = 0.4 * math.log(10) * phi * (10**(0.4*(mstar-m)))**(alpha +1) * exp(-10**(0.4*(mstar-m)))

    return s
    
    
def absoluteMag(m,d):
    # m is apparent magnitude, d is distance in Mpc
    
    M = float(m) - 5*log10((float(d)*10**6)/10)
    return M


def lstarValue(mstar,m):
    # calculate and return L/Lstar
    
    lratio = 10**(-0.4*(m-mstar))
    return lratio


def readlines_into_lists(file, header):
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
            splitLine = l.split()
            outList.append(splitLine)
        
    return outList


def main():

    # open and read the galaxy table    
    if getpass.getuser() == 'frenchd':
        inputFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full.csv'
        outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc.csv'
        extincTableDirectory = '/usr/data/moosejaw/frenchd/GT_update2/extincTables/'
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # open the galaxy file
    inputFile = open(inputFilename,'rU')
    reader = csv.DictReader(inputFile)
        
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
    'photometry',\
    'alternativeNames')

    writerOutFile = open(outputFilename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    
    header = ['\\','|']
    # open output
    file1 = open('{0}/0-20000_extinc.txt'.format(extincTableDirectory),'rU')
    lines1 = readlines_into_lists(file1, header)
#     print 'lines1: ',lines1
#     print
#     print
#     print
#     print
#     print
#     print
    
    file2 = open('{0}/20000-40000_extinc.txt'.format(extincTableDirectory),'rU')
    lines2 = readlines_into_lists(file2, header)
    
    file3 = open('{0}/40000-60000_extinc.txt'.format(extincTableDirectory),'rU')
    lines3 = readlines_into_lists(file3, header)
    
    file4 = open('{0}/60000-80000_extinc.txt'.format(extincTableDirectory),'rU')
    lines4 = readlines_into_lists(file4, header)

    file5 = open('{0}/80000-100000_extinc.txt'.format(extincTableDirectory),'rU')
    lines5 = readlines_into_lists(file5, header)
    
    file6 = open('{0}/100000-120000_extinc.txt'.format(extincTableDirectory),'rU')
    lines6 = readlines_into_lists(file6, header)
    
    file7 = open('{0}/120000-140000_extinc.txt'.format(extincTableDirectory),'rU')
    lines7 = readlines_into_lists(file7, header)
    
    allLines = lines1+lines2+lines3+lines4+lines5+lines6+lines7
    
    count = 0
#     for gline, eline in zip(reader, itertools.chain(lines1,lines2,lines3,lines4,lines5,lines6)):
    for gline, eline in zip(reader, allLines):
        ra,dec = eval(gline['degreesJ2000RA_Dec'])
        
        era,edec = eline[0],eline[1]
        EBminusV_SF = eline[3]
        meanEBminusV_SF = eline[4]
        
        # convert ra,dec to same format as the extinction values
        rat = trunc(ra,9)
        lenRat = len(rat)
        if lenRat == 3:
            rat = str(rat)+'000000'
        if lenRat == 4:
            rat = str(rat)+'00000'
        if lenRat ==5:
            rat = str(rat)+'0000'
        if lenRat == 6:
            rat = str(rat)+'000'
        if lenRat == 7:
            rat = str(rat)+'00'
        if lenRat == 8:
            rat = str(rat)+'0'
        else:
            rat = str(rat)
            
        dect = trunc(dec,9)
        dect2 = trunc(str(dect),5)
        rat2 = trunc(str(rat),5)
        
        erat = trunc(str(era),5)
        edect = trunc(str(edec),5)
        
#         if count <20000:
#             lines1
#         
#         if count>=20000 and count<40000:
#             file2.write(' {0}   {1}\n'.format(rat,dect))
# 
#         if count>=40000 and count<60000:
#             file3.write(' {0}   {1}\n'.format(rat,dect))
#             
#         if count>=60000 and count<80000:
#             file4.write(' {0}   {1}\n'.format(rat,dect))
#             
#         if count>=80000 and count <100000:
#             file5.write(' {0}   {1}\n'.format(rat,dect))
# 
#         if count>=100000:
#             file6.write(' {0}   {1}\n'.format(rat,dect))
            
        count+=1
        if erat == rat2 and edect == dect2:
#             print 'match ',count
            objectInfoList = (\
            gline['preferredName'],\
            gline['oldName'],\
            gline['redshift'],\
            gline['degreesJ2000RA_Dec'],\
            gline['J2000RA_Dec'],\
            gline['galacticLong_Lat'],\
            gline['rIndependentDistMean_sd_min_max (Mpc)'],\
            gline['morphology'],\
            gline['distanceIndicator'],\
            gline['luminosityClass'],\
            meanEBminusV_SF,\
            gline['radialVelocity (km/s)'],\
            gline['vcorr (km/s)'],\
            gline['angDiameters (arcsec)'],\
            gline['linDiameters (kpc)'],\
            gline['distvcorr (Mpc)'],\
            gline['inclination (deg)'],\
            gline['photometry'],\
            gline['alternativeNames'])
              
       
            row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
            writer.writerow(row)
        
            # update counter
            percentComplete = round((float(count)/130759)*100,2)
            sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
            sys.stdout.flush()
            
        else:
            print 'no match: {0},{1} != {2},{3}'.format(erat, edect, rat2, dect2)
        
    
    print 'count = ',count

    file1.close()
    file2.close()
    file3.close()
    file4.close()
    file5.close()
    file6.close()
    file7.close()


    inputFile.close()
    writerOutFile.close()


if __name__=="__main__":
    main()
    
    