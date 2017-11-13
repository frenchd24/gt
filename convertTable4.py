#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: convertTable4.py  v 4  05/25/17

Convert the master galaxy table from dictionary csv to tab separated csv with no quotes
or parenthesis

Convert new table: '/usr/users/frenchd/gt/NewGalaxyTable.csv'

v3: Convert table: '/usr/users/frenchd/gt/NewGalaxyTable3.csv' into 'NewGalaxyTable3_alternative.csv'
    redo this after fixing MESSIER and NGC names in NewGalaxyTable3.csv - 03/06/15
    
v3.1: Convert table: '/usr/users/frenchd/gt/NewGalaxyTable5.csv'  - 06/01/15

v4: Convert GT_update2 table: '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable.csv'

"""

import getpass
import sys
import os
import csv
import string
from math import *
from utilities import *

# programName = 'convertTable2'
# version = '2.0'
# date = '04/08/14'
    
################################################################

# def isNumber(s):
#     try:
#         float(s)
#         return True
#     except ValueError:
#         return False
# 
# def trunc(x,n):
#     # truncates a number, x, to n places without rounding
#     if isNumber(x):
#         x = float(x)
#         slen = len('%.*f' %(n,x))
#         return str(x)[:slen]
#     else:
#         return x
    

def main():
    
    # Decide which directory to use
    user = getpass.getuser()
    if user == 'frenchd':
        galaxyFileName = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable.csv'
        outputFileName = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable_alt.csv'
    
    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()    
        
        
    galaxyFile = open(galaxyFileName,'rU')
    galaxyReader = csv.DictReader(galaxyFile)
    
    fieldnames = (\
    'preferredName',\
    'NEDname',\
    'z',\
    'degreesJ2000RA',\
    'degreesJ2000Dec',\
    'J2000RA',\
    'J2000Dec',\
    'galacticLong',\
    'galacticLat',\
    'helioVelocity (km/s)',\
    'vcorr (km/s)',\
    'distvcorr (Mpc)',\
    'zIndependentDist_mean (Mpc)',\
    'zIndependentDist_std (Mpc)',\
    'zIndependentDist_min (Mpc)',\
    'zIndependentDist_max (Mpc)',\
    'bestDistance (Mpc)',\
    'distErr (Mpc)',\
    'angMajor (arcsec)',\
    'angMinor (arcsec)',\
    'angMajor_err (arcsec)',\
    'angMinor_err (arcsec)',\
    'linMajor (kpc)',\
    'linMinor (kpc)',\
    'linMajor_err (kpc)',\
    'linMinor_err (kpc)',\
    'inc (deg)',\
    'adjustedInc (deg)',\
    'incErr (deg)',\
    'PA (deg)',\
    'diameterKey',\
    'ratioKey',\
    'paKey',\
    'RC3_type',\
    'RC3_d25 (arcsec)',\
    'RC3_r25',\
    'RC3_pa (deg)',\
    'group_dist_std (Mpc)',\
    'groupInfo',\
    'morphology',\
    'distanceIndicator',\
    'luminosityClass',\
    'E(B-V)',\
    'B_median',\
    'B_median_key',\
    'B_max',\
    'B_max_key',\
    'B_min',\
    'B_min_key',\
    'B_sdss',\
    'g_sdss',\
    'r_sdss',
    'z_sdss',\
    'Lstar_median',\
    'Lstar_median_err',\
    'Lstar_max',\
    'Lstar_max_err',\
    'Lstar_min',\
    'Lstar_min_err',\
    'Lstar_sdss',\
    'Lstar_sdss_err',\
    'altNames')
        
    writerOutFile = open(outputFileName,'wb')
    writer = csv.writer(writerOutFile, delimiter="\t", quotechar='', quoting=csv.QUOTE_NONE)
    writer.writerow(fieldnames)
    
    for line in galaxyReader:
        
        degreesJ2000RA,degreesJ2000Dec = eval(line['degreesJ2000RA_Dec'])
        J2000RA,J2000Dec = eval(line['J2000RA_Dec'])
        galacticLong,galacticLat = eval(line['galacticLong_Lat'])
        distvcorr = line['distvcorr (Mpc)']
        RID_mean, RID_std, RID_min, RID_max = eval(line['zIndependentDist_mean_std_min_max (Mpc)'])
        angDiameters = eval(line['angDiameters (arcsec)'])
        angDiameters_err = eval(line['angDiameters_err (arcsec)'])
        linDiameters = eval(line['linDiameters (kpc)'])
        linDiameters_err = eval(line['linDiameters_err (kpc)'])
        diameterKey,ratioKey,paKey = eval(line['diameterKey'])

        # truncate some values because of Python's long floats
        dRA = trunc(degreesJ2000RA,5)
        dDec = trunc(degreesJ2000Dec,5)
#         ra = str(ra[0])+"h"+str(ra[1])+"m"+str(trunc(ra[2],5))+'s'
#         dec = str(dec[0])+'d'+str(dec[1])+'m'+str(trunc(dec[2],5))+'s'
        gLong = trunc(galacticLong,5)
        gLat = trunc(galacticLat,5)
        distvcorr = trunc(distvcorr,2)
#         bestDistance = trunc(bestDistance,3)
#         inc = trunc(inc,2)
#         EBminusV_new = trunc(EBminusV_new,3)
        
        
        objectInfoList = [line['preferredName'],\
        line['NEDname'],\
        line['z'],\
        dRA,dDec,\
        J2000RA,J2000Dec,\
        gLong,gLat,\
        line['helioVelocity (km/s)'],\
        line['vcorr (km/s)'],\
        distvcorr,\
        RID_mean, RID_std, RID_min, RID_max,\
        line['bestDistance (Mpc)'],\
        line['distErr (Mpc)'],\
        angDiameters[0],\
        angDiameters[1],\
        angDiameters_err[0],\
        angDiameters_err[1],\
        linDiameters[0],\
        linDiameters[1],\
        linDiameters_err[0],\
        linDiameters_err[1],\
        line['inc (deg)'],\
        line['adjustedInc (deg)'],\
        line['incErr (deg)'],\
        line['PA (deg)'],\
        diameterKey,\
        ratioKey,\
        paKey,\
        line['RC3_type'],\
        line['RC3_d25 (arcsec)'],\
        line['RC3_r25'],\
        line['RC3_pa (deg)'],\
        line['group_dist_std (Mpc)'],\
        line['groupInfo'],\
        line['morphology'],\
        line['distanceIndicator'],\
        line['luminosityClass'],\
        line['E(B-V)'],\
        line['B_median'],\
        line['B_median_key'],\
        line['B_max'],\
        line['B_max_key'],\
        line['B_min'],\
        line['B_min_key'],\
        line['B_sdss'],\
        line['g_sdss'],\
        line['r_sdss'],\
        line['z_sdss'],\
        line['Lstar_median'],\
        line['Lstar_median_err'],\
        line['Lstar_max'],\
        line['Lstar_max_err'],\
        line['Lstar_min'],\
        line['Lstar_min_err'],\
        line['Lstar_sdss'],\
        line['Lstar_sdss_err'],\
        line['altNames']]
        
        writer.writerow(objectInfoList)
                        
    galaxyFile.close()
    writerOutFile.close()


if __name__=="__main__":
    main()
