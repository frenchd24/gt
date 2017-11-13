#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_typeList.py, v 1 07/25/2017

print out the maximum length of every column in FinalGalaxyTable2_group1.csv

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

def choose_altNames(list):
    # pick out NGC, IC, UGC, MRK, SDSS, 2MASS names from the list
    
    newList = []
    haveSDSS = False
    have2MASS = False
    haveNGC = False
    haveUGC = False
    haveIC = False
    haveMRK = False
    for i in list:
        if i[:3] == 'NGC' and not haveNGC:
            newList.append(i)
            haveNGC=True
        elif i[:3] == 'UGC' and not haveUGC:
            newList.append(i)
            haveUGC = True
        elif i[:3] == 'MRK' and not haveMRK:
            newList.append(i)
            haveMRK = True
        elif i[:2] == 'IC' and not haveIC and i[:4] != 'ICRF':
            newList.append(i)
            haveIC = True
        elif i[:4] == 'SDSS' and not haveSDSS:
            newList.append(i)
            haveSDSS = True
        elif i[:4] == '2MAS' and not have2MASS:
            newList.append(i)
            have2MASS = True

    return newList


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
        filename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2_group1.csv'
    
    elif user =='David':
        filename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable2_group1.csv'

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
        preferredName = line['preferredName']
        oldName = line['NEDname']
        z = str(line['z'])
        degreesJ2000RA_Dec = str(line['degreesJ2000RA_Dec'])
        J2000RA_Dec = str(line['J2000RA_Dec'])
        galacticLong_Lat = str(line['galacticLong_Lat'])
        helioVel = str(line['helioVelocity (km/s)'])
        vcorr = str(line['vcorr (km/s)'])
        distvcorr = str(line['distvcorr (Mpc)'])
        RID = str(line['zIndependentDist_mean_std_min_max (Mpc)'])
    	bestDistance = str(line['bestDistance (Mpc)'])
    	distErr = str(line['distErr (Mpc)'])
    	angDiameters = str(line['angDiameters (arcsec)'])
    	angDiameters_err = str(line['angDiameters_err (arcsec)'])
    	linDiameters = str(line['linDiameters (kpc)'])
    	linDiameters_err = str(line['linDiameters_err (kpc)'])
        inc = str(line['inc (deg)'])
    	adjustedInc = str(line['adjustedInc (deg)'])
    	incErr = str(line['incErr (deg)'])
    	PA = str(line['PA (deg)'])
    	diameterKey = str(line['diameterKey'])
        RC3_type = str(line['RC3_type'])
        RC3_d25 = str(line['RC3_d25 (arcsec)'])
        RC3_r25 = str(line['RC3_r25'])
        RC3_pa = str(line['RC3_pa (deg)'])
        Tully2015_group_number = str(line['Tully2015_group_number'])
        group_members = str(line['group_members'])
        groupDist = str(line['groupDist (Mpc)'])
        morph = str(line['morphology'])
        distanceIndicator = str(line['distanceIndicator'])
        lumClass = str(line['luminosityClass'])
        EBminusV = str(line['E(B-V)'])
        B_median = str(line['B_median'])
        B_median_key = str(line['B_median_key'])
        B_max = str(line['B_max'])
        B_max_key = str(line['B_max_key'])
        B_min = str(line['B_min'])
        B_min_key = str(line['B_min_key'])
        B_sdss = str(line['B_sdss'])
        g_sdss = str(line['g_sdss'])
        r_sdss = str(line['r_sdss'])
        z_sdss = str(line['z_sdss'])
        Lstar_median = str(line['Lstar_median'])
        Lstar_median_err = str(line['Lstar_median_err'])
        Lstar_max = str(line['Lstar_max'])
        Lstar_max_err = str(line['Lstar_max_err'])
        Lstar_min = str(line['Lstar_min'])
        Lstar_min_err = str(line['Lstar_min_err'])
        Lstar_sdss = str(line['Lstar_sdss'])
        Lstar_sdss_err = str(line['Lstar_sdss_err'])
        altNames = str(line['altNames'])
        
        if d.has_key('preferredName'):
            if len(d['preferredName']) < len(preferredName):
                d['preferredName'] = preferredName
        else:
            d['preferredName'] = preferredName
            
        if d.has_key('oldName'):
            if len(d['oldName']) < len(oldName):
                d['oldName']=oldName
        else:
            d['oldName']=oldName

        if d.has_key('z'):
            if len(d['z']) < len(z):
                d['z']=z
        else:
            d['z']=z

        if d.has_key('distvcorr'):
            if len(d['distvcorr']) < len(distvcorr):
                d['distvcorr']=distvcorr
        else:
            d['distvcorr']=distvcorr

        if d.has_key('bestDistance'):
            if len(d['bestDistance']) < len(bestDistance):
                d['bestDistance']=bestDistance
        else:
            d['bestDistance']=bestDistance
            
        if d.has_key('distErr'):
            if len(d['distErr']) < len(distErr):
                d['distErr']=distErr
        else:
            d['distErr']=distErr

        Amajor, Aminor = eval(angDiameters)
        if d.has_key('Amajor'):
            if len(d['Amajor']) < len(str(Amajor)):
                d['Amajor']=str(Amajor)
        else:
            d['Amajor']=str(Amajor)

        major, minor = eval(linDiameters)
        if d.has_key('major'):
            if len(d['major']) < len(str(major)):
                d['major']=str(major)
        else:
            d['major']=str(major)
            
        if d.has_key('inc'):
            if len(d['inc']) < len(inc):
                d['inc']=inc
        else:
            d['inc']=inc

        if d.has_key('diameterKey'):
            if len(d['diameterKey']) < len(diameterKey):
                d['diameterKey']=diameterKey
        else:
            d['diameterKey']=diameterKey
            
        if d.has_key('RC3_type'):
            if len(d['RC3_type']) < len(RC3_type):
                d['RC3_type']=RC3_type
        else:
            d['RC3_type']=RC3_type
            
        if d.has_key('RC3_d25'):
            if len(d['RC3_d25']) < len(RC3_d25):
                d['RC3_d25']=RC3_d25
        else:
            d['RC3_d25']=RC3_d25
            
            
        if d.has_key('RC3_r25'):
            if len(d['RC3_r25']) < len(RC3_r25):
                d['RC3_r25']=RC3_r25
        else:
            d['RC3_r25']=RC3_r25
            
        if d.has_key('RC3_pa'):
            if len(d['RC3_pa']) < len(RC3_pa):
                d['RC3_pa']=RC3_pa
        else:
            d['RC3_pa']=RC3_pa
            
        if d.has_key('Tully2015_group_number'):
            if len(d['Tully2015_group_number']) < len(Tully2015_group_number):
                d['Tully2015_group_number']=Tully2015_group_number
        else:
            d['Tully2015_group_number']=Tully2015_group_number
            
        if d.has_key('group_members'):
            if len(d['group_members']) < len(group_members):
                d['group_members']=group_members
        else:
            d['group_members']=group_members
            
        if d.has_key('groupDist'):
            if len(d['groupDist']) < len(groupDist):
                d['groupDist']=groupDist
        else:
            d['groupDist']=groupDist

        if d.has_key('morph'):
            if len(d['morph']) < len(morph):
                d['morph']=morph
        else:
            d['morph']=morph
            
        if d.has_key('lumClass'):
            if len(d['lumClass']) < len(lumClass):
                d['lumClass']=lumClass
        else:
            d['lumClass']=lumClass
    
        if d.has_key('distanceIndicator'):
            if len(d['distanceIndicator']) < len(distanceIndicator):
                d['distanceIndicator']=distanceIndicator
        else:
            d['distanceIndicator']=distanceIndicator
            
        if d.has_key('B_median_key'):
            if len(d['B_median_key']) < len(B_median_key):
                d['B_median_key']=B_median_key
        else:
            d['B_median_key']=B_median_key

        if d.has_key('B_max_key'):
            if len(d['B_max_key']) < len(B_max_key):
                d['B_max_key']=B_max_key
        else:
            d['B_max_key']=B_max_key
            
        if d.has_key('B_min_key'):
            if len(d['B_min_key']) < len(B_min_key):
                d['B_min_key']=B_min_key
        else:
            d['B_min_key']=B_min_key
        
        altNames = str(altNames)
        if d.has_key('altNames'):
            if len(d['altNames']) < len(altNames):
                d['altNames']=altNames
        else:
            d['altNames']=altNames
            
        newAltNames = str(choose_altNames(eval(altNames)))
        
        if d.has_key('newAltNames'):
            if len(d['newAltNames']) < len(newAltNames):
                d['newAltNames']=newAltNames
        else:
            d['newAltNames']=newAltNames
            
    print 'Results: ',
    for i in d:
        print i,' : ',d[i], ' len = ',len(d[i])
        
        
    basicFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
