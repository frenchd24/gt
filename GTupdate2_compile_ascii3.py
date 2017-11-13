#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_compile_ascii.py, v1 07/10/2017

Made FinalGalaxyTable2.csv on (6/29/17)

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

            
    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc_rc3.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedDiams.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedPhot2_extinc2.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable.csv'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2_group1.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable3.csv'

        inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2_group1.csv'
        
        outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable3_3.dat'
        outFilename2 = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable3_3_altnames.dat'

    elif user == 'David':
        inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable2_group1.csv'
        
        outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable3_3.dat'
        outFilename2 = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable3_3_altnames.dat'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the files
    csv.field_size_limit(sys.maxsize)
    
    # basic info file
    inFile = open(inFilename,'rU')
    reader = csv.DictReader(inFile)

    nullFloat = -99.99
    nullInt = -99
    nullStr = 'x'
    
    # the header/column names
    fieldnames = (\
    'Name',\
    'NEDname',\
    'z',\
    'RAdeg',\
    'DEdeg',\
    'RAh',\
    'RAm',\
    'RAs',\
    'DE-',\
    'DEd',\
    'DEm',\
    'DEs',\
    'GLON',\
    'GLAT',\
    'HRV',\
    'vcorr',\
    'distvcorr',\
    'RID_mean',\
    'RID_std',\
    'RID_min',\
    'RID_max',\
    'bestDist',\
    'e_bestDist',\
    'MajDiam_ang',\
    'MinDiam_ang',\
    'e_MajDiam_ang',\
    'e_MinDiam_ang',\
    'MajDiam',\
    'MinDiam',\
    'e_MajDiam',\
    'e_MinDiam',\
    'inc',\
    'adjustedInc',\
    'e_inc',\
    'PA',\
    'diameterKey',\
    'RC3_type',\
    'RC3_d25',\
    'RC3_r25',\
    'RC3_pa',\
    'group_num',\
    'group_mem',\
    'group_dist',\
    'MType',\
    'distIndicator',\
    'lumClass',\
    'E(B-V)',\
    'Bmag',\
    'Bmag_key',\
    'Bmag_max',\
    'Bmag_max_key',\
    'Bmag_min',\
    'Bmag_min_key',\
    'Bmag_sdss',\
    'gmag_sdss',\
    'rmag_sdss',
    'zmag_sdss',\
    'Lstar_med',\
    'e_Lstar_med',\
    'Lstar_max',\
    'e_Lstar_max',\
    'Lstar_min',\
    'e_Lstar_min',\
    'Lstar_sdss',\
    'e_Lstar_sdss',\
    'altNames')
    
    fieldnames2 = ('Name','altNames')
    
    
    d = {'Name':'U30',\
    'NEDname':'U30',\
    'z':'f8',\
    'RAdeg':'f8',\
    'DEdeg':'f8',\
    'RAh':'i4',\
    'RAm':'i4',\
    'RAs':'f8',\
    'DE-':'U1',\
    'DEd':'i4',\
    'DEm':'i4',\
    'DEs':'f8',\
    'GLON':'f8',\
    'GLAT':'f8',\
    'HRV':'i4',\
    'vcorr':'i4',\
    'distvcorr':'f4',\
    'RID_mean':'f4',\
    'RID_std':'f4',\
    'RID_min':'f4',\
    'RID_max':'f4',\
    'bestDist':'f4',\
    'e_bestDist':'f4',\
    'MajDiam_ang':'f4',\
    'MinDiam_ang':'f4',\
    'e_MajDiam_ang':'f4',\
    'e_MinDiam_ang':'f4',\
    'MajDiam':'f4',\
    'MinDiam':'f4',\
    'e_MajDiam':'f4',\
    'e_MinDiam':'f4',\
    'inc':'i4',\
    'adjustedInc':'i4',\
    'e_inc':'i4',\
    'PA':'i4',\
    'diameterKey':'U102',\
    'RC3_type':'U8',\
    'RC3_d25':'f4',\
    'RC3_r25':'f4',\
    'RC3_pa':'i4',\
    'group_num':'i4',\
    'group_mem':'i4',\
    'group_dist':'f4',\
    'MType':'U20',\
    'distIndicator':'U20',\
    'lumClass':'U10',\
    'E(B-V)':'f4',\
    'Bmag':'f4',\
    'Bmag_key':'U20',\
    'Bmag_max':'f4',\
    'Bmag_max_key':'U20',\
    'Bmag_min':'f4',\
    'Bmag_min_key':'U20',\
    'Bmag_sdss':'f4',\
    'gmag_sdss':'f4',\
    'rmag_sdss':'f4',
    'zmag_sdss':'f4',\
    'Lstar_med':'f4',\
    'e_Lstar_med':'f4',\
    'Lstar_max':'f4',\
    'e_Lstar_max':'f4',\
    'Lstar_min':'f4',\
    'e_Lstar_min':'f4',\
    'Lstar_sdss':'f4',\
    'e_Lstar_sdss':'f4',\
    'altNames':'U110'}
    
    
    # initiate the table
#     table = ascii.Table(names=fieldnames, dtype=('f4', 'i4', 'S2'))
    t = table.Table(names=fieldnames,dtype=(d['Name'],\
    d['NEDname'],\
    d['z'],\
    d['RAdeg'],\
    d['DEdeg'],\
    d['RAh'],\
    d['RAm'],\
    d['RAs'],\
    d['DE-'],\
    d['DEd'],\
    d['DEm'],\
    d['DEs'],\
    d['GLON'],\
    d['GLAT'],\
    d['HRV'],\
    d['vcorr'],\
    d['distvcorr'],\
    d['RID_mean'],\
    d['RID_std'],\
    d['RID_min'],\
    d['RID_max'],\
    d['bestDist'],\
    d['e_bestDist'],\
    d['MajDiam_ang'],\
    d['MinDiam_ang'],\
    d['e_MajDiam_ang'],\
    d['e_MinDiam_ang'],\
    d['MajDiam'],\
    d['MinDiam'],\
    d['e_MajDiam'],\
    d['e_MinDiam'],\
    d['inc'],\
    d['adjustedInc'],\
    d['e_inc'],\
    d['PA'],\
    d['diameterKey'],\
    d['RC3_type'],\
    d['RC3_d25'],\
    d['RC3_r25'],\
    d['RC3_pa'],\
    d['group_num'],\
    d['group_mem'],\
    d['group_dist'],\
    d['MType'],\
    d['distIndicator'],\
    d['lumClass'],\
    d['E(B-V)'],\
    d['Bmag'],\
    d['Bmag_key'],\
    d['Bmag_max'],\
    d['Bmag_max_key'],\
    d['Bmag_min'],\
    d['Bmag_min_key'],\
    d['Bmag_sdss'],\
    d['gmag_sdss'],\
    d['rmag_sdss'],
    d['zmag_sdss'],\
    d['Lstar_med'],\
    d['e_Lstar_med'],\
    d['Lstar_max'],\
    d['e_Lstar_max'],\
    d['Lstar_min'],\
    d['e_Lstar_min'],\
    d['Lstar_sdss'],\
    d['e_Lstar_sdss'],\
    d['altNames']))

    
    t2 = table.Table(names=fieldnames2,dtype=(d['Name'],'U2038'))
    
    
    # output file
#     outFile = open(outFilename,'wt')
#     writer = csv.DictWriter(outFile, fieldnames=fieldnames)
#     headers = dict((n,n) for n in fieldnames)
#     writer.writerow(headers)
    
    count = 0
    # loop through and do the work
    for line in reader:
        # update the counter
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
                
        # stuff in basicReader:
        preferredName = line['preferredName'].strip()
        NEDname = line['NEDname'].strip()
        z = line['z']
        degreesJ2000RA_Dec = line['degreesJ2000RA_Dec']
        J2000RA_Dec = line['J2000RA_Dec']
        galacticLong_Lat = line['galacticLong_Lat']
        helioVel = line['helioVelocity (km/s)']
        vcorr = line['vcorr (km/s)']
        distvcorr = line['distvcorr (Mpc)']
        RID = line['zIndependentDist_mean_std_min_max (Mpc)']
        bestDistance = line['bestDistance (Mpc)']
        distErr = line['distErr (Mpc)']
        angDiameters = line['angDiameters (arcsec)']
        angDiameters_err = line['angDiameters_err (arcsec)']
        linDiameters = line['linDiameters (kpc)']
        linDiameters_err = line['linDiameters_err (kpc)']
        inc = line['inc (deg)']
        adjustedInc = line['adjustedInc (deg)']
        incErr = line['incErr (deg)']
        PA = line['PA (deg)']
        diameterKey = line['diameterKey'].strip()
        RC3_type = line['RC3_type'].strip()
        RC3_d25 = line['RC3_d25 (arcsec)']
        RC3_r25 = line['RC3_r25']
        RC3_pa = line['RC3_pa (deg)']
        group_number = line['Tully2015_group_number']
        group_members = line['group_members']
        group_dist = line['groupDist (Mpc)']
        morph = line['morphology'].strip()
        distanceIndicator = line['distanceIndicator'].strip()
        lumClass = line['luminosityClass'].strip()
        EBminusV = line['E(B-V)']
        B_median = line['B_median']
        B_median_key = line['B_median_key'].strip()
        B_max = line['B_max']
        B_max_key = line['B_max_key'].strip()
        B_min = line['B_min']
        B_min_key = line['B_min_key'].strip()
        B_sdss = line['B_sdss']
        g_sdss = line['g_sdss']
        r_sdss = line['r_sdss']
        z_sdss = line['z_sdss']
        Lstar_median = line['Lstar_median']
        Lstar_median_err = line['Lstar_median_err']
        Lstar_max = line['Lstar_max']
        Lstar_max_err = line['Lstar_max_err']
        Lstar_min = line['Lstar_min']
        Lstar_min_err = line['Lstar_min_err']
        Lstar_sdss = line['Lstar_sdss']
        Lstar_sdss_err = line['Lstar_sdss_err']
        altNames = line['altNames']


        # now adjust some shit
        RAdeg,DEdeg = eval(degreesJ2000RA_Dec)
        
        # separate coordinates
        ra,dec = eval(J2000RA_Dec)
        RAh, RAm, RAs = ra[0],ra[1],ra[2]
        
        # figure out the sign
        DE = str(dec[0])[0]
        
        if isNumber(DE):
            DEd = DE
            DE = '+'
        else:
            DEd = int(str(dec[0])[1:])
            
        DEm = dec[1]
        DEs = dec[2]
        
        # galactic coords
        GLON,GLAT = eval(galacticLong_Lat)
        GLON = round(GLON,6)
        GLAT = round(GLAT,6)

        # vcorr
        if isNumber(vcorr):
            vcorr = int(round(float(vcorr),0))
        else:
            vcorr = nullFloat
        
        # distvcorr
        if isNumber(distvcorr):
            distvcorr = round(float(distvcorr),3)
        else:
            distvcorr = nullFloat
            
        # RID
        RID_mean, RID_std, RID_min, RID_max = eval(RID)
        
        if isNumber(RID_mean):
            RID_mean = float(RID_mean)
        else:
            RID_mean = nullFloat
        
        if isNumber(RID_std):
            RID_std = float(RID_std)
        else:
            RID_std = nullFloat
        
        if isNumber(RID_min):
            RID_min = float(RID_min)
        else:
            RID_min = nullFloat
        
        if isNumber(RID_max):
            RID_max = float(RID_max)
        else:
            RID_max = nullFloat
        
        # angDiameters
        MajDiam_ang, MinDiam_ang = eval(angDiameters)
        e_MajDiam_ang, e_MinDiam_ang = eval(angDiameters_err)
        
        if not isNumber(MajDiam_ang):
            MajDiam_ang = nullFloat
        
        if not isNumber(MinDiam_ang):
            MinDiam_ang = nullFloat
            
        if not isNumber(e_MajDiam_ang):
            e_MajDiam_ang = nullFloat
        
        if not isNumber(e_MinDiam_ang):
            e_MinDiam_ang = nullFloat
            
        
        # linDiameters
        MajDiam,MinDiam = eval(linDiameters)
        e_MajDiam,e_MinDiam = eval(linDiameters_err)
        
        if not isNumber(MajDiam):
            MajDiam = nullFloat
        
        if not isNumber(MinDiam):
            MinDiam = nullFloat
            
        if not isNumber(e_MajDiam):
            e_MajDiam = nullFloat
        
        if not isNumber(e_MinDiam):
            e_MinDiam = nullFloat
        
        # inc
        if not isNumber(inc):
            inc = nullInt
        else:
            inc = int(float(inc))
        
        # adjustedInc
        if not isNumber(adjustedInc):
            adjustedInc = nullInt
        else:
            adjustedInc = int(float(adjustedInc))
            
        # incErr
        if not isNumber(incErr):
            incErr = nullInt
        else:
            incErr = int(float(incErr))
            
        # PA
        if not isNumber(PA):
            PA = nullInt
        else:
            PA = int(float(PA))
            
        # diameterKeys - already using 'x'
#       diameterKeys
        
        # RC3_type - already using 'x'
#         RC3_type
        
        # RC3_d25 - given to 2 decimal places in RC3
        if isNumber(RC3_d25):
            RC3_d25 = round(float(RC3_d25),2)
        else:
            RC3_d25 = nullFloat

        # RC3_r25 - given to 2 decimal places in RC3
        if isNumber(RC3_r25):
            RC3_r25 = round(float(RC3_r25),2)
        else:
            RC3_r25 = nullFloat
        
        # RC3_pa - integer in RC3
        if isNumber(RC3_pa):
            RC3_pa = int(round(float(RC3_pa),0))
        else:
            RC3_pa = nullInt
        
        # Tully2015_group_number
        if not isNumber(group_number):
            group_number = nullInt
            
        if not isNumber(group_members):
            group_members = nullInt
        
        if not isNumber(group_dist):
            group_dist = nullFloat
        
        # morph - already using 'x'
#         morph
        
        # distanceIndicator - already using 'x'
#         distanceIndicator

        # lumClass - already using 'x'
#         lumClass
        
        # EBminusV
        if not isNumber(EBminusV):
            EBminusV = nullFloat
        
        # B_median
        if not isNumber(B_median):
            B_median = nullFloat
        
        # B_median_key - already using 'x'
#         B_median_key

        # B_max
        if not isNumber(B_max):
            B_max = nullFloat
        
        # B_max_key - already using 'x'
#         B_max_key
        
        # B_min
        if not isNumber(B_min):
            B_min = nullFloat
        
        # B_min_key - already using 'x'
#         B_min_key
        
        # B_sdss
        if not isNumber(B_sdss):
            B_sdss = nullFloat
            
        # g_sdss
        if not isNumber(g_sdss):
            g_sdss = nullFloat
            
        # r_sdss
        if not isNumber(r_sdss):
            r_sdss = nullFloat
                
        # z_sdss
        if not isNumber(z_sdss):
            z_sdss = nullFloat
                    
        # Lstar_median
        if not isNumber(Lstar_median):
            Lstar_median = nullFloat
        else:
            Lstar_median = round_to_sig(float(Lstar_median),sig=2)
            
        # Lstar_median_err
        if not isNumber(Lstar_median_err):
            Lstar_median_err = nullFloat
        else:
            Lstar_median_err = round_to_sig(float(Lstar_median_err),sig=2)
        
        # Lstar_max
        if not isNumber(Lstar_max):
            Lstar_max = nullFloat
        else:
            Lstar_max = round_to_sig(float(Lstar_max),sig=2)
            
        # Lstar_max_err
        if not isNumber(Lstar_max_err):
            Lstar_max_err = nullFloat
        else:
            Lstar_max_err = round_to_sig(float(Lstar_max_err),sig=2)

      # Lstar_min
        if not isNumber(Lstar_min):
            Lstar_min = nullFloat
        else:
            Lstar_min = round_to_sig(float(Lstar_min),sig=2)
            
      # Lstar_min_err
        if not isNumber(Lstar_min_err):
            Lstar_min_err = nullFloat
        else:
            Lstar_min_err = round_to_sig(float(Lstar_min_err),sig=2)
            
      # Lstar_sdss
        if not isNumber(Lstar_sdss):
            Lstar_sdss = nullFloat
        else:
            Lstar_sdss = round_to_sig(float(Lstar_sdss),sig=2)
            
      # Lstar_sdss_err
        if not isNumber(Lstar_sdss_err):
            Lstar_sdss_err = nullFloat
        else:
            Lstar_sdss_err = round_to_sig(float(Lstar_sdss_err),sig=2)
            
        # altNames - choose only a few to include, include all in a separate table
        newAltNames = choose_altNames(eval(altNames))
        print 'newAltNames: ',newAltNames
        

##########################################################################################
##########################################################################################
        # now write it to the new file
        
        dat = [preferredName,\
        NEDname,\
        z,\
        RAdeg,\
        DEdeg,\
        RAh,\
        RAm,\
        RAs,\
        DE,\
        DEd,\
        DEm,\
        DEs,\
        GLON,\
        GLAT,\
        helioVel,\
        vcorr,\
        distvcorr,\
        RID_mean,\
        RID_std,\
        RID_min,\
        RID_max,\
        bestDistance,\
        distErr,\
        MajDiam_ang,\
        MinDiam_ang,\
        e_MajDiam_ang,\
        e_MinDiam_ang,\
        MajDiam,\
        MinDiam,\
        e_MajDiam,\
        e_MinDiam,\
        inc,\
        adjustedInc,\
        incErr,\
        PA,\
        diameterKey,\
        RC3_type,\
        RC3_d25,\
        RC3_r25,\
        RC3_pa,\
        group_number,\
        group_members,\
        group_dist,\
        morph,\
        distanceIndicator,\
        lumClass,\
        EBminusV,\
        B_median,\
        B_median_key,\
        B_max,\
        B_max_key,\
        B_min,\
        B_min_key,\
        B_sdss,\
        g_sdss,\
        r_sdss,\
        z_sdss,\
        Lstar_median,\
        Lstar_median_err,\
        Lstar_max,\
        Lstar_max_err,\
        Lstar_min,\
        Lstar_min_err,\
        Lstar_sdss,\
        Lstar_sdss_err,\
        str(newAltNames)]
        
#         print 'dat: ',dat
        
        dat2 = [preferredName,altNames]
    
        # add this row to the table
        t.add_row(dat)
        t2.add_row(dat2)
        
        if not count % 1000:
            # write it to file
#             ascii.write(t, outFilename,format='fixed_width', bookend=False,\
#             delimiter=' ',overwrite=True)
            ascii.write(t, outFilename,format='fixed_width',overwrite=True,delimiter='',bookend=False)
            ascii.write(t2, outFilename2,format='fixed_width',overwrite=True,delimiter='',bookend=False)
         
    
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
