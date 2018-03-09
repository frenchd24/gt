#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_finalCSV.py, v2.1 10/19/17

Read in: FinalGalaxyTable2_group2.csv

Separate everything into it's own column and round/trunc numbers as needed

Add a R_vir column
Add RID_median column

Make: FinalGalaxyTable7.csv (9/25/17)

v2: Redo "Name" column moving "MESSIER" farther down the list. Also shorten distIndicator
names as follows:

                   x   =>  x
                  FP   =>  FP      (what is this?)
                 CMD   =>  CMD
                 SBF   =>  SBF
                SNIa   =>  SNIa
                TRGB   =>  TRGB
               Maser   =>  Maser
              Sosies   =>  Sosie   (what is this?)
             D-Sigma   =>  Dsigm
            Cepheids   =>  Ceph
            RR Lyrae   =>  RRLyr
           Tully est   =>  TFest
          SNII radio   =>  SNIIr
        SNII optical   =>  SNIIo
        Tully-Fisher   =>  TF
       GeV TeV ratio   =>  gamma
       Ring Diameter   =>  DRing
       distIndicator   =>  ?????????? (which one?)
     Brightest Stars   =>  BrStr
    Type II Cepheids   =>  CepII
   BL Lac Luminosity   =>  BLLum
   Horizontal Branch   =>  HB
 HII region diameter   =>  dHII
 
 (10/06/17)

v2.1: run on FinalGalaxyTable9_med_groups.csv, rejected_final_combined4_groups.csv
Which now have better distIndicators and proper median RID.
Made: FinalGalaxyTable10.csv, FinalGalaxyTable10_altNames.csv

No longer uses the distIndicator shortener, as this is done in GTupdate2_RIDcombine3.py
instead (10/19/17)

Remake because of the median_low() error (11/03/17)
-> FinalGalaxyTable11.csv, FinalGalaxyTable11_altNames.csv

Remake because of little photometry updates (11/06/17)
-> FinalGalaxyTable11.csv, FinalGalaxyTable11_altNames.csv (duplicate the previous ones)

Remake because of little diameter updates (01/03/18)
-> FinalGalaxyTable12.csv, FinalGalaxyTable12_altNames.csv

Reorder the galaxy names (02/01/18)
-> Make FinalGalaxyTable13.csv, FinalGalaxyTable13_altNames.csv


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
from astropy.table import Table
import pandas as pd

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


def pickPreferredName3(altNames, oldName):
    # oldname is the original name, and possibly not part of altnames
    # altnames is the list of alternative names (probably returned from NED)

#     order = ['NGC','IC','MRK','UGC','PHL','3C','SBS','MCG','ISO','TON','PGC',\
#     'PG','PB','FGC','HS','HE','KUG','IRAS','RX','CGCG','FBQS','LBQS','SDSS','VCC',\
#     '2MASS','2DF','6DF','HIPASS','2MASX']

    order = ['NGC','IC','MRK','UGC','UGCA','PHL','3C','SBS','MCG','ESO','TON','TONS',\
    'PGC','PG','PB','FGC','HS','HE','KUG','IRAS','RX','CGCG','KAZ','FCC','FAIRALL',\
    'HOLM','IZw','IIZw','IIIZw','IVZw','VZw','VIZw','VIIZw','VIIIZw','IRAS','IRASF',\
    'KISS','KISSR','FBQS','LBQS','PKS','SDSS','VCC','2MASS','2DF','6DF','HIPASS','2MASX']
    
    # add oldName to the list of alternate names if it is not already there
    if len(altNames) >=1:
        if not bfind(str(altNames),str(oldName)):
            try:
                altNames.append(oldName)
            except Exception,e:
                sys.stdout.write("Issue with picking preferred name: {0}".format(e))
                sys.stdout.write("altNames list was: {0}".format(altNames))
    else:
        altNames = [oldName]
    
    found = False
    for i in order:
        for n in altNames:
            # n is input string, i is what i'm looking for
            if bfind(n,i) and not found and not bfind(n,':') and not bfind(n,'['):
                finalName = n
                found = True
                break

    if not found:
        finalName = oldName

    return finalName


    
def shorten_diameterkeys(key):
    # returns a shortened version of key for the table. 
    # E.g., 'r (SDSS Isophotal)' becomes r_sdss_iso
    
    d = {'r (SDSS Isophotal)':'r_sdss_iso',\
    'r (SDSS Exponential)':'r_sdss_exp',\
    'r (SDSS de Vaucouleurs)':'r_sdss_dev',\
    'r (SDSS deVaucouleurs)':'r_sdss_dev',\
    'r (SDSS Petrosian)':'r_sdss_pet',\
    'ESO-Uppsala "Quick Blue" IIa-O':'eso_upp',\
    'ESO-LV "Quick Blue" IIa-O':'eso-lv',\
    'K_s (LGA/2MASS "total")':'K_lga2mass_tot',\
    'K_s (2MASS "total")':'K_2mass_tot',\
    'K_s (LGA/2MASS "isophotal")':'K_lga2mass_iso',\
    'K_s (2MASS isophotal)':'K_2mass_iso',\
    'POSS1 103a-E':'poss_103a-E',\
    'POSS1 103a-O':'poss_103a-O',\
    'R (Kron-Cousins)':'R_kron_cousins',\
    'RC3 D_25, R_25 (blue)':'rc3_dr_25',\
    'RC3 D_0 (blue)':'rc3_d0',\
    'x':'x',\
    'None':'x'}
    
    newKey = d[key]
    return newKey


def shorten_distIndicator(key):
    # returns a shortened version of distIndicator for the table
    
    d = {'x':'x',\
    'FP':'FP',\
    'CMD':'CMD',\
    'SBF':'SBF',\
    'SNIa':'SNIa',\
    'TRGB':'TRGB',\
    'Maser':'Maser',\
    'Sosies':'Sosie',\
    'D-Sigma':'Dsigm',\
    'Cepheids':'Ceph',\
    'RR Lyrae':'RRLyr',\
    'Tully est':'TFest',\
    'SNII radio':'SNIIr',\
    'SNII optical':'SNIIo',\
    'Tully-Fisher':'TF',\
    'GeV TeV ratio':'gamma',\
    'Ring Diameter':'DRing',\
    'Brightest Stars':'BrStr',\
    'Type II Cepheids':'CepII',\
    'BL Lac Luminosity':'BLLum',\
    'Horizontal Branch':'HB',\
    'HII region diameter':'dHII'}

    newKey = d[key]
    return newKey
    
    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':        
#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2_group2.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5.csv'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable6_groups.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7.csv'
#         outFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_altNames.csv'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable6_groups_plusRejected.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7.csv'
#         outFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_altNames.csv'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable8_groups_plusRejected.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9.csv'
#         outFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9_altNames.csv'

#         inFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable11_med_groups_plusRejected.csv'
#         outFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable12.csv'
#         outFilename_altNames = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable12_altNames.csv'

        inFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable11_med_groups_plusRejected.csv'
        outFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable13.csv'
        outFilename_altNames = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable13_altNames.csv'
        
    elif user =='David':
#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_med_groups_plusRejected.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11.csv'
#         outFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11_altNames.csv'
        pass

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
        
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
    'Vhel',\
    'vcorr',\
    'distvcorr',\
    'RID_mean',\
    'RID_median',\
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
    'R_vir',\
    'inc',\
    'adjustedInc',\
    'e_inc',\
    'PA',\
    'diameterKey',\
    'ratioKey',\
    'paKey',\
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
    
    fieldnames_altNames = (\
    'Name',\
    'altNames')
        
    # open the file

    nullFloat = -99.99
    nullInt = -99
    nullStr = 'x'


    # open the files
    csv.field_size_limit(sys.maxsize)
    
    # basic info file
    inFile = open(inFilename,'rU')
    reader = csv.DictReader(inFile)

    # open the output file
    writerOutFile = open(outFilename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    # altNames table
    writerOutFile_altNames = open(outFilename_altNames,'wt')
    writer_altNames = csv.DictWriter(writerOutFile_altNames, fieldnames=fieldnames_altNames)
    headers_altNames = dict((n,n) for n in fieldnames_altNames)
    writer_altNames.writerow(headers_altNames)
    
    count = 0
    # loop through and do the work
    for l in reader:
        count +=1
        
        # stuff in basicReader:
        preferredName = l['preferredName']
        NEDname = l['NEDname']
        z = l['z']
        degreesJ2000RA_Dec = l['degreesJ2000RA_Dec']
        J2000RA_Dec = l['J2000RA_Dec']
        galacticLong_Lat = l['galacticLong_Lat']
        Vhel = l['helioVelocity (km/s)']
        vcorr = l['vcorr (km/s)']
        distvcorr = l['distvcorr (Mpc)']
        RID = l['zIndependentDist_mean_median_std_min_max (Mpc)']
        bestDistance = l['bestDistance (Mpc)']
        e_bestDistance = l['distErr (Mpc)']
        angDiameters = l['angDiameters (arcsec)']
        angDiameters_err = l['angDiameters_err (arcsec)']
        linDiameters = l['linDiameters (kpc)']
        linDiameters_err = l['linDiameters_err (kpc)']
        inc = l['inc (deg)']
        adjustedInc = l['adjustedInc (deg)']
        incErr = l['incErr (deg)']
        PA = l['PA (deg)']
        diameterKey = l['diameterKey']
        RC3_type = l['RC3_type']
        RC3_d25 = l['RC3_d25 (arcsec)']
        RC3_r25 = l['RC3_r25']
        RC3_pa = l['RC3_pa (deg)']
        group_num = l['Tully2015_group_number']
        group_mem = l['group_members']
        group_dist = l['groupDist (Mpc)']
        morph = l['morphology']
        distanceIndicator = l['distanceIndicator']
        lumClass = l['luminosityClass']
        EBminusV = l['E(B-V)']
        B_median = l['B_median']
        B_median_key = l['B_median_key']
        B_max = l['B_max']
        B_max_key = l['B_max_key']
        B_min = l['B_min']
        B_min_key = l['B_min_key']
        B_sdss = l['B_sdss']
        g_sdss = l['g_sdss']
        r_sdss = l['r_sdss']
        z_sdss = l['z_sdss']
        Lstar_median = l['Lstar_median']
        Lstar_median_err = l['Lstar_median_err']
        Lstar_max = l['Lstar_max']
        Lstar_max_err = l['Lstar_max_err']
        Lstar_min = l['Lstar_min']
        Lstar_min_err = l['Lstar_min_err']
        Lstar_sdss = l['Lstar_sdss']
        Lstar_sdss_err = l['Lstar_sdss_err']
        altNames = l['altNames']

        # process stuff
        preferredName = NEDname
        Name = pickPreferredName3(eval(altNames),preferredName)
        
#         NEDname 
#         z

        # now adjust some shit
        RAdeg,DEdeg = eval(degreesJ2000RA_Dec)
        
        # separate coordinates
        ra,dec = eval(J2000RA_Dec)
            
        RAh, RAm, RAs = ra[0],ra[1],ra[2]
        
        DEd = int(str(dec[0]).replace('-',''))
        DEm = dec[1]
        DEs = dec[2]
        
        # figure out the sign
        if DEdeg <0:
            DE = '-'
        else:
            DE = '+'
        
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
        RID_mean, RID_median, RID_std, RID_min, RID_max = eval(RID)
        
        if isNumber(RID_mean):
            RID_mean = float(RID_mean)
        else:
            RID_mean = nullFloat
            
        if isNumber(RID_median):
            RID_median = float(RID_median)
        else:
            RID_median = nullFloat
        
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
            R_vir = nullFloat
        else:
            # virial radius
            try:
                R_vir = round(calculateVirialRadius(MajDiam),2)
            except Exception,e:
                print "Couldn't calculate R_vir with MajDiam = {0}".format(MajDiam)
                print "Error was: ",e
                print
                R_vir = nullFloat
        
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
            
        # diameterKeys - split into diameterKey, ratioKey, paKey
        dKeys = eval(diameterKey)
        diameterKey,ratioKey,paKey = dKeys
        diameterKey,ratioKey,paKey = str(diameterKey),str(ratioKey),str(paKey)
        
        diameterKey2 = shorten_diameterkeys(diameterKey)
        ratioKey2 = shorten_diameterkeys(ratioKey)
        paKey2 = shorten_diameterkeys(paKey)
        
        
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
        if not isNumber(group_num):
            group_num = nullInt
            
        if not isNumber(group_mem):
            group_mem = nullInt
        
        if not isNumber(group_dist):
            group_dist = nullFloat
        
        # morph - already using 'x'
        if morph != 'x':
            morph = morph.strip().replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
            morph = morph.replace(' ','_')
        
        # distanceIndicator - already using 'x'
#         distIndicator = shorten_distIndicator(distanceIndicator)
        distIndicator = distanceIndicator
        
        
#         if distanceIndicator != 'x':
#             distanceIndicator = distanceIndicator.strip().replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
#             distanceIndicator = distanceIndicator.replace(' ','_')

        # lumClass - already using 'x'
        if lumClass != 'x':
            lumClass = lumClass.strip().replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
            lumClass = lumClass.replace(' ','_')
        
        # EBminusV
        if not isNumber(EBminusV):
            EBminusV = nullFloat
        
        # B_median
        if not isNumber(B_median):
            B_median = nullFloat
        
        # B_median_key - already using 'x'
        if B_median_key != 'x':
            B_median_key = B_median_key.strip().replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
            B_median_key = B_median_key.replace(' ','_')

        # B_max
        if not isNumber(B_max):
            B_max = nullFloat
        
        # B_max_key - already using 'x'
        if B_max_key != 'x':
            B_max_key = B_max_key.strip().replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
            B_max_key = B_max_key.replace(' ','_')
            
        # B_min
        if not isNumber(B_min):
            B_min = nullFloat
        
        # B_min_key - already using 'x'
        if B_min_key != 'x':
            B_min_key = B_min_key.strip().replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
            B_min_key = B_min_key.replace(' ','_')
        
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
                    
        # back to the level of the full galaxy table - write to file
        outputList = [\
        Name,\
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
        Vhel,\
        vcorr,\
        distvcorr,\
        RID_mean,\
        RID_median,\
        RID_std,\
        RID_min,\
        RID_max,\
        bestDistance,\
        e_bestDistance,\
        MajDiam_ang,\
        MinDiam_ang,\
        e_MajDiam_ang,\
        e_MinDiam_ang,\
        MajDiam,\
        MinDiam,\
        e_MajDiam,\
        e_MinDiam,\
        R_vir,\
        inc,\
        adjustedInc,\
        incErr,\
        PA,\
        diameterKey2,\
        ratioKey2,\
        paKey2,\
        RC3_type,\
        RC3_d25,\
        RC3_r25,\
        RC3_pa,\
        group_num,\
        group_mem,\
        group_dist,\
        morph,\
        distIndicator,\
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
        r_sdss,
        z_sdss,\
        Lstar_median,\
        Lstar_median_err,\
        Lstar_max,\
        Lstar_max_err,\
        Lstar_min,\
        Lstar_min_err,\
        Lstar_sdss,\
        Lstar_sdss_err,\
        newAltNames]
                      
        row = dict((f,o) for f,o in zip(fieldnames,outputList))
        writer.writerow(row)
        
        # now do the altNames
        outputList_altNames = [\
        Name,\
        altNames]
        
        row_altNames = dict((f,o) for f,o in zip(fieldnames_altNames,outputList_altNames))
        writer_altNames.writerow(row_altNames)

        # update counter
        percentComplete = round((float(count)/130759)*100,2)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        sys.stdout.flush()

    inFile.close()
    writerOutFile.close()
    writerOutFile_altNames.close()
    print "Done."

    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
