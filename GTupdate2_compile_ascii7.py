#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_compile_ascii7.py, v3 09/26/17

Write FinalGalaxyTable5.dat - change HRV to Vhel, adjust spacing and left-justify text

Starts with FinalGalaxyTable2_group2.csv, which already had the rejected results in it.
(09/13/17)


Based on:
Id: GTupdate2_compile_ascii4_rejected.py, v2.1 09/11/17

Does the below except on the list of rejected galaxies (09/11/17)
 - rejected_final_combined_group.csv
 - produces the following: 
        outFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group_ascii.dat'
        outFilename2 = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group_ascii_altnames.dat'
        outFilename3 = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group_ascii_altnamesLEFT.dat'
 

Based on:
Id: GTupdate2_compile_ascii4.py, v2.1 08/25/17

Made FinalGalaxyTable2.csv on (6/29/17)

Made it work (07/10/2017)

v2: Update on laptop - left justify altnames table using savetxt instead of ascii.table
    (08/25/17)
    
v2.1: Changed name from GTupdate2_compile_ascii_lap.py on (09/11/17) to '4', because this 
makes the '4' version of the tables.

v3: total revamp: use savetxt only, read in from:
    FinalGalaxyTable7.csv and 
    FinalGalaxyTable7_altNames.csv
    
    No longer separate altNames out here. (09/26/17)
    
    -------
    Read in:
        FinalGalaxyTable10_filtered.csv and 
        FinalGalaxyTable10_filtered_altNames.csv
    
    Made:
        FinalGalaxyTable10_filtered.dat
        FinalGalaxyTable10_altNames.dat
    
    (10/19/2017)

    ------
    Read in:
        FinalGalaxyTable11_filtered.csv and 
        FinalGalaxyTable11_filtered_altNames.csv
    
    Made:
        FinalGalaxyTable11_filtered.dat
        FinalGalaxyTable11_altNames.dat
    
    (11/06/2017)
    
'''

import sys
import os
import csv
# import string
import warnings
import numpy as np
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

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group.csv'
#         
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group_ascii2.dat'
#         outFilename2 = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group_ascii_altnames2.dat'
#         outFilename3 = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group_ascii_altnamesLEFT2.dat'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2_group2.csv'
#         
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5.dat'
#         outFilename2 = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5_altnames.dat'
#         outFilename3 = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable5_altNamesLEFT.dat'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7.csv'
#         inFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_altNames.csv'
# 
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7.dat'
#         outFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_altNames.dat'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_filtered.csv'
#         inFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_altNames.csv'
# 
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_filtered.dat'
#         outFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_altNames.dat'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9_filtered.csv'
#         inFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9_altNames.csv'
# 
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9_filtered.dat'
#         outFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9_altNames.dat'
        pass
        
    elif user == 'David':
#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable7_filtered.csv'
#         inFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable7_altNames.csv'
# 
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable7_filtered.dat'
#         outFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable7_altNames.dat'

#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_filtered.csv'
#         inFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_altNames.csv'
# 
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_filtered.dat'
#         outFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_altNames.dat'

#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_filtered.csv'
#         inFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_altNames.csv'
# 
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_filtered.dat'
#         outFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_altNames.dat'

        inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11_filtered.csv'
        inFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11_altNames.csv'

        outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11_filtered.dat'
        outFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11_altNames.dat'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the files
    csv.field_size_limit(sys.maxsize)
    
    # basic info file
    inFile = open(inFilename,'rU')
    reader = csv.DictReader(inFile)
    
    # altNames file
    inFile_altNames = open(inFilename_altNames,'rU')
    reader_altNames = csv.DictReader(inFile_altNames)

    fullPreferred = []
    fullAltNames = []
    
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
    'flag',\
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
    
    fieldnames_altNames = ('Name','altNames')
        
#     'z':'f8',\

    d = {'Name':'U29',\
    'NEDname':'U29',\
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
    'Vhel':'i4',\
    'vcorr':'i4',\
    'distvcorr':'f4',\
    'RID_mean':'f4',\
    'RID_median':'f4',\
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
    'R_vir':'f4',\
    'inc':'i4',\
    'adjustedInc':'i4',\
    'e_inc':'i4',\
    'PA':'i4',\
    'diameterKey':'U14',\
    'ratioKey':'U14',\
    'paKey':'U14',\
    'RC3_type':'U8',\
    'RC3_d25':'f4',\
    'RC3_r25':'f4',\
    'RC3_pa':'i4',\
    'group_num':'i4',\
    'group_mem':'i4',\
    'group_dist':'f4',\
    'MType':'U20',\
    'flag':'i1',\
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
#     table = ascii.Table(list(reader),names=fieldnames, dtype=('f4', 'i4', 'S2'))
    t = table.Table(list(reader), names=fieldnames,dtype=(d['Name'],\
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
    d['Vhel'],\
    d['vcorr'],\
    d['distvcorr'],\
    d['RID_mean'],\
    d['RID_median'],\
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
    d['R_vir'],\
    d['inc'],\
    d['adjustedInc'],\
    d['e_inc'],\
    d['PA'],\
    d['diameterKey'],\
    d['ratioKey'],\
    d['paKey'],\
    d['RC3_type'],\
    d['RC3_d25'],\
    d['RC3_r25'],\
    d['RC3_pa'],\
    d['group_num'],\
    d['group_mem'],\
    d['group_dist'],\
    d['MType'],\
    d['flag'],\
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
    
    
    ascii.write(t, outFilename,format='fixed_width',overwrite=True,delimiter=' ',\
    bookend=False,delimiter_pad=None, formats={'Name':'%-29s',\
    'NEDname':'%-29s',\
    'z':'-.6f',\
    'RAdeg':'-.5f',\
    'DEdeg':'-.5f',\
    'RAh':'-4d',\
    'RAm':'-4d',\
    'RAs':'%-.4f',\
    'DE-':'%-1s',\
    'DEd':'-4d',\
    'DEm':'-4d',\
    'DEs':'-.4f',\
    'GLON':'-.6f',\
    'GLAT':'-.6f',\
    'Vhel':'-4d',\
    'vcorr':'-4d',\
    'distvcorr':'-.2f',\
    'RID_mean':'-.2f',\
    'RID_median':'-.2f',\
    'RID_std':'-.2f',\
    'RID_min':'-.2f',\
    'RID_max':'-.2f',\
    'bestDist':'-.2f',\
    'e_bestDist':'-.2f',\
    'MajDiam_ang':'-.2f',\
    'MinDiam_ang':'-.2f',\
    'e_MajDiam_ang':'-.2f',\
    'e_MinDiam_ang':'-.2f',\
    'MajDiam':'-.2f',\
    'MinDiam':'-.2f',\
    'e_MajDiam':'-.2f',\
    'e_MinDiam':'-.2f',\
    'R_vir':'-.2f',\
    'inc':'-4d',\
    'adjustedInc':'-4d',\
    'e_inc':'-4d',\
    'PA':'-4d',\
    'diameterKey':'%-14s',\
    'ratioKey':'%-14s',\
    'paKey':'%-14s',\
    'RC3_type':'%-8s',\
    'RC3_d25':'-.2f',\
    'RC3_r25':'-.2f',\
    'RC3_pa':'-4d',\
    'group_num':'-6d',\
    'group_mem':'-4d',\
    'group_dist':'-.2f',\
    'MType':'%-20s',\
    'flag':'-1d',\
    'distIndicator':'%-20s',\
    'lumClass':'%-10s',\
    'E(B-V)':'-.4f',\
    'Bmag':'-.2f',\
    'Bmag_key':'%-20s',\
    'Bmag_max':'-.2f',\
    'Bmag_max_key':'%-20s',\
    'Bmag_min':'-.2f',\
    'Bmag_min_key':'%-20s',\
    'Bmag_sdss':'-.2f',\
    'gmag_sdss':'-.2f',\
    'rmag_sdss':'-.2f',
    'zmag_sdss':'-.2f',\
    'Lstar_med':'-.2f',\
    'e_Lstar_med':'-.2f',\
    'Lstar_max':'-.2f',\
    'e_Lstar_max':'-.2f',\
    'Lstar_min':'-.2f',\
    'e_Lstar_min':'-.2f',\
    'Lstar_sdss':'-.2f',\
    'e_Lstar_sdss':'-.2f',\
    'altNames':'%-110s'})
    
    t_altNames = table.Table(list(reader_altNames), names=fieldnames_altNames,dtype=(d['Name'],'U2038'))

    ascii.write(t_altNames, outFilename_altNames,format='fixed_width',overwrite=True,delimiter=' ',\
    bookend=False,delimiter_pad=None,formats={'Name':'%-30s','altNames':'%-2038s'})


#     t2 = table.Table(names=fieldnames2,dtype=(d['Name'],'U2038'))
#     t2 = table.Table(names=fieldnames2,dtype=(d['Name'],'%-2038s'))
    
    
#     fullTable = np.loadtxt(inFile, delimiter=',',dtype={'names': fieldnames,\
#     'formats': ('%-29s',\
#     '%-29s',\
#     '-.6f',\
#     '-.5f',\
#     '-.5f',\
#     '-4d',\
#     '-4d',\
#     '%-.4f',\
#     '%-1s',\
#     '-4d',\
#     '-4d',\
#     '-.4f',\
#     '-.6f',\
#     '-.6f',\
#     '-4d',\
#     '-4d',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-4d',\
#     '-4d',\
#     '-4d',\
#     '-4d',\
#     '%-14s',\
#     '%-14s',\
#     '%-14s',\
#     '%-8s',\
#     '-.2f',\
#     '-.2f',\
#     '-4d',\
#     '-6d',\
#     '-4d',\
#     '-.2f',\
#     '%-20s',\
#     '%-20s',\
#     '%-10s',\
#     '-.4f',\
#     '-.2f',\
#     '%-20s',\
#     '-.2f',\
#     '%-20s',\
#     '-.2f',\
#     '%-20s',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '-.2f',\
#     '%-110s')})

    
#     count = 0
#     # loop through and do the work
#     for line in reader:
#         # update the counter
#         count+=1
#         sys.stdout.write("\r Count = {0}".format(count))
#         sys.stdout.flush()
                
        # stuff in basicReader:
#         preferredName = line['preferredName'].strip()
#         NEDname = line['NEDname'].strip()
#         z = line['z']
#         degreesJ2000RA_Dec = line['degreesJ2000RA_Dec']
#         J2000RA_Dec = line['J2000RA_Dec']
#         galacticLong_Lat = line['galacticLong_Lat']
#         helioVel = line['helioVelocity (km/s)']
#         vcorr = line['vcorr (km/s)']
#         distvcorr = line['distvcorr (Mpc)']
#         RID = line['zIndependentDist_mean_std_min_max (Mpc)']
#         bestDistance = line['bestDistance (Mpc)']
#         distErr = line['distErr (Mpc)']
#         angDiameters = line['angDiameters (arcsec)']
#         angDiameters_err = line['angDiameters_err (arcsec)']
#         linDiameters = line['linDiameters (kpc)']
#         linDiameters_err = line['linDiameters_err (kpc)']
#         inc = line['inc (deg)']
#         adjustedInc = line['adjustedInc (deg)']
#         incErr = line['incErr (deg)']
#         PA = line['PA (deg)']
#         diameterKey = line['diameterKey'].strip()
#         RC3_type = line['RC3_type'].strip()
#         RC3_d25 = line['RC3_d25 (arcsec)']
#         RC3_r25 = line['RC3_r25']
#         RC3_pa = line['RC3_pa (deg)']
#         group_number = line['Tully2015_group_number']
#         group_members = line['group_members']
#         group_dist = line['groupDist (Mpc)']
#         morph = line['morphology'].strip()
#         distanceIndicator = line['distanceIndicator'].strip()
#         lumClass = line['luminosityClass'].strip()
#         EBminusV = line['E(B-V)']
#         B_median = line['B_median']
#         B_median_key = line['B_median_key'].strip()
#         B_max = line['B_max']
#         B_max_key = line['B_max_key'].strip()
#         B_min = line['B_min']
#         B_min_key = line['B_min_key'].strip()
#         B_sdss = line['B_sdss']
#         g_sdss = line['g_sdss']
#         r_sdss = line['r_sdss']
#         z_sdss = line['z_sdss']
#         Lstar_median = line['Lstar_median']
#         Lstar_median_err = line['Lstar_median_err']
#         Lstar_max = line['Lstar_max']
#         Lstar_max_err = line['Lstar_max_err']
#         Lstar_min = line['Lstar_min']
#         Lstar_min_err = line['Lstar_min_err']
#         Lstar_sdss = line['Lstar_sdss']
#         Lstar_sdss_err = line['Lstar_sdss_err']
#         altNames = line['altNames']

        

##########################################################################################
##########################################################################################
        # now write it to the new file
        
#         dat = [preferredName,\
#         NEDname,\
#         z,\
#         RAdeg,\
#         DEdeg,\
#         RAh,\
#         RAm,\
#         RAs,\
#         DE,\
#         DEd,\
#         DEm,\
#         DEs,\
#         GLON,\
#         GLAT,\
#         helioVel,\
#         vcorr,\
#         distvcorr,\
#         RID_mean,\
#         RID_std,\
#         RID_min,\
#         RID_max,\
#         bestDistance,\
#         distErr,\
#         MajDiam_ang,\
#         MinDiam_ang,\
#         e_MajDiam_ang,\
#         e_MinDiam_ang,\
#         MajDiam,\
#         MinDiam,\
#         e_MajDiam,\
#         e_MinDiam,\
#         inc,\
#         adjustedInc,\
#         incErr,\
#         PA,\
#         diameterKey2,\
#         ratioKey2,\
#         paKey2,\
#         RC3_type,\
#         RC3_d25,\
#         RC3_r25,\
#         RC3_pa,\
#         group_number,\
#         group_members,\
#         group_dist,\
#         morph,\
#         distanceIndicator,\
#         lumClass,\
#         EBminusV,\
#         B_median,\
#         B_median_key,\
#         B_max,\
#         B_max_key,\
#         B_min,\
#         B_min_key,\
#         B_sdss,\
#         g_sdss,\
#         r_sdss,\
#         z_sdss,\
#         Lstar_median,\
#         Lstar_median_err,\
#         Lstar_max,\
#         Lstar_max_err,\
#         Lstar_min,\
#         Lstar_min_err,\
#         Lstar_sdss,\
#         Lstar_sdss_err,\
#         str(newAltNames)]
        
        
#         dat2 = [preferredName,altNames]
#         fullPreferred.append(preferredName)
#         fullAltNames.append(str(altNames))
        
        
        # add this row to the table
#         t.add_row(dat)
#         t2.add_row(dat2)
#         
#         if not count % 5000:
            
#             ascii.write(t, outFilename,format='fixed_width',overwrite=True,delimiter=' ',\
#             bookend=False,delimiter_pad=None, formats={'Name':'%-29s',\
#             'NEDname':'%-29s',\
#             'z':'-.6f',\
#             'RAdeg':'-.5f',\
#             'DEdeg':'-.5f',\
#             'RAh':'-4d',\
#             'RAm':'-4d',\
#             'RAs':'%-.4f',\
#             'DE-':'%-1s',\
#             'DEd':'-4d',\
#             'DEm':'-4d',\
#             'DEs':'-.4f',\
#             'GLON':'-.6f',\
#             'GLAT':'-.6f',\
#             'Vhel':'-4d',\
#             'vcorr':'-4d',\
#             'distvcorr':'-.2f',\
#             'RID_mean':'-.2f',\
#             'RID_std':'-.2f',\
#             'RID_min':'-.2f',\
#             'RID_max':'-.2f',\
#             'bestDist':'-.2f',\
#             'e_bestDist':'-.2f',\
#             'MajDiam_ang':'-.2f',\
#             'MinDiam_ang':'-.2f',\
#             'e_MajDiam_ang':'-.2f',\
#             'e_MinDiam_ang':'-.2f',\
#             'MajDiam':'-.2f',\
#             'MinDiam':'-.2f',\
#             'e_MajDiam':'-.2f',\
#             'e_MinDiam':'-.2f',\
#             'inc':'-4d',\
#             'adjustedInc':'-4d',\
#             'e_inc':'-4d',\
#             'PA':'-4d',\
#             'diameterKey':'%-14s',\
#             'ratioKey':'%-14s',\
#             'paKey':'%-14s',\
#             'RC3_type':'%-8s',\
#             'RC3_d25':'-.2f',\
#             'RC3_r25':'-.2f',\
#             'RC3_pa':'-4d',\
#             'group_num':'-6d',\
#             'group_mem':'-4d',\
#             'group_dist':'-.2f',\
#             'MType':'%-20s',\
#             'distIndicator':'%-20s',\
#             'lumClass':'%-10s',\
#             'E(B-V)':'-.4f',\
#             'Bmag':'-.2f',\
#             'Bmag_key':'%-20s',\
#             'Bmag_max':'-.2f',\
#             'Bmag_max_key':'%-20s',\
#             'Bmag_min':'-.2f',\
#             'Bmag_min_key':'%-20s',\
#             'Bmag_sdss':'-.2f',\
#             'gmag_sdss':'-.2f',\
#             'rmag_sdss':'-.2f',
#             'zmag_sdss':'-.2f',\
#             'Lstar_med':'-.2f',\
#             'e_Lstar_med':'-.2f',\
#             'Lstar_max':'-.2f',\
#             'e_Lstar_max':'-.2f',\
#             'Lstar_min':'-.2f',\
#             'e_Lstar_min':'-.2f',\
#             'Lstar_sdss':'-.2f',\
#             'e_Lstar_sdss':'-.2f',\
#             'altNames':'%-110s'})
#             
#             
#             ascii.write(t2, outFilename2,format='fixed_width',overwrite=True,delimiter=' ',\
#             bookend=False,delimiter_pad=None,formats={'Name':'%-30s','altNames':'%-2038s'})

    
    
    # close the files
    inFile.close()
    inFile_altNames.close()
#     outFile.close()
#     outFile2.close()

    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
