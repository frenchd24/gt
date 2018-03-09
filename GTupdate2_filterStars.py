#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_filterStars.py, v1.2 01/03/18

v1: Include a new column called 'flag'. flag = 0 for normal galaxies, flag = 1 for 
suspected stars (09/13/17)

v2: also set flag = 2 for objects with:
    abs(RID_median * 71 - vcorr) > 1500  ---> or objects with a weird velocity/distance
    
    makes: FinalGalaxyTable10_filtered.csv
    (10/19/2017)
    
Remade due to distance, diameter, and photometry adjustments (11/06/17)
Makes: FinalGalaxyTable11_filtered

Remade due to diameter adjustments. Added more types to filter out (01/03/18)
Makes: FinalGalaxyTable12_filtered

Reran on FinalGalaxyTable13.csv because of a naming order change (02/01/18)

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

    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':        
#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group_ascii2_1.dat'
#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7.csv'
#         inFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_altNames.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_filtered.csv'
#         outFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7_filtered_altNames.csv'

#         inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9.csv'
#         inFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9_altNames.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9_filtered.csv'
#         outFilename_altNames = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable9_filtered_altNames.csv'

#         inFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable12.csv'
#         inFilename_altNames = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable12_altNames.csv'
#         outFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable12_filtered.csv'
#         outFilename_altNames = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable12_filtered_altNames.csv'

        inFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable13.csv'
        inFilename_altNames = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable13_altNames.csv'
        outFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable13_filtered.csv'
        outFilename_altNames = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable13_filtered_altNames.csv'
        
    elif user =='David':
#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable7.csv'
#         inFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable7_altNames.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable7_filtered.csv'
#         outFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable7_filtered_altNames.csv'

#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9.csv'
#         inFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_altNames.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_filtered.csv'
#         outFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_filtered_altNames.csv'

#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10.csv'
#         inFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_altNames.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_filtered.csv'
#         outFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_filtered_altNames.csv'

#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11.csv'
#         inFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11_altNames.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11_filtered.csv'
#         outFilename_altNames = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable11_filtered_altNames.csv'
        pass
        
    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
        
    # open the files
        # open the files
    csv.field_size_limit(sys.maxsize)
    
    # read in the tables
    inFile = open(inFilename,'rU')
    reader = csv.DictReader(inFile)
        
    inFile_altNames = open(inFilename_altNames,'rU')
    reader_altNames = csv.DictReader(inFile_altNames)
    
    # what are the null values equal to?
    nullFloat = -99.99
    nullInt = -99
    nullStr = 'x'
    
    # hubble constant
    hubbleC = 71.0

    
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
    
    fieldnames_altNames = (\
    'Name',\
    'altNames')
    
            
    # output file
    writerOutFile = open(outFilename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    writerOutFile_altNames = open(outFilename_altNames,'wt')
    writer_altNames = csv.DictWriter(writerOutFile_altNames, fieldnames=fieldnames_altNames)
    headers_altNames = dict((n,n) for n in fieldnames_altNames)
    writer_altNames.writerow(headers_altNames)

#     data = pd.read_csv(inFilename, header=0)
#     print data.head()
    
#     print
#     print 'colnames: ',inFile.colnames
#     print
#     print "inFile['NEDname'][2] = ",inFile['NEDname'][2]

    
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
    
    excludeList = [\
    'M-star',\
    'M_star',\
    'Opt.var.',\
    'K4-K5;Candidate_WD',\
    'F6-F8;Candidate_WD',\
    'A',\
    'Candidate_AGN',\
    'M1',\
    'star??',\
    'O',\
    'K_Star',\
    'PN?',\
    'K1',\
    'M0',\
    'M0V',\
    'A0',\
    'DA-star',\
    'High_vel._cloud',\
    'O',\
    'Carbon',\
    'Point_Src_[SDSS]',\
    'Possible_star',\
    'Planetary_nebula',\
    'M3-M4',\
    'F2',\
    'A-star'\
    'PN:',\
    'Cand._glob._cluster',\
    'Candidate_PN',\
    'F',\
    '0.9',\
    '0.92',\
    '14',\
    '14.247',\
    '14.632',\
    '14.728',\
    '14.818',\
    '14.998',\
    '15.159',\
    '15.171',\
    '15.242',\
    '15.341',\
    '15.458',\
    '15.79',\
    '15.819',\
    '16',\
    '16.281',\
    '16.309',\
    '16.348',\
    '16.394',\
    '16.556',\
    '16.736',\
    '16.764',\
    '16.783',\
    '16.981',\
    '17.012',\
    '17.039',\
    '17.441',\
    '17.597',\
    '2E',\
    '2MASS_Extended_Ver.2',\
    '2_S0_galaxies',\
    '2_S0_pec_galaxies',\
    '2_SB0?_pec_galaxies',\
    '2_Spec?',\
    '2_compacts',\
    '2_or_3?_spirals',\
    '2_spirals',\
    '2_symm.sp.arms',\
    '3_S0_galaxies',\
    ':',\
    'A',\
    'A-star',\
    'A0',\
    'A3_HII',\
    'AGN',\
    'AGN+SF',\
    'AGN1',\
    'AGN2',\
    'AGN:',\
    'AGN?',\
    'ALG',\
    'Amorphous',\
    'B',\
    'B...',\
    'D',\
    'DA',\
    'DA+M',\
    'DA+M:',\
    'DA+M:;_Cand._QSO',\
    'DA-star',\
    'DA',\
    'DA+M',\
    'DA+M:',\
    'DA+M:;_Cand._QSO',\
    'DA-star',\
    'DA:',\
    'DANS',\
    'DANS?',\
    'DANS?_Sbrst',\
    'DANS_WR?',\
    'DA_auto',\
    'DBA',\
    'DC:',\
    'DGTO',\
    'DISRPTD',\
    'DISTRBD',\
    'DQ:',\
    'DQ;_Cand._QSO',\
    'DSa',\
    'K1',\
    'K4-K5;Candidate_WD',\
    'K_Star',\
    'M',\
    'M-star',\
    'M0',\
    'M0V',\
    'M1',\
    'M3-M4',\
    'M_Star',\
    'M_star',\
    'Planetary,_or_galaxy',\
    'Planetary?',\
    'Planetary_nebula',\
    'Possible_*Cl',\
    'Possible_star',\
    'bright_near*',\
    'star:',\
    'star?',\
    'star??',\
    'stellar',\
    'stellar-like',\
    'stellar:']
    
    
##########################################################################################
##########################################################################################    
    # do the work
    count = 0
    for l,lalt in zip(reader,reader_altNames):
        count +=1
        
        Name = l['Name']
        Vhel = l['Vhel']
        vcorr = float(l['vcorr'])
        RID_median = float(l['RID_median'])
        RID_mean = float(l['RID_mean'])
        MType = l['MType'].strip()
        MajDiam = l['MajDiam']
        
        Name_alt = lalt['Name']
        
        flag = 0
        # first check for stars
        if float(Vhel) <= 500:
            if float(MajDiam) == nullFloat:
                if MType == nullStr:
                    flag = 1
                elif MType in excludeList:
                    flag = 1
                else:
                    flag = 0
                    
            elif MType in excludeList:
                flag = 1
        elif MType in excludeList:
            flag = 1
        
        # now check RID and vcorr
        if RID_mean != nullFloat:
            rmeanv = float(RID_mean) * hubbleC
            rmedianv = float(RID_median) * hubbleC
        
            if abs(rmedianv - vcorr) >= 1500:
                print 'Name = ',Name, 'RID_median = {0} vs vcorr/H_0 = {1}'.format(RID_median,vcorr/hubbleC)
                
                if flag ==1:
                    print 'oh shit! flag ==1 already!'
                flag = 2
                
            if abs(rmedianv - vcorr) >= 3000:
                print 'Name = ',Name, 'RID_median = {0} vs vcorr/H_0 = {1}'.format(RID_median,vcorr/hubbleC)
                
                flag = 3
                
##########################################################################################
##########################################################################################
        # now write it to the new file

        outputList = [\
        l['Name'],\
        l['NEDname'],\
        l['z'],\
        l['RAdeg'],\
        l['DEdeg'],\
        l['RAh'],\
        l['RAm'],\
        l['RAs'],\
        l['DE-'],\
        l['DEd'],\
        l['DEm'],\
        l['DEs'],\
        l['GLON'],\
        l['GLAT'],\
        l['Vhel'],\
        l['vcorr'],\
        l['distvcorr'],\
        l['RID_mean'],\
        l['RID_median'],\
        l['RID_std'],\
        l['RID_min'],\
        l['RID_max'],\
        l['bestDist'],\
        l['e_bestDist'],\
        l['MajDiam_ang'],\
        l['MinDiam_ang'],\
        l['e_MajDiam_ang'],\
        l['e_MinDiam_ang'],\
        l['MajDiam'],\
        l['MinDiam'],\
        l['e_MajDiam'],\
        l['e_MinDiam'],\
        l['R_vir'],\
        l['inc'],\
        l['adjustedInc'],\
        l['e_inc'],\
        l['PA'],\
        l['diameterKey'],\
        l['ratioKey'],\
        l['paKey'],\
        l['RC3_type'],\
        l['RC3_d25'],\
        l['RC3_r25'],\
        l['RC3_pa'],\
        l['group_num'],\
        l['group_mem'],\
        l['group_dist'],\
        l['MType'],\
        flag,\
        l['distIndicator'],\
        l['lumClass'],\
        l['E(B-V)'],\
        l['Bmag'],\
        l['Bmag_key'],\
        l['Bmag_max'],\
        l['Bmag_max_key'],\
        l['Bmag_min'],\
        l['Bmag_min_key'],\
        l['Bmag_sdss'],\
        l['gmag_sdss'],\
        l['rmag_sdss'],
        l['zmag_sdss'],\
        l['Lstar_med'],\
        l['e_Lstar_med'],\
        l['Lstar_max'],\
        l['e_Lstar_max'],\
        l['Lstar_min'],\
        l['e_Lstar_min'],\
        l['Lstar_sdss'],\
        l['e_Lstar_sdss'],\
        l['altNames']]
        
        
        outputList_altNames = [lalt['Name'],lalt['altNames']]
                      
        row = dict((f,o) for f,o in zip(fieldnames,outputList))
        writer.writerow(row)
        
        row_altNames = dict((f,o) for f,o in zip(fieldnames_altNames,outputList_altNames))
        writer_altNames.writerow(row_altNames)
        
        
        # update the counter
        percentComplete = round((float(count)/130759)*100,2)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        sys.stdout.flush()
            
    
    # close the files
    inFile.close()
    inFile_altNames.close()
    writerOutFile.close()
    writerOutFile_altNames.close()

    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
