#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_tully2015combine.py, v 1 06/29/2017

Add group info from Tully 2015 to the galaxy table (FinalGalaxyTable2.csv)

Output: FinalGalaxyTable3.csv

ran on rejected_final_combined.csv (09/11/17)

ran on FinalGalaxyTable6.csv -> makes FinalGalaxyTable6_groups.csv

ran on FinalGalaxyTable9_med.csv and rejected_final_combined4.csv to make:
    FinalGalaxyTable9_med_groups.csv, rejected_final_combined4_groups.csv (10/19/17)

ran on FinalGalaxyTable10_med.csv and rejected_final_combined10.csv to make:
    FinalGalaxyTable10_med_groups.csv, rejected_final_combined10_groups.csv (11/03/17)

ran on FinalGalaxyTable11_med.csv and rejected_final_combined11.csv to make:
    FinalGalaxyTable11_med_groups.csv, rejected_final_combined11_groups.csv (01/03/18)

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

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2.csv'
#         groupFilename = '/usr/data/moosejaw/frenchd/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2_group1.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined.csv'
#         groupFilename = '/usr/data/moosejaw/frenchd/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable6.csv'
#         groupFilename = '/usr/data/moosejaw/frenchd/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable6_groups.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined2.csv'
#         groupFilename = '/usr/data/moosejaw/frenchd/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined2_groups.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined3.csv'
#         groupFilename = '/usr/data/moosejaw/frenchd/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined3_groups.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable8.csv'
#         groupFilename = '/usr/data/moosejaw/frenchd/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable8_groups.csv'

#         basicFilename = '/Users/frenchd/Research/GT_update2_files/rejected_final_combined11.csv'
#         groupFilename = '/Users/frenchd/Research/GT_update2_files/tully2015_groups_table5_2.tsv'
#         outputFilename = '/Users/frenchd/Research/GT_update2_files/rejected_final_combined11_groups.csv'
        
        basicFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable11_med.csv'
        groupFilename = '/Users/frenchd/Research/GT_update2_files/tully2015_groups_table5_2.tsv'
        outputFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable11_med_groups.csv'
        

    elif user =='David':
#         basicFilename = '/Users/David/Research_Documents/GT_update2/rejected_final_combined4.csv'
#         groupFilename = '/Users/David/Research_Documents/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/Users/David/Research_Documents/GT_update2/rejected_final_combined4_groups.csv'

#         basicFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_med.csv'
#         groupFilename = '/Users/David/Research_Documents/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_med_groups.csv'

#         basicFilename = '/Users/David/Research_Documents/GT_update2/rejected_final_combined10.csv'
#         groupFilename = '/Users/David/Research_Documents/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/Users/David/Research_Documents/GT_update2/rejected_final_combined10_groups.csv'

#         basicFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_med.csv'
#         groupFilename = '/Users/David/Research_Documents/GT_update2/tully2015_groups_table5_2.tsv'
#         outputFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_med_groups.csv'
        pass
        
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
    fieldnames = (\
    'preferredName',\
    'NEDname',\
    'z',\
    'degreesJ2000RA_Dec',\
    'J2000RA_Dec',\
    'galacticLong_Lat',\
    'helioVelocity (km/s)',\
    'vcorr (km/s)',\
    'distvcorr (Mpc)',\
    'zIndependentDist_mean_median_std_min_max (Mpc)',\
    'bestDistance (Mpc)',\
    'distErr (Mpc)',\
    'angDiameters (arcsec)',\
    'angDiameters_err (arcsec)',\
    'linDiameters (kpc)',\
    'linDiameters_err (kpc)',\
    'inc (deg)',\
    'adjustedInc (deg)',\
    'incErr (deg)',\
    'PA (deg)',\
    'diameterKey',\
    'RC3_type',\
    'RC3_d25 (arcsec)',\
    'RC3_r25',\
    'RC3_pa (deg)',\
    'Tully2015_group_number',\
    'group_members',\
    'groupDist (Mpc)',\
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

    writerOutFile = open(outputFilename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    # group file
    header = ['#','\n','-','PGC',' ; ']
    groupFile = open(groupFilename,'rU')
    groupLines = readlines_into_lists(groupFile, header, ';')
    groupFile.close()
    
    count = 0
    # loop through and do the work
    for bline in basicReader:
        count +=1
        # stuff in basicReader:

        degreesJ2000RA_Dec = bline['degreesJ2000RA_Dec']
        J2000RA_Dec = bline['J2000RA_Dec']
        galacticLong_Lat = bline['galacticLong_Lat']

        altNames = eval(bline['altNames'])
        
        nest = 'x'
        nmb = 'x'
        groupDist = 'x'    
        
        for altName in altNames:
            if bfind(altName,'PGC'):
#                 print 'found one: ',altName
                for gline in groupLines:
                    # order in the group file
                    # PGC;LEDA;GLON;GLAT;SGLON;SGLAT;MType;HV;Vls;Vcmba;J-H;J-K;Kmag;logLKi;logpK;f_Nest;Nest;Nmb;PGC1;Dist;DM;GSGLON;GSGLAT;logLK;CF;sigP;R2t;<Vcmba>;Vbw;e_Vbw;sigbw;sigV;Rbw;e_Rbw;Mvir;Mlum;HDC;LDC;2M++;SGX;SGY;SGZ;Simbad;_RA.icrs;_DE.icrs
                    
                    pgcname = gline[0].strip()
                    
                    if len(pgcname) == 1:
                        newPGCname = 'PGC00000'+pgcname
        
                    elif len(pgcname) == 2:
                        newPGCname = 'PGC0000'+pgcname

                    elif len(pgcname) == 3:
                        newPGCname = 'PGC000'+pgcname
            
                    elif len(pgcname) == 4:
                        newPGCname = 'PGC00'+pgcname
            
                    elif len(pgcname) == 5:
                        newPGCname = 'PGC0'+pgcname
            
                    elif len(pgcname) == 6:
                        newPGCname = 'PGC'+pgcname
                        
                    else:
                        newPGCname = 'PGC'+pgcname
                    
#                     print newPGCname
        
                    if altName == newPGCname:
                        print 'found: ',newPGCname
                        
                        nest = gline[16]
                        nmb = gline[17]
                        groupDist = gline[19]
                
                        break
                    
        # back to the level of the full galaxy table - write to file
        outputList = [\
        bline['preferredName'],\
        bline['NEDname'],\
        bline['z'],\
        bline['degreesJ2000RA_Dec'],\
        bline['J2000RA_Dec'],\
        bline['galacticLong_Lat'],\
        bline['helioVelocity (km/s)'],\
        bline['vcorr (km/s)'],\
        bline['distvcorr (Mpc)'],\
        bline['zIndependentDist_mean_median_std_min_max (Mpc)'],\
        bline['bestDistance (Mpc)'],\
        bline['distErr (Mpc)'],\
        bline['angDiameters (arcsec)'],\
        bline['angDiameters_err (arcsec)'],\
        bline['linDiameters (kpc)'],\
        bline['linDiameters_err (kpc)'],\
        bline['inc (deg)'],\
        bline['adjustedInc (deg)'],\
        bline['incErr (deg)'],\
        bline['PA (deg)'],\
        bline['diameterKey'],\
        bline['RC3_type'],\
        bline['RC3_d25 (arcsec)'],\
        bline['RC3_r25'],\
        bline['RC3_pa (deg)'],\
        nest,\
        nmb,\
        groupDist,\
        bline['morphology'],\
        bline['distanceIndicator'],\
        bline['luminosityClass'],\
        bline['E(B-V)'],\
        bline['B_median'],\
        bline['B_median_key'],\
        bline['B_max'],\
        bline['B_max_key'],\
        bline['B_min'],\
        bline['B_min_key'],\
        bline['B_sdss'],\
        bline['g_sdss'],\
        bline['r_sdss'],\
        bline['z_sdss'],\
        bline['Lstar_median'],\
        bline['Lstar_median_err'],\
        bline['Lstar_max'],\
        bline['Lstar_max_err'],\
        bline['Lstar_min'],\
        bline['Lstar_min_err'],\
        bline['Lstar_sdss'],\
        bline['Lstar_sdss_err'],\
        bline['altNames']]
                      
        row = dict((f,o) for f,o in zip(fieldnames,outputList))
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
