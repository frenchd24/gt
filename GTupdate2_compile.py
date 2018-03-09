#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_compile.py, v2.3 01/03/17

Combine the following files:
return_basic_full_extinc.csv
processedPhot2_extinc.csv
processedDiams.csv

Into one galaxy table. 

Made FinalGalaxyTable2.csv on (6/29/17)


v2: After updating a bunch of stuff, combine these:
return_basic_full_extinc_rc3_rid2.csv
processedPhot3.csv
processedDiams7.csv

(09/25/17)

v2.1: After updating a bunch of stuff, combine these:
return_basic_full_extinc_rc3_rid4.csv
processedPhot3.csv
processedDiams9.csv

Makes: FinalGalaxyTable9_med.csv

(10/19/17)

v2.2: After fixing median_low(), combine these:
return_basic_full_extinc_rc3_rid5.csv
processedPhot3.csv
processedDiams10.csv

Makes: FinalGalaxyTable10_med.csv

and also:
rejected_results_redo_extinc_rc3_rid5.csv
rejected_processedDiams_rid10.csv
rejected_processedPhot2.csv

Makes: rejected_final_combined10.csv
Note - I did this twice, as I updated the photometry to v4 in the meantime.

(11/03/17)


v2.3: After fixing choose_diameter_measurement.py, combine these:
return_basic_full_extinc_rc3_rid5.csv
processedPhot4.csv
processedDiams11.csv

Makes: FinalGalaxyTable11_med.csv

and also:
rejected_results_redo_extinc_rc3_rid5.csv
rejected_processedDiams_rid11.csv
rejected_processedPhot2.csv

Makes: rejected_final_combined11.csv

(01/03/17)


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
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc_rc3.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedDiams.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedPhot2_extinc2.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc_rc3.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedDiams2.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedPhot2_extinc3.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable2.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc_rc3_rid2.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedDiams7.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedPhot3.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable6.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_extinc_rc3_rid.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_processedDiams_rid2.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_processedPhot2.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined2.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_extinc_rc3_rid2.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_processedDiams_rid3.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_processedPhot2.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined3.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc_rc3_rid3.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedDiams8.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedPhot3.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable8.csv'

#         basicFilename = '/Users/frenchd/Research/GT_update2_files/rejected_results_redo_extinc_rc3_rid5.csv'
#         paFilename = '/Users/frenchd/Research/GT_update2_files/rejected_processedDiams_rid11.csv'
#         photFilename = '/Users/frenchd/Research/GT_update2_files/rejected_processedPhot4.csv'
#         outFilename = '/Users/frenchd/Research/GT_update2_files/rejected_final_combined11.csv'

        basicFilename = '/Users/frenchd/Research/GT_update2_files/return_basic_full_extinc_rc3_rid5.csv'
        paFilename = '/Users/frenchd/Research/GT_update2_files/processedDiams11.csv'
        photFilename = '/Users/frenchd/Research/GT_update2_files/processedPhot4.csv'
        outFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable11_med.csv'


        
    elif user =='David':
#         basicFilename = '/Users/David/Research_Documents/GT_update2/rejected_results_redo_extinc_rc3_rid3.csv'
#         paFilename = '/Users/David/Research_Documents/GT_update2/rejected_processedDiams_rid4.csv'
#         photFilename = '/Users/David/Research_Documents/GT_update2/rejected_processedPhot2.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/rejected_final_combined4.csv'

#         basicFilename = '/Users/David/Research_Documents/GT_update2/return_basic_full_extinc_rc3_rid4.csv'
#         paFilename = '/Users/David/Research_Documents/GT_update2/processedDiams9.csv'
#         photFilename = '/Users/David/Research_Documents/GT_update2/processedPhot3.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable9_med.csv'


#         basicFilename = '/Users/David/Research_Documents/GT_update2/rejected_results_redo_extinc_rc3_rid5.csv'
#         paFilename = '/Users/David/Research_Documents/GT_update2/rejected_processedDiams_rid10.csv'
#         photFilename = '/Users/David/Research_Documents/GT_update2/rejected_processedPhot4.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/rejected_final_combined10.csv'

#         basicFilename = '/Users/David/Research_Documents/GT_update2/return_basic_full_extinc_rc3_rid5.csv'
#         paFilename = '/Users/David/Research_Documents/GT_update2/processedDiams10.csv'
#         photFilename = '/Users/David/Research_Documents/GT_update2/processedPhot4.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_med.csv'

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

    # PA and diameter file
    paFile = open(paFilename,'rU')
    paReader = csv.DictReader(paFile)
    
    # Photometry file
    photFile = open(photFilename,'rU')
    photReader = csv.DictReader(photFile)
    
    # the header/column names
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
    
    # output file
    outFile = open(outFilename,'wt')
    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    count = 0
    # loop through and do the work
    for bline,paline,photline in zip(basicReader,paReader,photReader):
        # update the counter
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
        
        # stuff in basicReader:
        bpreferredName = bline['preferredName']
        boldName = bline['oldName']
        z = bline['redshift']
        degreesJ2000RA_Dec = bline['degreesJ2000RA_Dec']
        J2000RA_Dec = bline['J2000RA_Dec']
        galacticLong_Lat = bline['galacticLong_Lat']
        RID_old = bline['rIndependentDistMean_sd_min_max (Mpc)']
        RID = bline['RID']
        morph = bline['morphology']
        distanceIndicator = bline['distanceIndicator']
        lumClass = bline['luminosityClass']
        EBminusV = bline['EBminusV']
        helioVel = bline['radialVelocity (km/s)']
        vcorr = bline['vcorr (km/s)']
#         angDiameters = bline['angDiameters (arcsec)']
#         linDiameters = bline['linDiameters (kpc)']
        distvcorr = bline['distvcorr (Mpc)']
        inc = bline['inclination (deg)']
        altNames = bline['alternativeNames']
        RC3_type = bline['RC3_type']
        RC3_d25 = bline['RC3_d25 (arcsec)']
        RC3_r25 = bline['RC3_r25']
        RC3_pa = bline['RC3_pa (deg)']
    
        # stuff in paReader
        papreferredName = paline['preferredName']
    	vcorr = paline['vcorr (km/s)']
    	bestDistance = paline['bestDistance (Mpc)']
    	distErr = paline['distErr (Mpc)']
    	angDiameters = paline['angDiameters (arcsec)']
    	angDiameters_err = paline['angDiameter errors (arcsec)']
    	linDiameters = paline['linDiameters (kpc)']
    	linDiameters_err = paline['linDiameter errors (kpc)']
    	inc = paline['inc (deg)']
    	adjustedInc = paline['adjustedInc (deg)']
    	incErr = paline['incErr (deg)']
    	PA = paline['PA (deg)']
    	diameterKeys = paline['diameterKeys']
        
        # stuff in photReader
        photgalaxyName = photline['oldName']
        B_median = photline['B_median']
        B_median_key = photline['B_median_key']
        B_max = photline['B_max']
        B_max_key = photline['B_max_key']
        B_min = photline['B_min']
        B_min_key = photline['B_min_key']
        B_sdss = photline['B_sdss']
        g_sdss = photline['g_sdss']
        r_sdss = photline['r_sdss']
        z_sdss = photline['z_sdss']
        Lstar_median = photline['Lstar_median']
        Lstar_median_err = photline['Lstar_median_err']
        Lstar_max = photline['Lstar_max']
        Lstar_max_err = photline['Lstar_max_err']
        Lstar_min = photline['Lstar_min']
        Lstar_min_err = photline['Lstar_min_err']
        Lstar_sdss = photline['Lstar_sdss']
        Lstar_sdss_err = photline['Lstar_sdss_err']
        
        # check to make sure all the rows match. They should always match
        if bpreferredName == papreferredName and bpreferredName == photgalaxyName:
        
            group_dist_std = []
            groupInfo = []

            # now write it to file
            outputList = [bpreferredName,\
            boldName,\
            z,\
            degreesJ2000RA_Dec,\
            J2000RA_Dec,\
            galacticLong_Lat,\
            helioVel,\
            vcorr,\
            distvcorr,\
            RID,\
            bestDistance,\
            distErr,\
            angDiameters,\
            angDiameters_err,\
            linDiameters,\
            linDiameters_err,\
            inc,\
            adjustedInc,\
            incErr,\
            PA,\
            diameterKeys,\
            RC3_type,\
            RC3_d25,\
            RC3_r25,\
            RC3_pa,\
            group_dist_std,\
            groupInfo,\
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
            altNames]

            row = dict((f,o) for f,o in zip(fieldnames,outputList))
            writer.writerow(row)
        
        else:
            sys.stdout.write('Something has gone terribly wrong.')
            sys.stdout.write('bpreferredName = {0}, paoldName = {1}, photgalaxyName = {2}'.format(bpreferredName,paoldName,photgalaxyName))
            sys.exit()
    

    outFile.close()
    basicFile.close()
    paFile.close()
    photFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
