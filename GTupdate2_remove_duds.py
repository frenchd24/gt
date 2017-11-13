#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_remove_duds.py, v 1 03/27/2017

Search through the following:
return111111111111111.csv
returnPhot11111111111111.csv
returnPA11111111111111.csv

and remove any targets with velocity or redshift = 'x'

Write out:
return_basic_full.csv
return_pa_full.csv
return_phot_full.csv

rejected_results.csv contains all the results with velocity or redshift ='x'

'''

import sys
import os
import csv
# import string
import warnings
from pylab import *
import numpy
# import atpy
import getpass
from utilities import *
import math


# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################

def median_value(x):
    # returns the value closest to the middle of the sorted list, biased towards the 
    # smaller (i.e., high magnitude) number
    xsorted = sorted(x,reverse=True)
    lx = len(x)
    lx2 = lx/2
    m = xsorted[lx2]
    
    return m



def calculate_absoluteMag_noExtinc(m,dm,d,dd):
    # m is apparent magnitude, d is distance in Mpc
    #
    # dm is the error in apparent magnitude, dd is the error in distance
    
#     print 'm={0}, dm={1}, d={2}, dd={3}'.format(m,dm,d,dd)
    M = float(m) - 5*math.log10((float(d)*10**6)/10)
    
    # now do the error
    dM = math.sqrt(float(dm)**2 + ((5*float(dd))**2 / (float(d) * math.log10(10))**2))
    
    return M,dM


def calculate_absoluteMag(m,d,e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)
    return M


def calculate_lstar(mstar,m,dm):
    # calculate and return L/Lstar
    # mstar is the average M_star value, m is the absolute magnitude in question
    #
    # dm is the error in the absolute magnitude
    
    lratio = 10**(-0.4*(m-mstar))
    
    # now the error
    dLratio = 0.921034 * math.sqrt(10**(-0.8*(m-mstar)) * dm**2)
    
    return lratio,dLratio


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
    
    pass
    
    
def convert_sdss_to_B_jester(g,r):
    # convert SDSS photometry to B-band Johnson following Jester et al. 2005
    
    g = float(g)
    r = float(r)
    b = g + 0.39*(g - r) + 0.21
    return b
    
    
def convert_sdss_to_B_lupton(g,dg,r,dr):
    # convert SDSS photometry to B-band Johnson following Lupton 2005
    # Find it here: http://www.sdss.org/dr12/algorithms/sdssUBVRITransform/#Lupton2005
    #
    # sigma = 0.0107 for this conversion, but the measurement errors are assumed
    # to dominate
    # dg and dr are the errors in the g and r measurements
    
    g = float(g)
    r = float(r)
    dg = float(dg)
    dr = float(dr)
    
    # convert to B
    b = g + 0.3130*(g - r) + 0.2271
    
    # calculate the error
    # the derivatives come out as:
    dbdr = -0.313
    dbdg = 1.313
    
    # now the error formula
    berr = math.sqrt((dbdr*dr)**2 + (dbdg*dg)**2)
    
    return b,berr
    


def pick_photometry(photList):
    # this function picks the best photometry values from the input list of options
    # 
    # input photList should be all of a single photometry band
    #
    # output is either:
    # [min, median, and max] band values, along with the [keys for each] 
    #
    # or:
    # output is [min_val, min_err, min_key], [med_val, med_err, med_key], [max_val, max_err, max_key]

    
    vals = []
    errs = []
    keys = []
    for i in photList:
        key = i[0]
        val = i[1]
        err = i[2]
        unit = i[3]
        
        if unit == 'mag' and isNumber(val):
            # remove the '+/-' from the error value
            try:
                err = str(err).replace('+/-','')
                if isNumber(err):
                    err = float(err)
                else:
                    err = float(val)*0.001
            except Exception, e:
                sys.stdout.write("\n Error: {0}".format(e))
                err = float(val)*0.001
                    
            # add them to the appropriate lists
        
            vals.append(val)
            errs.append(err)
            keys.append(key)

        elif unit == 'microJy':
            # convert to mags
            if isNumber(val):
                if float(val) >0:
                    mag_val = 23.9 - 2.5*math.log10(float(val))
                    
    #                     
    #             except Exception, e:
    #                 sys.stdout.write("Error: {0}".format(e))
    #                 sys.stdout.write("val was = {0}".format(val))
    #                 sys.exit()
            
                    # then remove the '+/-' from the error value
                    try:
                        err = str(err).replace('+/-','')
                        if isNumber(err):
                            mag_err = 23.9 - 2.5*math.log10(float(err))
                        else:
                            mag_err = float(mag_val)*0.001
                    except Exception, e:
                        sys.stdout.write("\n Error: {0}".format(e))
                        mag_err = float(mag_val)*0.001
                
                    # add them to the appropriate lists
    
                    vals.append(mag_val)
                    errs.append(mag_err)
                    keys.append(key)
            
        elif unit == 'Jy':
            # convert to mags
            if isNumber(val):
                if float(val) >0:
                    mag_val = 23.9 - 2.5*math.log10(float(val)*10**6)
            
                    # then remove the '+/-' from the error value
                    try:
                        err = str(err).replace('+/-','')
                        if isNumber(err):
                            mag_err = 23.9 - 2.5*math.log10(float(err)*10**6)
                        else:
                            mag_err = float(mag_val)*0.001
                    except Exception, e:
                        sys.stdout.write("\n Error: {0}".format(e))
                        mag_err = float(mag_val)*0.001
                
                    # add them to the appropriate lists
    
                    vals.append(mag_val)
                    errs.append(mag_err)
                    keys.append(key)
            
    if len(vals) >0:
        # now sort them, largest first
        all = zip(vals,errs,keys)
        all_sorted = sorted(all,reverse=True)

        # unzip them now that they are all sorted
        vals_sorted,errs_sorted,keys_sorted = zip(*all_sorted)
    #     print 'vals_sorted: ',vals_sorted
    #     print 'vals: ',vals

        # grab the median, min and max values
        # MIN MEANS BRIGHTER!
        val_med = median_value(vals_sorted)
        val_min = min(vals_sorted)
        val_max = max(vals_sorted)

    #     print 'val_med: ',val_med
    
        # find the index of the median within the list
        val_med_index = vals_sorted.index(val_med)

        # now use that index to grab the corresponding key and error
        err_med = errs_sorted[val_med_index]
        key_med = keys_sorted[val_med_index]

        # the min and max are easier:
        err_min = errs_sorted[-1]
        err_max = errs_sorted[0]

        key_min = keys_sorted[-1]
        key_max = keys_sorted[0]
    
        results = [val_max,err_max,key_max],[val_med,err_med,key_med],[val_min,err_min,key_min]

    else:
        results = [val_max,err_max,key_max],[val_med,err_med,key_med],[val_min,err_min,key_min] = ['x','x','x'],['x','x','x'],['x','x','x']
    # now return the results
    return results
        
    
def pick_sdss_photometry(photList):
    # this function picks the best sdss photometry values from the input list of options
    # 
    # if there are multiple measurements of a particular type, returns the brightest
    # 
    # input photList should be all of a single photometry band
    #
    # output is [val, err, and key] for best available of cmodel, petrosian, model 
    # or other SDSS magnitudes
    
    vals = []
    errs = []
    keys = []
    
    d = {}
    for i in photList:
        key = i[0]
        val = i[1]
        err = i[2]
        unit = i[3]
        
        if unit == 'mag':
            # remove the '+/-' from the error value
            try:
                err = str(err).replace('+/-','')
                if isNumber(err):
                    err = float(err)
                else:
                    err = float(val)*0.001
            except Exception, e:
                sys.stdout.write("\n Error: {0}".format(e))
                err = float(val)*0.001
            
        # 1st choice is Petrosian
        if bfind(key.lower(),'petrosian'):
            if d.has_key('petrosian'):
                # check if this new measurement is brighter than the last
                # take the brightest (i.e., lowest)
                [val_old,err_old,key_old] = d['petrosian']
                if val < val_old:
                    d['petrosian']=[val,err,key]
            else:
                d['petrosian']=[val,err,key]  
            
        # 2nd choice is model
        elif bfind(key.lower(),'model'):
            if d.has_key('model'):
                # check if this new measurement is brighter than the last
                # take the brightest (i.e., lowest)
                [val_old,err_old,key_old] = d['model']
                if val < val_old:
                    d['model']=[val,err,key]
            else:
                d['model']=[val,err,key]
                
        # 3rd choice is cModel values
        if bfind(key.lower(),'cmodel'):
            if d.has_key('cmodel'):
                # check if this new measurement is brighter than the last
                # take the brightest (i.e., lowest)
                [val_old,err_old,key_old] = d['cmodel']
                if val < val_old:
                    d['cmodel']=[val,err,key]
            else:
                d['cmodel']=[val,err,key]
                
        else:
            if not d.has_key('else'):
                d['else']=[val,err,key]
                
    # best choice
    if d.has_key('cmodel'):
        best = d['cmodel']
        
    # second best choice
    elif d.has_key('petrosian'):
        best = d['petrosian']
        
    # third best choice
    elif d.has_key('model'):
        best = d['model']
        
    # anything SDSS at all
    elif d.has_key('else'):
        best = d['else']
    
    # if something went wrong
    else:
        best = 'x'
              
    # now return the best result
    return best
    
    
def main():
    
    # Mstar to use throughout
    averageBstar = -19.57
    
    # Hubble constant used throughout
    hubbleConstant = 71.0
    
    # what velocity to assign to galaxies with negative or 0 redshifts?
    # initially 110 -> leads to a distance of ~1.6 Mpc, the radius of the local group 
    zeroVelocity = 110
    
    # check which computer we're on, and grab the appropriate file
    user = getpass.getuser()
    if user == 'frenchd':
        basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return111111111111111.csv'
        paFilename = '/usr/data/moosejaw/frenchd/GT_update2/returnPA11111111111111.csv'
        photFilename = '/usr/data/moosejaw/frenchd/GT_update2/returnPhot11111111111111.csv'
        
        outputFilename_basic = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full.csv'
        outputFilename_pa = '/usr/data/moosejaw/frenchd/GT_update2/return_pa_full.csv'
        outputFilename_phot = '/usr/data/moosejaw/frenchd/GT_update2/return_phot_full.csv'
        
        outputFilename_rejected = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results.csv'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the files
    
    # set the field size to really big
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
    fieldnames_basic = ('preferredName',\
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
    
    fieldnames_pa = ('oldName','freqTargeted','diameters','dRatio','pa','completeList')
    fieldnames_phot = ('galaxyName','B','u','g','r','i','z','J','H','K','all')
    
    fieldnames_rejected = ('preferredName',\
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
    'alternativeNames',
    'galaxyName',\
    'B',\
    'u',\
    'g',\
    'r',\
    'i',\
    'z',\
    'J',\
    'H',\
    'K',\
    'all',\
    'oldName',\
    'freqTargeted',\
    'diameters',\
    'dRatio',\
    'pa',\
    'completeList')

    
    # basic output file
    writerOutFile_basic = open(outputFilename_basic,'wt')
    writer_basic = csv.DictWriter(writerOutFile_basic, fieldnames=fieldnames_basic)
    headers_basic = dict((n,n) for n in fieldnames_basic)
    writer_basic.writerow(headers_basic)
    
    # pa output file
    writerOutFile_pa = open(outputFilename_pa,'wt')
    writer_pa = csv.DictWriter(writerOutFile_pa, fieldnames=fieldnames_pa)
    headers_pa = dict((n,n) for n in fieldnames_pa)
    writer_pa.writerow(headers_pa)
    
    # phot output file
    writerOutFile_phot = open(outputFilename_phot,'wt')
    writer_phot = csv.DictWriter(writerOutFile_phot, fieldnames=fieldnames_phot)
    headers_phot = dict((n,n) for n in fieldnames_phot)
    writer_phot.writerow(headers_phot)
    
    # rejected output file
    writerOutFile_rejected = open(outputFilename_rejected,'wt')
    writer_rejected = csv.DictWriter(writerOutFile_rejected, fieldnames=fieldnames_rejected)
    headers_rejected = dict((n,n) for n in fieldnames_rejected)
    writer_rejected.writerow(headers_rejected)
    
    
    # loop through and do the work
    count = -1
    for bline,photline,paline in zip(basicReader,photReader,paReader):
    
        # update the counter and print the number
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
        
        # stuff in basicReader:
        basic_preferredName = bline['preferredName']
        basic_oldName = bline['oldName']
        z = bline['redshift']
        degreesJ2000RA_Dec = bline['degreesJ2000RA_Dec']
        J2000RA_Dec = bline['J2000RA_Dec']
        galacticLong_Lat = bline['galacticLong_Lat']
        RID = eval(bline['rIndependentDistMean_sd_min_max (Mpc)'])
        morphology = bline['morphology']
        distanceIndicator = bline['distanceIndicator']
        lumClass = bline['luminosityClass']
        EBminusV = bline['EBminusV']
        helioVel = bline['radialVelocity (km/s)']
        vcorr = bline['vcorr (km/s)']
        angDiameters = bline['angDiameters (arcsec)']
        linDiameters = bline['linDiameters (kpc)']
        distvcorr = bline['distvcorr (Mpc)']
        inc = bline['inclination (deg)']
        altNames = bline['alternativeNames']
        
        # stuff in photReader
        phot_galaxyName = photline['galaxyName']
        photB = photline['B']
        photu = photline['u']
        photg = photline['g']
        photr = photline['r']
        photi = photline['i']
        photz = photline['z']
        photJ = photline['J']
        photH = photline['H']
        photK = photline['K']
        all = photline['all']
        
        # stuff in paReader
        pa_oldName = paline['oldName']
        freqTargeted = paline['freqTargeted']
        diameters = paline['diameters']
        dRatio = paline['dRatio']
        pa = paline['pa']
        completeList = paline['completeList']
        

##########################################################################################
        # define the best distance
        # use the redshift independent distance if available
        if isNumber(helioVel):

##########################################################################################
            # now write it all back to file
            # basic one
            writeOutList_basic = [\
            basic_preferredName,\
            basic_oldName,\
            z,\
            degreesJ2000RA_Dec,\
            J2000RA_Dec,\
            galacticLong_Lat,\
            RID,\
            morphology,\
            distanceIndicator,\
            lumClass,\
            EBminusV,\
            helioVel,\
            vcorr,\
            angDiameters,\
            linDiameters,\
            distvcorr,\
            inc,\
            altNames]
            
            row_basic = dict((f,o) for f,o in zip(fieldnames_basic,writeOutList_basic))
            writer_basic.writerow(row_basic)
            
            # pa one
            writeOutList_pa = [\
            pa_oldName,\
            freqTargeted,\
            diameters,\
            dRatio,\
            pa,\
            completeList]
            
            row_pa = dict((f,o) for f,o in zip(fieldnames_pa,writeOutList_pa))
            writer_pa.writerow(row_pa)
            
            # phot one
            writeOutList_phot = [\
            phot_galaxyName,\
            photB,\
            photu,\
            photg,\
            photr,\
            photi,\
            photz,\
            photJ,\
            photH,\
            photK,\
            all]
            
            row_phot = dict((f,o) for f,o in zip(fieldnames_phot,writeOutList_phot))
            writer_phot.writerow(row_phot)

        else:
            writeOutList_rejected = [\
            basic_preferredName,\
            basic_oldName,\
            z,\
            degreesJ2000RA_Dec,\
            J2000RA_Dec,\
            galacticLong_Lat,\
            RID,\
            morphology,\
            distanceIndicator,\
            lumClass,\
            EBminusV,\
            helioVel,\
            vcorr,\
            angDiameters,\
            linDiameters,\
            distvcorr,\
            inc,\
            altNames,
            phot_galaxyName,\
            photB,\
            photu,\
            photg,\
            photr,\
            photi,\
            photz,\
            photJ,\
            photH,\
            photK,\
            all,\
            pa_oldName,\
            freqTargeted,\
            diameters,\
            dRatio,\
            pa,\
            completeList]

            row_rejected = dict((f,o) for f,o in zip(fieldnames_rejected,writeOutList_rejected))
            writer_rejected.writerow(row_rejected)
    
    # close the files
    basicFile.close()
    photFile.close()
    paFile.close()
    writerOutFile_basic.close()
    writerOutFile_pa.close()
    writerOutFile_phot.close()
    writerOutFile_rejected.close()
    
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
