#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_photometry_rejected.py, v 2.2 09/11/2017

Does the below for the rejected results: rejected_results_redo_phot.csv


Based on:
Id: GTupdate2_photometry.py, v 2.1 04/25/2017

Calculate best B band magnitudes and Lstar estimates from: return_phot_full.csv

combined with distance data from: return_basic_full.csv

Reminder: return_basic_full.csv, return_phot_full.csv, return_pa_full.csv come from:

return111111111111111.csv
returnPhot11111111111111.csv
returnPA11111111111111.csv

and are the result of taking out galaxies with vel or redshift = 'x' from the original
results.

Created processedPhot2.csv on (3/28/17), but have not included E(B-V) yet.

v2.1: fixed to include E(B-V) - made processedPhot2_extinc.csv (4/25/17)

reran to make processedPhot2_extinc3.csv - to insure everything is consistent
with GTupdate2_diameters.py (06/20/17)

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


def calculate_absoluteMag(m,dm,d,dd,e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    #
    # dm is the error in apparent magnitude, dd is the error in distance
    
    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)
    
    # now do the error
    # error is dominated by dd, so don't even bother with an extinction error

    dM = math.sqrt(float(dm)**2 + ((5*float(dd))**2 / (float(d) * math.log10(10))**2))
    
    return M, dM


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
    # reset to 71.0 -> leads to a distance of 1.0 Mpc
    zeroVelocity = 71.
    
    # check which computer we're on, and grab the appropriate file
    user = getpass.getuser()
    if user == 'frenchd':
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return111111111111111.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/returnPA11111111111111.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/returnPhot11111111111111.csv'
        
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_pa_full.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_phot_full.csv'
#         
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedPhot2.csv'

#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_pa_full.csv'
#         photFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_phot_full.csv'
#         
#         outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedPhot2_extinc3.csv'

        basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_extinc_rc3.csv'
        paFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_pa2.csv'
        photFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_phot.csv'
        
        outputFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_processedPhot_extinc.csv'

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
#     paFile = open(paFilename,'rU')
#     paReader = csv.DictReader(paFile)
    
    # Photometry file
    photFile = open(photFilename,'rU')
    photReader = csv.DictReader(photFile)

    
    # the header/column names
    fieldnames = (\
    'oldName',\
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
    'Lstar_sdss_err')
    
    # output file
    writerOutFile = open(outputFilename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    
    # loop through and do the work
    count = -1
    for bline,photline in zip(basicReader,photReader):
    
        # update the counter and print the number
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
        
#         if count >=50:
#             sys.exit()
        
        # stuff in basicReader:
        bpreferredName = bline['preferredName']
        boldName = bline['oldName']
        z = bline['redshift']
        degreesJ2000RA_Dec = bline['degreesJ2000RA_Dec']
        J2000RA_Dec = bline['J2000RA_Dec']
        galacticLong_Lat = bline['galacticLong_Lat']
        RID = eval(bline['rIndependentDistMean_sd_min_max (Mpc)'])
        morph = bline['morphology']
        distanceIndicator = bline['distanceIndicator']
        lumClass = bline['luminosityClass']
        extinc = float(bline['EBminusV'])
        helioVel = bline['radialVelocity (km/s)']
        vcorr = bline['vcorr (km/s)']
        angDiameters = bline['angDiameters (arcsec)']
        linDiameters = bline['linDiameters (kpc)']
        distvcorr = bline['distvcorr (Mpc)']
        inc = bline['inclination (deg)']
        altNames = bline['alternativeNames']
        
##########################################################################################
        # define the best distance
        # use the redshift independent distance if available
        if isNumber(RID[0]):
            bestDistance = float(RID[0])
            
            # use the given error, or assume 10% error
            if isNumber(RID[1]):
                if RID[1] == 0:
                    distErr = bestDistance*0.1
                else:
                    distErr = RID[1]
            else:
                distErr = bestDistance*0.1
        else:
            # otherwise use Hubble's law to define a distance
            if float(vcorr) >0:
                bestDistance = float(vcorr)/hubbleConstant
                distErr = bestDistance*0.1
            else:
                bestDistance = float(zeroVelocity)/hubbleConstant
                distErr = bestDistance*0.2
                        
##########################################################################################
        # stuff in photReader
        # .replace("', ","', '").replace(", '","', '").replace("''","'").replace('"',"'").replace("''","'")
        
        photgalaxyName = photline['galaxyName']
        try:
            photB = eval(photline['B'].replace('masked',"''"))
            photu = eval(photline['u'].replace('masked',"''"))
            photg = eval(photline['g'].replace('masked',"''"))
            photr = eval(photline['r'].replace('masked',"''"))
            photi = eval(photline['i'].replace('masked',"''"))
            photz = eval(photline['z'].replace('masked',"''"))
            photJ = eval(photline['J'].replace('masked',"''"))
            photH = eval(photline['H'].replace('masked',"''"))
            photK = eval(photline['K'].replace('masked',"''"))
        except Exception, e:
            sys.stdout.write("Error: {0}".format(e))
            sys.stdout.write("\n Here was the row: {0}".format(photline))
            sys.exit()
            
            
#         photall = str(photline['all'])
        
        # fix the missing quotes in the lists - fuck it this thing sucks
#         photall = eval(photall.replace("', ","', '").replace(", '","', '").replace("''","'").replace('"',"'").strip())
            

        # set everytyhing to 'x' first in case things don't get set right later
        b_med_val = 'x'
        b_meds = ['x','x','x']
        b_max_val = 'x'
        b_maxes = ['x','x','x']
        b_min_val = 'x'
        b_mins = ['x','x','x']
        sdss_b_val = 'x'
        g_val = 'x'
        r_val = 'x'
        z_val = 'x'
        Lstar_med_val = 'x'
        Lstar_med_err_val = 'x'
        Lstar_max_val = 'x'
        Lstar_max_err_val = 'x'
        Lstar_min_val = 'x'
        Lstar_min_err_val = 'x'
        Lstar_sdss_val = 'x'
        Lstar_sdss_err_val = 'x'
        
        # check to make sure all the rows match. They should always match
        if bpreferredName == photgalaxyName:
            
            # first deal with  magnitudes
            # is there a B band available? If so, choose the median, min, max values
            if len(photB) >=1:
                # b_mins = brightest, b_maxes = dimmest
                b_mins,b_meds,b_maxes = pick_photometry(photB)
            else:
                b_mins,bmeds,bmaxes = ['x','x','x'],['x','x','x'],['x','x','x']
                
            if not isNumber(b_mins[0]) and isNumber(b_meds[0]):
                print
                print '1.5 = ', b_mins,b_meds,b_maxes
                print
                print 'photB: ',photB

            # best SDSS u-values now
            if len(photu) >=1:
                u_val,u_err,u_key = pick_sdss_photometry(photu)
            else:
                u_val,u_err,u_key = 'x','x','x'

            # best SDSS g-values now
            if len(photg) >=1:
                g_val,g_err,g_key = pick_sdss_photometry(photg)
            else:
                g_val,g_err,g_key = 'x','x','x'
                
            # best SDSS r-values now
            if len(photr) >=1:
                r_val,r_err,r_key = pick_sdss_photometry(photr)
            else:
                r_val,r_err,r_key = 'x','x','x'
                
            # best SDSS i-values now
            if len(photi) >=1:
                i_val,i_err,i_key = pick_sdss_photometry(photi)
            else:
                i_val,i_err,i_key = 'x','x','x'
                
            # best SDSS z-values now
            if len(photz) >=1:
                z_val,z_err,z_key = pick_sdss_photometry(photz)
            else:
                z_val,z_err,z_key = 'x','x','x'
                
##########################################################################################
            # now convert to B from SDSS
            if isNumber(g_val) and isNumber(r_val):
                # if there are no errors, set them to 1%
                if not isNumber(g_err):
                    g_err = float(g_val)*0.01
                if not isNumber(r_err):
                    r_err = float(r_val)*0.01
                
                # now do it - using Lupton's 2005 tranformation
                sdss_b, sdss_b_err = convert_sdss_to_B_lupton(g_val,g_err,r_val,r_err)
            
            else:
                sdss_b, sdss_b_err = 'x','x'

##########################################################################################
            # calculate absolute magnitude
            if isNumber(b_mins[0]):
#                 M_b_min, M_b_min_err = calculate_absoluteMag_noExtinc(b_mins[0],b_mins[1],bestDistance,distErr)
#                 M_b_med, M_b_med_err = calculate_absoluteMag_noExtinc(b_meds[0],b_meds[1],bestDistance,distErr)
#                 M_b_max, M_b_max_err = calculate_absoluteMag_noExtinc(b_maxes[0],b_maxes[1],bestDistance,distErr)
                M_b_min, M_b_min_err = calculate_absoluteMag(b_mins[0],b_mins[1],bestDistance,distErr,extinc)
                M_b_med, M_b_med_err = calculate_absoluteMag(b_meds[0],b_meds[1],bestDistance,distErr,extinc)
                M_b_max, M_b_max_err = calculate_absoluteMag(b_maxes[0],b_maxes[1],bestDistance,distErr,extinc)
            else:
                M_b_min, M_b_min_err = 'x','x'
                M_b_med, M_b_med_err = 'x','x'
                M_b_max, M_b_max_err = 'x','x'
                
            # absolute magnitude from SDSS data
            if isNumber(sdss_b):
#                 M_sdss_b, M_sdss_b_err = calculate_absoluteMag_noExtinc(sdss_b,sdss_b_err,bestDistance,distErr)
                M_sdss_b, M_sdss_b_err = calculate_absoluteMag(sdss_b,sdss_b_err,bestDistance,distErr,extinc)
            else:
                M_sdss_b, M_sdss_b_err = 'x','x'
        
            # now calculate Lstar values
            if isNumber(M_b_min):
                Lstar_min, Lstar_min_err = calculate_lstar(averageBstar,M_b_min,M_b_min_err)
                Lstar_med, Lstar_med_err = calculate_lstar(averageBstar,M_b_med,M_b_med_err)
                Lstar_max, Lstar_max_err = calculate_lstar(averageBstar,M_b_max,M_b_max_err)
            else:
                Lstar_min, Lstar_min_err = 'x','x'
                Lstar_med, Lstar_med_err = 'x','x'
                Lstar_max, Lstar_max_err = 'x','x'
            
            if isNumber(M_sdss_b):
                Lstar_sdss, Lstar_sdss_err = calculate_lstar(averageBstar,M_sdss_b,M_sdss_b_err)
            else:
                Lstar_sdss, Lstar_sdss_err = 'x', 'x'

##########################################################################################
            # round things
            if isNumber(b_meds[0]):
                b_med_val = round_to_sig(b_meds[0],sig=3)
            else:
                b_med_val = 'x'
            
            if isNumber(b_maxes[0]):
                b_max_val = round_to_sig(b_maxes[0],sig=3)
            else:
                b_max_val = 'x'
            
            if isNumber(b_mins[0]):
                b_min_val = round_to_sig(b_mins[0],sig=3)
            else:
                b_min_val = 'x'
                
            if isNumber(sdss_b):
                sdss_b_val = round_to_sig(sdss_b,sig=3)
            else:
                sdss_b_val = 'x'
            
            if isNumber(Lstar_med):
                Lstar_med_val = round_to_sig(Lstar_med,sig=3)
                Lstar_med_err_val = round_to_sig(Lstar_med_err,sig=3)
            else:
                Lstar_med_val = 'x'
                Lstar_med_err_val = 'x'

            if isNumber(Lstar_max):
                Lstar_max_val = round_to_sig(Lstar_max,sig=3)
                Lstar_max_err_val = round_to_sig(Lstar_max_err,sig=3)
            else:
                Lstar_max_val = 'x'
                Lstar_max_err_val = 'x'
                
            if isNumber(Lstar_min):
                Lstar_min_val = round_to_sig(Lstar_min,sig=3)
                Lstar_min_err_val = round_to_sig(Lstar_min_err,sig=3)
            else:
                Lstar_min_val = 'x'
                Lstar_min_err_val = 'x'
                
            if isNumber(Lstar_sdss):
                Lstar_sdss_val = round_to_sig(Lstar_sdss,sig=3)
                Lstar_sdss_err_val = round_to_sig(Lstar_sdss_err,sig=3)
            else:
                Lstar_sdss_val = 'x'
                Lstar_sdss_err_val = 'x'

            # now write it all to file
            writeOutList = [\
            bpreferredName,\
            b_med_val,\
            b_meds[2],\
            b_max_val,\
            b_maxes[2],\
            b_min_val,\
            b_mins[2],\
            sdss_b_val,\
            g_val,\
            r_val,\
            z_val,\
            Lstar_med_val,\
            Lstar_med_err_val,\
            Lstar_max_val,\
            Lstar_max_err_val,\
            Lstar_min_val,\
            Lstar_min_err_val,\
            Lstar_sdss_val,\
            Lstar_sdss_err_val]
            
            row = dict((f,o) for f,o in zip(fieldnames,writeOutList))
            writer.writerow(row)
        
        else:
            sys.stdout.write('Something has gone terribly wrong.')
            sys.stdout.write('bpreferredName = {0}, photgalaxyName = {1}'.format(bpreferredName,photgalaxyName))
            sys.exit()
    

    basicFile.close()
    photFile.close()
    writerOutFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
