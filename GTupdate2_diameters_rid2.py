#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_diameters_rid2.py, v 2.4 01/03/17

Calculate diameters by normalizing to 2MASS K_s "Total" values:
Uses:

return_pa_full.csv

combined with distance data from: return_basic_full.csv

Reminder: return_basic_full.csv, return_phot_full.csv, return_pa_full.csv come from:

return111111111111111.csv
returnPhot11111111111111.csv
returnPA11111111111111.csv

and are the result of taking out galaxies with vel or redshift = 'x' from the original
results.

Made processedDiams.csv (5/23/17) and processedDiams2.csv on (6/21/17)

v1.1 reran on (09/19/17) using the updated choose_diameter_measurement3.py module
    made processedDiams3.csv 
    
v2: reran and adapted for use with return_basic_full_extinc_rc3_rid2.csv and the new
    redshift independent distances therein (09/22/17)
    made processedDiams4.csv, processedDiams5.csv, processedDiams6.csv
    
v2.1: fixed a few errors - ESO-Uppsala and ESO-LV had the same fits accidentally
    made processedDiams7.csv, rejected_processedDiams_rid2.csv (09/25/17)
    
    made processedDiams8.csv and rejected_processedDiams_rid3.csv sometime after this
    
v2.2: reran and adapted for use with return_basic_full_extinc_rc3_rid4.csv and the new
    redshift independent distances therein (10/19/17)
    made: rejected_processedDiams_rid4.csv, processedDiams9.csv
    
v2.3: reran for use with return_basic_full_extinc_rc3_rid5.csv. Fixed a typo in defining
    bestDistance (11/02/17)
    made: rejected_processedDiams_rid10.csv, processedDiams10.csv
    
v2.4: reran with drastically relaxed majorThresh (= 14.0) because some 2MASS_total diams
    were getting filtered out erroneously (01/03/17)
    made: rejected_processedDiams_rid11.csv, processedDiams11.csv
    
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
from choose_diameter_measurement4 import choose_diameter_measurement

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################

def absoluteMag_noExtinc(m,d):
    # m is apparent magnitude, d is distance in Mpc
    M = float(m) - 5*math.log10((float(d)*10**6)/10)
    return M


def absoluteMag(m,d,e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)
    return M


def return_best_fit(survey,major,ratio,bestDistance,distErr):
    # given a survey name and major and minor axis, returns the best fit 2MASS converted
    # value
    
    def func(p, x):
        b, c = p
        return b*x + c
    
    d ={\
    'K_s (2MASS "total")':[[1,0],[0,0]],\
    'K_s (LGA/2MASS "total")':[[1,0],[0,0]],\
    'K_s (2MASS isophotal)':[[1.76504577,1.31252501],[0.00279687,0.06004326]],\
    'K_s (LGA/2MASS "isophotal")':[[1.76504577,1.31252501],[0.00279687,0.06004326]],\
    'POSS1 103a-O':[[0.86894407,17.60096996],[0.00722472,0.34873717]],\
    'POSS1 103a-E':[[1.05459186,26.22334691],[0.04321869,1.98116379]],\
    'r (SDSS Isophotal)':[[1.03338444,0.84324824],[0.00537146,0.1721808]],\
    'ESO-LV "Quick Blue" IIa-O':[[0.8062246,-9.72598723],[0.01547763,1.3711384]],\
    'RC3 D_0 (blue)':[[1.03983059,-1.28989475],[0.00869384,0.58225039]],\
    'RC3 D_25, R_25 (blue)':[[1.10673687,-3.08535534],[0.00924231,0.59831013]],\
    'r (SDSS Petrosian)':[[4.7280691,3.37797516],[0.03181351,0.20621408]],\
    'r (SDSS deVaucouleurs)':[[3.49110412,14.45098353],[0.0389602,0.20565904]],\
    'r (SDSS de Vaucouleurs)':[[2.69532537,15.64423927],[0.03605362,0.22265487]],\
    'RC3 A_e (Johnson B)':[[2.24618406,19.19528334],[0.05678933,1.78234475]],\
    'r (SDSS Exponential)':[[8.23675224,7.72878368],[0.0632389,0.19645029]],\
    'B (Johnson)':[[1.23514983,-25.31353422],[0.09173294,9.72601421]],\
    'R (Kron-Cousins)':[[1.46643311, -35.89213248],[0.13853774, 14.4131571]],\
    'ESO-Uppsala "Quick Blue" IIa-O':[[1.06022324, -13.38635197],[0.01899924, 1.37611021]],\
    'ESO-LV IIIa-F':[[4.6539302,23.02358588],[0.12613219,0.9216136]]
    }
    
    if d.has_key(survey):
        # get the fits and errors
        fit,errs = d[survey]
        
        minor = float(major)*float(ratio)
        
        # apply the fits
#         newMajor,newMinor = func(fit,major),func(fit,minor)
        newMajor = func(fit,major)
        newMinor = newMajor*float(ratio)

        # if something goes wrong, default to zero
        if newMinor <= 0:
            newMajor = func(fit,0.0)
            newMinor = newMajor * float(ratio)
    
        # add and subtract the errors
        upErr = np.array(fit)+np.array(errs)
        downErr = np.array(fit)-np.array(errs)

        # get the upper and lower bounds
        majorUpErr,minorUpErr = func(upErr,major),func(upErr,minor)
        majorDownErr,minorDownErr = func(downErr,major),func(downErr,minor)
        
        # take the larger error (up vs down) as the single error value
        if majorUpErr - newMajor >= newMajor - majorDownErr:
            finalMajorErr = majorUpErr - newMajor
        else:
            finalMajorErr = newMajor - majorDownErr
        
        if minorUpErr - newMinor >= newMinor - minorDownErr:
            finalMinorErr = minorUpErr - newMinor
        else:
            finalMinorErr = newMinor - minorDownErr
        
        # best linear diameters
        linMajor,linMinor = calculateLinearDiameters(newMajor,newMinor,bestDistance)
        
        # maximum linear diameters with error
        # could be as big as:
        linMajorErrUp,linMinorErrUp = calculateLinearDiameters(newMajor+finalMajorErr,newMinor+finalMinorErr,bestDistance+distErr)
        
        # or as small as:
        linMajorErrDown,linMinorErrDown = calculateLinearDiameters(newMajor-finalMajorErr,newMinor-finalMinorErr,bestDistance-distErr)
        
        print
        print 'linMajor = {0} + {1}'.format(linMajor,linMajorErrUp)
        print 'newMajor= {0}, finalMajorErr = {1}, bestDistance= {2}, distErr={3}'.format(newMajor,finalMajorErr,bestDistance,distErr)


        # just take the biggest:
        if linMajorErrUp - linMajor >= linMajor - linMajorErrDown:
            finalLinMajorErr = linMajorErrUp - linMajor
        else:
            finalLinMajorErr = linMajor - linMajorErrDown
        
        if linMinorErrUp - linMinor >= linMinor - linMinorErrDown:
            finalLinMinorErr = linMinorErrUp - linMinor
        else:
            finalLinMinorErr = linMinor - linMinorErrDown
            
#         finalMajorErr = round_to_sig(finalMajorErr,sig=errSig)
#         finalMinorErr = round_to_sig(finalMinorErr,sig=errSig)
#         finalLinMajorErr = round_to_sig(finalLinMajorErr,sig=errSig)
#         finalLinMinorErr = round_to_sig(finalLinMinorErr,sig=errSig)

        
    else:
        print 'survey not found: ',survey
        newMajor,newMinor = 'x','x'
        finalMajorErr,finalMinorErr = 'x','x'
        
        linMajor,linMinor = 'x','x'
        finalLinMajorErr,finalLinMinorErr = 'x','x'
        
        
    return [newMajor,newMinor],[finalMajorErr,finalMinorErr],[linMajor,linMinor],[finalLinMajorErr,finalLinMinorErr]
    
    
def calcIncError(major,minor,dmaj,dmin):
    major = float(major)
    minor = float(minor)
    dmaj = float(dmaj)
    dmin = float(dmin)
    
    if minor != major:
        # minor term
        minErrTerm = (180 * dmin)**2 / ((1 - (minor/major)**2) * major**2 * math.pi**2)

        # major term
        majErrTerm = (180 * dmaj * minor)**2 / ((1 - (minor/major)**2) * major**4 * math.pi**2)
        

    else:
        minErrTerm = 0.0
        majErrTerm = 0.0
        
    incError = math.sqrt(majErrTerm + minErrTerm)
    return incError
    
    
def main():

    # define a 'zeroVelocity' to use if vcorr<=0
    zeroVelocity = 71.
    
    # hubble constant to use throughout
    hubbleConstant = 71.
    
    # sig figs for diameter rounding
    diamSig = 4
    
    # sig figs for inclination rounding
    incSig = 2
    
    # sig figs for diameter errors
    errSig = 2
    
    # check which computer we're on, and grab the appropriate file
    user = getpass.getuser()

    if user == 'frenchd':
        # the full table
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc_rc3_rid2.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_pa_full.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedDiams7.csv'

        # the initially rejected set
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_extinc_rc3_rid.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_pa2.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_processedDiams_rid2.csv'

        # the full table, again
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_basic_full_extinc_rc3_rid3.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/return_pa_full2.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/processedDiams8.csv'

        # the initially rejected set, again
#         basicFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_results_redo_extinc_rc3_rid2.csv'
#         paFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_PA2.csv'
#         outFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_processedDiams_rid3.csv'

        # the full table, again
#         basicFilename = '/Users/frenchd/Research/GT_update2_files/return_basic_full_extinc_rc3_rid5.csv'
#         paFilename = '/Users/frenchd/Research/GT_update2_files/return_pa_full2.csv'
#         outFilename = '/Users/frenchd/Research/GT_update2_files/processedDiams11.csv'

        # the initially rejected set, again
        basicFilename = '/Users/frenchd/Research/GT_update2_files/rejected_results_redo_extinc_rc3_rid5.csv'
        paFilename = '/Users/frenchd/Research/GT_update2_files/rejected_PA2.csv'
        outFilename = '/Users/frenchd/Research/GT_update2_files/rejected_processedDiams_rid11.csv'


    elif user == 'David':
#         basicFilename = '/Users/David/Research_Documents/GT_update2/return_basic_full_extinc_rc3_rid3.csv'
#         paFilename = '/Users/David/Research_Documents/GT_update2/return_pa_full2.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/processedDiams9.csv'

#         basicFilename = '/Users/David/Research_Documents/GT_update2/rejected_results_redo_extinc_rc3_rid3.csv'
#         paFilename = '/Users/David/Research_Documents/GT_update2/rejected_PA2.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/rejected_processedDiams_rid4.csv'
        
        
#         basicFilename = '/Users/David/Research_Documents/GT_update2/return_basic_full_extinc_rc3_rid5.csv'
#         paFilename = '/Users/David/Research_Documents/GT_update2/return_pa_full2.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/processedDiams10.csv'
        
        pass
        
#         basicFilename = '/Users/David/Research_Documents/GT_update2/rejected_results_redo_extinc_rc3_rid5.csv'
#         paFilename = '/Users/David/Research_Documents/GT_update2/rejected_PA2.csv'
#         outFilename = '/Users/David/Research_Documents/GT_update2/rejected_processedDiams_rid10.csv'

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
    
    # the header/column names
    fieldnames = (\
    'preferredName',\
    'vcorr (km/s)',\
    'bestDistance (Mpc)',\
    'distErr (Mpc)',\
    'angDiameters (arcsec)',\
    'angDiameter errors (arcsec)',\
    'linDiameters (kpc)',\
    'linDiameter errors (kpc)',\
    'inc (deg)',\
    'adjustedInc (deg)',\
    'incErr (deg)',\
    'PA (deg)',\
    'diameterKeys')
    
    # output file
    outFile = open(outFilename,'wt')
    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    count = -1
    
    # loop through and do the work
    for bline,paline in zip(basicReader,paReader):
    
        # update the counter and print the number
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
        
        # stuff in basicReader:
        bpreferredName = bline['preferredName']
        boldName = bline['oldName']
        z = bline['redshift']
        degreesJ2000RA_Dec = eval(bline['degreesJ2000RA_Dec'])
        J2000RA_Dec = bline['J2000RA_Dec']
        galacticLong_Lat = bline['galacticLong_Lat']
        RID = eval(bline['rIndependentDistMean_sd_min_max (Mpc)'])
        RID2 = eval(bline['RID'])
        morph = bline['morphology']
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
    
        # stuff in paReader
        paoldName = paline['oldName']
        freqTargeted = paline['freqTargeted']
        paDiameters = eval(paline['diameters'])
        dRatio = paline['dRatio']
        pa = paline['pa']
        completeList = paline['completeList']
        
##########################################################################################
        # define the best distance
        # use the redshift independent distance if available
    
        if isNumber(RID2[1]):
            bestDistance = float(RID2[1])
            
            # use the given error, or assume 10% error
            if isNumber(RID2[2]):
                if RID2[2] == 0:
                    distErr = bestDistance*0.1
                else:
                    distErr = float(RID2[2])
            else:
                distErr = bestDistance*0.1
        else:
            # otherwise use Hubble's law to define a distance
            if float(vcorr) >0:
                bestDistance = float(vcorr)/hubbleConstant
                distErr = bestDistance*0.1
            else:
                bestDistance = float(zeroVelocity)/hubbleConstant
                distErr = bestDistance*0.5
                
##########################################################################################

        # check to make sure all the rows match. They should always match
        if bpreferredName == paoldName:
            major, minor = paDiameters
            
#             if isNumber(major) and isNumber(minor):

            # choose which measurement to use
            bestDiameters, bestKeys = choose_diameter_measurement(completeList)
                
            # the first key is for the diameter fits, the second for the PA
            diameterKey,ratioKey, paKey = bestKeys
            
            # grab the ratio
            bestRatio = bestDiameters[2]

            # grab the PA
            bestPA = bestDiameters[3]
            
            # switch the order of major to minor if they are backwards
            bestMajor,bestMinor = bestDiameters[1]
            
            if isNumber(bestMajor) and isNumber(bestMinor):
                if bestMinor > bestMajor:
                    bestMajor,bestMinor = bestMinor,bestMajor
            
                # make fits for new diameters
                newDiameters,errors,linDiameters,linErrors = return_best_fit(diameterKey,bestMajor,bestRatio,bestDistance,distErr)
#                 print 'errors: ',errors
                
                newDiameters = (round_to_sig(newDiameters[0],sig=diamSig),round_to_sig(newDiameters[1],sig=diamSig))
                linDiameters = (round_to_sig(linDiameters[0],sig=diamSig),round_to_sig(linDiameters[1],sig=diamSig))

                # round the errors
                if errors[0]>0.0 and errors[1]>0.0:
                    print 'inside error rounding',errors,linErrors
                    errors = (round_to_sig(errors[0],sig=errSig),round_to_sig(errors[1],sig=errSig))
                else:
                    errors = (round_to_sig(newDiameters[0]*0.05,sig=errSig),round_to_sig(newDiameters[1]*0.05,sig=errSig))
                
                if linErrors[0] >0.0 and linErrors[1]>0.0:
                    linErrors = (round_to_sig(linErrors[0],sig=errSig),round_to_sig(linErrors[1],sig=errSig))
                else:
                    linErrors = (round_to_sig(linDiameters[0]*0.05,sig=errSig),round_to_sig(linDiameters[1]*0.05,sig=errSig))
        
                if isNumber(major) and not isNumber(newDiameters[0]):
                    print 'problem: major = {0} and new = {1}'.format(major,newDiameters[0])
                
#                 if count >12211:
#                     print
#                     print 'newDiameters: ',newDiameters
#                     print 'bestDiameters: ',bestDiameters
                    
                
                # calculate inclination
                if isNumber(newDiameters[0]):
                    inc = calculateInclination(newDiameters[0],newDiameters[1])
                    incErr = calcIncError(newDiameters[0],newDiameters[1],errors[0],errors[1])
                    
                    if inc>0:
                        inc = round_to_sig(inc,sig=incSig)
                        incErr = round_to_sig(incErr,sig=incSig-1)
        
                    # calculate adjusted inclination with q0=0.2
                    q0 = 0.2
                    adjustedInc = calculateFancyInclination(newDiameters[0],newDiameters[1],q0)
                    
                    if adjustedInc >0:
                        adjustedInc = round_to_sig(adjustedInc,sig=incSig)
                else:
                    inc = 'x'
                    adjustedInc = 'x'
                    incErr = 'x'
                    
            elif isNumber(major):
                # make fits for new diameters - set ratio = 1.0 when the minor is not available, discard later
                newDiameters,errors,linDiameters,linErrors = return_best_fit(diameterKey,bestMajor,1.0,bestDistance,distErr)
        
                if isNumber(major) and not isNumber(newDiameters[0]):
                    print 'problem: major = {0} and new = {1}'.format(major,newDiameters[0])
                    
                newDiameters = (round_to_sig(newDiameters[0],sig=diamSig),'x')
                linDiameters = (round_to_sig(linDiameters[0],sig=diamSig),'x')
                
                if errors[0]>0.0:
                    errors = (round_to_sig(errors[0],sig=errSig),'x')
                    linErrors = (round_to_sig(linErrors[0],sig=errSig),'x')
                    
            else:
                # set things to 'x' if there is no data
                newDiameters = ('x','x')
                errors = ('x','x')
                linDiameters = ('x','x')
                linErrors = ('x','x')
                inc = 'x'
                adjustedInc = 'x'
                incErr = 'x'
                diameterKey = 'x'
                bestPA = 'x'
                bestKeys = ['x','x','x']
                
                # print something out if this shit happens
                if isNumber(major) and not isNumber(minor):
                    print 'WTF: major = {0}, minor = {1}'.format(major,minor)
                    
                    
##########################################################################################
            if not isNumber(bestPA):
                bestPA = 'x'

            # now write it to file
            outputList = [\
            bpreferredName,\
            trunc(vcorr,2),\
            trunc(bestDistance,2),\
            trunc(distErr,2),\
            newDiameters,\
            errors,\
            linDiameters,\
            linErrors,\
            inc,\
            adjustedInc,\
            incErr,\
            bestPA,\
            bestKeys]

            row = dict((f,o) for f,o in zip(fieldnames,outputList))
            writer.writerow(row)
            
        else:
            sys.stdout.write('Something has gone terribly wrong.')
            sys.stdout.write('bpreferredName = {0}, paoldName = {1}'.format(bpreferredName,paoldName))
            sys.exit()
    

    paFile.close()
    basicFile.close()
    outFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
