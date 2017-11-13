#!/usr/bin/env python

"""
    By: David M. French (frenchd@astro.wisc.edu)
    
    $Id: pickGalaxiesNames4.py, v 4.0  08/30/16
    
    This code decides which galaxies to query ned with and compiles of list of them.
    OUTPUT: galaxyNamesUpdate.csv
    
    v2: try to speed it up with multiprocessing
    
    v4: back to basics, trying again
    
    
"""

import csv
import sys
import os
# import multiprocessing
import time
import getpass
from utilities import *

#########################################################################################

def main():
    # Decide which directory to use
    user = getpass.getuser()
    
    if user == 'frenchd':
        oldFileName = '/usr/users/frenchd/gt/NewGalaxyTable5.csv'
        newFileName = '/usr/users/frenchd/galaxy_table_code/GT_update2/compiledGalaxyNames_typeCut.csv'
        outputFileName = '/usr/users/frenchd/galaxy_table_code/GT_update2/newGalaxyNames.csv'
    
    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    fieldnames = ('name','ra','dec')

    # open the files
    oldGalaxyTable = open(oldFileName,'rU')
    allNewGalaxyNames = open(newFileName,'rU')
    writerOutFile = open(outputFileName,'wt')

    # read in tables
    oldGalaxyReader = csv.DictReader(oldGalaxyTable)
    allNewGalaxyReader = csv.DictReader(allNewGalaxyNames)
    
    # write fieldname headers to the output file
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    preferredNameList = []
    oldNameList = []
    altNamesList = []
    coordsList = []
    
    print 'Populating the old galaxy table lists...'
    
    for l in oldGalaxyReader:
        name = l['preferredName']
        oldName = l['oldName']
        altNames = l['alternativeNames']
        coords = l['degreesJ2000RA_Dec']
        
        preferredNameList.append(name)
        oldNameList.append(oldName)
        coordsList.append(coords)
        
        for a in altNames:
            altNamesList.append(a)
        
        
    print 'Finished, now commencing search...'
    print
    
    count = 0
    total = 122000
    for i in allNewGalaxyReader:
        count +=1
        percentComplete = round(float(count/total * 100),2)
        sys.stdout.write("Percent Complete %s %%. \r" %percentComplete)
        sys.stdout.flush()
        
        name = i['name']
        ra = i['ra']
        dec = i['dec']
        
        ra = ra.replace('h',':').replace('m',':').replace('s','').split(':')
        dec = dec.replace('d',':').replace('m',':').replace('s','').split(':')
        
        raDeg, decDeg = convertRAandDec(ra,dec,'degree')
        
        found = False
        if name in preferredNameList:
            found = True
        
        elif name in oldNameList:
            found = True
            
        elif name in altNamesList:
            found = True
        
        else:
            for c in coordsList:
                ra,dec = eval(c)
                
                if abs(float(ra) - float(raDeg)) <0.0001:
                    if abs(float(dec) - float(decDeg)) <0.0001:
                        found = True
                        print 'ra == raDeg -> {0} == {1}'.format(ra,raDeg)
                        print 'dec == decDeg -> {0} == {1}'.format(dec,decDeg)
                        print name
                        print            
            
        if not found:
            outputList = [name,raDeg,decDeg]
            row = dict((f,o) for f,o in zip(fieldnames,outputList))
            writer.writerow(row)
            
    oldGalaxyTable.close()
    allNewGalaxyNames.close()
    writerOutFile.close()
         
        
if __name__ == '__main__':
    main()
    
    