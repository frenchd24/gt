#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  nedphotometryretriever2.py, v 2.0 2/15/17

A program to grab photometric galaxy data from the NED server by means of xml-based VOTable
distributions and html parsing, and recompile specified data into a new csv file(s)

returns two files. One contains all useful photometry data, the second contains specific mined data
consisting of one chosen measurement in each of B, V, H, J, K, g, r, i bands, if it exists

v2: much simplified version (2/15/17)

'''

import optparse
import sys
import os
import csv
import string
import warnings
import urllib
import numpy
# import atpy

from utilities import *

from astropy.io.votable import parse, tree
# from vo.table import parse
# import vo.tree

import urllib2
import time

__author__ = "David M. French - frenchd@astro.wisc.edu"
__version__="2.0"


def parse_commandline():
    # Parse the options given on the command-line.

    parser = optparse.OptionParser(usage=__doc__,\
                                   version="$Id: nedphotometryretriever2,v 2.0 02/15/2017")
    parser.add_option("-f", "--filename", dest='filename', action='store', type='string',\
            help="The name of the file containing the galaxy names to be queried")
    parser.add_option("-o", "--output", dest='output',action='store',\
        default="/usr/users/frenchd/Desktop/NEDResults/",\
        help="Specify the directory in which to save newly generated csv table")
    parser.add_option("-n", "--outname", dest='outname',action='store',\
        help="Specify the output filename (not including extension) for the newly generated csv table")
    parser.add_option("-v","--verbose",action="store_true",\
        default=False, help="Run verbosely")
    parser.add_option("-t","--test",action="store_true",\
        default=False, help="Run in test mode;do not write results to disk")

    opts,args = parser.parse_args()
    
    # check if necessary input exists

    if opts.filename is None and opts.outname is None:
        print
        print "nedphotometryretriever2.py  version: 2.0   02/15/2017"
        print 
        print "A program to grab galaxy data from the NED server by means of xml-basedVOTable"
        print "distributions and html parsing, and recompile specified data into a new csv file(s)"
        print
        print "Normal usage example: python nedphotometryretriever2.py -f names.txt -n outFile.csv -o /Users/me/ -v"
        print
        print " -f is a text file of names of objects, one per line, for which data will be retrieved"
        print " -n is the name of the csv file that will be created and the NED data written to (don't include filename extension)"
        print " -o is the full pathname describing where the file should be saved"
        print " -v tells it to run verbosely"
        print " -t runs in test mode, does not write anything to disk "
        print
        sys.exit(1)
    
    test = True
    if opts.filename is None:
        print "--filename is a required parameter"
        test = False
    if opts.outname is None:
        print "--outname is a required parameter"
        test = False
    
    if test == False:
        sys.exit(1)

    # make an output directory if it doesn't exist yet
    if not os.path.exists(opts.output): os.makedirs(opts.output)
    
    # show parameters
    if opts.verbose:
        print
        print "********************** PARAMETERS ******************************"
        print 'filename: ', opts.filename
        print 'output directory: '; print opts.output;
        print 'output filename: ',opts.outname
        print
        print "************************ RESULTS *******************************"
    
    return opts

    
class galaxyClass(object):
    # a galaxy object containing it's information and methods of returning it
    
    def __init__(self,name,votable):
        self.name = name
        self.votable = votable
        
    def returnPhotometryArray(self):
        # find and return passband of photometry for each object
        
        totalArray = []
        warnings.simplefilter("ignore")
        try:
            table = self.votable.get_table_by_id("NED_PhotometricData")
            passbands = table.array["photo_col2"]
            photometryMeasurement = table.array["photo_col3"]
            uncertainty = table.array["photo_col4"]
            units = table.array["photo_col5"]
            frequency = table.array["photo_col6"]
            nedPM = table.array["photo_col7"]
            nedUncertainty = table.array["photo_col8"]
            nedUnits = table.array["photo_col9"]
            refcode = table.array["photo_col10"]
            significance = table.array["photo_col11"]
            publishedFrequency = table.array["photo_col12"]
            frequencyMode = table.array["photo_col13"]
            spatialMode = table.array["photo_col15"]
            qualifiers = table.array["photo_col16"]
            comments = table.array["photo_col17"]
            
            totalArray = zip(passbands,photometryMeasurement,uncertainty,units,frequency,nedPM,nedUncertainty,nedUnits,refcode,\
            significance,publishedFrequency,frequencyMode,spatialMode,qualifiers,comments)
                    
        except Exception, e:
            sys.stderr.write("\n Unable to return passband. Here is the error message built into the exception: \n %s \n" %e)
        
        if not totalArray:
            totalArray.append('x')
            
        return totalArray
        

def parseGalaxyNames(nameList):
    # format galaxy names for the url
    
    newNameList = []
    for name in nameList:
        nname = name.strip()
        nname = urllib.quote_plus(nname)
        nname = nname.replace('\n','')
        newNameList.append(nname)
    return newNameList
    

def createCSVTable(outFile,fieldnames):
    # creates and returns a DictReader object populated with header names
    
    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    return writer

###########################################################################

    
def main(opts):
    # assuming 'theFile' contains one name per line, read the file
    
#   theFile = open(opts.filename,'rU')
#   fileLines = theFile.readlines()
#   theFile.close()

    # assuming the input file is a standard galaxy table:
    theFile = open(opts.filename,'rU')
    reader = csv.DictReader(theFile)
    
    # quit after retrieving a certain number of results?
    maxRetrieve = 10000
    
    # a list to put just the names in
    fileLines = []
    for line in reader:
        preferredName = line['preferredName'].replace('\n\r','\n')
        fileLines.append(preferredName)

##########################################################################################
    
    # format the names in the list and return a new list
    newNameList = parseGalaxyNames(fileLines)
    
    # create and write output table header and fieldnames
    fieldnames = ('galaxyName','B','u','g','r','i','z','J','H','K','all')
    
    # depreciated; include the 'all' at the end of just one file
#     fieldnamesComplete = ('galaxyName','photometry')
#     writerOutFileComplete = open(opts.output+opts.outname+'Complete.csv','wt')
#     writerComplete = createCSVTable(writerOutFileComplete,fieldnamesComplete)

    # open a csv file to write out to - actually don't yet here
#     writerOutFile = open(opts.output+opts.outname+'.csv','wt')
#     writer = createCSVTable(writerOutFile,fieldnames)

    # check to see if file already exists
    fullPath = opts.output+opts.outname+'.csv'
    if os.path.exists(fullPath):
        # It exists. Ask what to do.
        answer = raw_input("%s already exists. Append results to this file? [y,n]" %opts.outname)
        while answer != 'y' and answer != 'n':
            answer = raw_input("Please answer with a 'y' or 'n': ")
        if answer == 'y':
            # open that file and read it as a dictionary csv, allowing it to be written to
            readerFile = open(fullPath,'rU')
            print 'opening: ',fullPath
            reader = csv.DictReader(readerFile)
            writerOutFile = open(opts.output+opts.outname+'1'+'.csv','wt')
            opts.outname = opts.outname+'1'
            writer = createCSVTable(writerOutFile,fieldnames)
            
            # increase the field size limit
            csv.field_size_limit(sys.maxsize)
            
            # also determine if the last name written to this file shows up in the 
            # names file. If so, assume that we should start from there and ignore
            # all names coming before. Else, start from the top
            oldName = 'x'
            for line in reader:
                oldName = line['galaxyName']
                
                lineList = [line['galaxyName'],\
                line['B'],\
                line['u'],\
                line['g'],\
                line['r'],\
                line['i'],\
                line['z'],\
                line['J'],\
                line['H'],\
                line['K'],\
                line['all']]
                
                lineRow = dict((f,o) for f,o in zip(fieldnames,lineList))
                writer.writerow(lineRow)
            
            readerFile.close()
            print 'Last entry in %s: '%opts.outname, oldName
            try:
                nameListIndex = 0
                index = 0
                for n in newNameList:
                    n = urllib.unquote_plus(n).replace('\n','').replace(' ','').strip()
                    n = n.replace('*','')
                    if n == oldName:
#                         print 'nname: ',nname
                        nameListIndex = index
                        break
                    else:
                        index +=1
#                 nameListIndex = newNameList.index(oldName)
            except Exception,e:
                # name not found so an exception is raised
                nameListIndex = False
                print 'e: ',e
            
            print 'nameListIndex: ',nameListIndex
            
            if nameListIndex:
                print "len - index: ",len(newNameList),', ',nameListIndex
                if len(newNameList) - nameListIndex > 2:
                    print 'Starting from {0} in object name list'.format(newNameList[nameListIndex+1:nameListIndex+2][0])
                    newNameList = newNameList  [nameListIndex+1:]
                else:
                    print 'Object {0} is the last object in the name list. Please rerun with a'.format(oldName)
                    print 'new list of objects to search for. Exiting...'
                    sys.exit()
            
            else:
                print 'Could not find {0} in object name list. Starting from the beginning.'.format(oldName)
                
        if answer == 'n':
            # ask for a new filename, or give the option of quitting
            aTwo = raw_input("Please enter new filename, or 'q' to quit: ")
            if aTwo == 'q':
                # quit the program
                sys.exit()
            else:
                if not opts.test:
                    # create a new file based on the newly entered file name
                    opts.outname = aTwo
                    writerOutFile = open(opts.output+opts.outname+'.csv','wt')
                    writer = createCSVTable(writerOutFile,fieldnames)

    else:
        if not opts.test:
            # It does not already exist. Create it as a dictionary csv file
            writerOutFile = open(opts.output+opts.outname+'.csv','wt')
            writer = createCSVTable(writerOutFile,fieldnames)
        
        
##########################################################################################
    
    theFile.close()
    
    # begin main data retrieving loop
    total = len(newNameList)
    counter = maxRetrieve
    totalStart = time.time()
    for name in newNameList:
        # write out stats on progress
        percentComplete = round(float((total-counter))/total *100,1)
        sys.stdout.write("\n")
        sys.stdout.write("Percent complete: %s %% \n" %percentComplete)
        sys.stdout.write("Galaxies left: %s \n" %counter)
        sys.stdout.write("\n")
        sys.stdout.write("Starting: %s \n" %name)
        
        startTime = time.time()
        
        hostURL = 'http://ned.ipac.caltech.edu/cgi-bin/nph-datasearch?objname=%s&meas_type=bot&ebars_spec=ebars&label_spec=no&x_spec=freq&y_spec=Fnu_jy&xr=-1&of=xml_main&search_type=Photometry' %name
        if opts.verbose: print 'hostURL: ',hostURL
        
        try:
            url = urllib2.urlopen(hostURL)
            if opts.verbose: print 'got url'
        except Exception, e:
            sys.stderr.write("\n Unable to return url file object. Here is the error message built into the exception: \n %s \n" %e)
            
        votable = 'x'
        warnings.simplefilter('ignore')
        try:
            if opts.verbose: print 'parsing votable'
            votable = parse(url,pedantic=False)
        except Exception,e:
            sys.stderr.write("\n Unable to parse voTable. Here is the message built into the exception: \n %s \n" %e)
        warnings.resetwarnings()
#       url.close()
        
        galaxy = galaxyClass(name,votable)
        sys.stdout.write("1...")
        sys.stdout.flush()
        
        totalArray = galaxy.returnPhotometryArray()
        sys.stdout.write("2...")
        sys.stdout.flush()
        
        
        # sort photometry information
        acceptables = ['b','u','g','r','i','z','j','h','k']
        
        photDict = {'b':[],'u':[],'g':[],'r':[],'i':[],'z':[],'j':[],'h':[],'k':[],'all':[]}
        completePhotometryArray = []
        if votable != 'x':
            for measurement in totalArray:
                # this is just the name of the measurement, e.g., 'B (m_B)' or 'g (SDSS PSF) AB'
                passband = str(measurement[0]).lower().strip()
                print 'passband: ',passband
                
                # now check if the first letter matches any of the specific bands I'm
                # looking for -> the ones in 'acceptables'
                
                if photDict.has_key(passband[0]):
                    # if yes, then add it to the corresponding dictionary array
                    
                    go = True
                    # if the first letter is h, j, k but is not '2MASS', don't put
                    # it in the corresponding bin
                    if passband[0] == 'h' and not bfind(passband,'2mass'):
                        go = False
                    if passband[0] == 'j' and not bfind(passband,'2mass'):
                        go = False
                    if passband[0] == 'k' and not bfind(passband,'2mass'):
                        go = False
                    
                    # if the first letter is u,g,r,i,z but is not 'sdss', don't
                    # put it in the corresponding bin
                    if passband[0] == 'u' and not bfind(passband,'sdss'):
                        go = False
                    if passband[0] == 'g' and not bfind(passband,'sdss'):
                        go = False
                    if passband[0] == 'r' and not bfind(passband,'sdss'):
                        go = False
                    if passband[0] == 'i' and not bfind(passband,'sdss'):
                        go = False
                    if passband[0] == 'z' and not bfind(passband,'sdss'):
                        go = False
                    
                    if go:
                        photDict[passband[0]].append(measurement)
                        
                # then add all of them into the 'all' array
                photDict['all'].append(measurement)

                
#         # round two of photometry elimination
#         masterPhotometryArray = []
        
        
#         # round one of photometry elimination
#         acceptables = ['b','v','h','j']
#         completePhotometryArray = []
#         if votable != 'x':
#             for measurement in totalArray:
#                 passband = str(measurement[0])[0].lower()
#                 for band in acceptables:
#                     if passband == band:
#                         if len(str(measurement[0])) > 1:
#                             if measurement[0][1] == ' ' or measurement[0][1] == '_':
#                                 completePhotometryArray.append(measurement)
#                         elif len(str(measurement[0])) ==1:
#                             completePhotometryArray.append(measurement)
#                         
#                 if passband == 'g' and 'sdss' in measurement[0].lower():
#                     completePhotometryArray.append(measurement)
#                     
#                 if passband == 'r' and 'sdss' in measurement[0].lower():
#                     completePhotometryArray.append(measurement)
#                     
#                 if passband == 'i' and 'sdss' in measurement[0].lower():
#                     completePhotometryArray.append(measurement)
#                 
#                 if passband == 'k':
#                     if 'sdss' in measurement[0].lower() or '2mass' in measurement[0].lower():
#                         completePhotometryArray.append(measurement)
#                 
#         # round two of photometry elimination
#         masterPhotometryArray = []
        
#         if not completePhotometryArray:
#             completePhotometryArray.append('x')
#             masterPhotometryArray.append('x')



        # remove URL formatting from the name
        oldName = urllib.unquote_plus(name).replace('\n','').strip()
        
        # make the list of entries for this row/galaxy
        rowValues = [oldName,\
        photDict['b'],\
        photDict['u'],\
        photDict['g'],\
        photDict['r'],\
        photDict['i'],\
        photDict['z'],\
        photDict['j'],\
        photDict['h'],\
        photDict['k'],\
        photDict['all']]
        
        # zip the values with the fieldnames
        completeRow = dict((f,v) for f,v in zip(fieldnames, rowValues))
        
        # write it to the file
        writer.writerow(completeRow)
            
        # advance the counter and print out the time status
        sys.stdout.write("\n")
        counter -=1
        
        queryTime = time.time() - startTime
        print 'query time: ',queryTime
        if queryTime <= 1:
            time.sleep(1 - queryTime)
    
        # quit if the max retrieval number is reached
        if counter == 0:
            print 'Reached maxRetrieve = {0}. Exiting...'.format(maxRetrieve)
            break
    
    print 'Done.'
    print
    print 'Total time to run: {0}'.format(time.time() - totalStart)
    print 'Results can be found in {0}{1}'.format(opts.output,opts.outname)
    print
    writerOutFile.close()
        
###############################################################################

if __name__=="__main__":
    # parse commandline
    commandlineOptions = parse_commandline()
    # do the work
    main(commandlineOptions)