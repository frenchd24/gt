#!/usr/bin/env python

'''
A program to grab galaxy data from the NED server by means of xml-based VOTable
distributions and html parsing, and recompile specified data into a new csv file(s)

Based on: NEDcompilerPA4.py, only modified  to read names out of csv type files instead
of txt files.

v5: updated for the galaxy table update (11/30/16)

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
# from vo.table import parse
# import vo.tree
from astropy.io.votable import parse,tree
import math

import Queue
import threading
import urllib2
import time

from utilities import *

# import glob
# import cPickle
# import shelve
# import tempfile


__author__ = "David M. French - frenchd@astro.wisc.edu"
__version__="5.0"


def parse_commandline():
    # Parse the options given on the command-line.

    parser = optparse.OptionParser(usage=__doc__,\
        version="$Id: NEDcompilerPA5,v 5.0 11/30/2016")
    parser.add_option("-f", "--filename", dest='filename', action='store', type='string',\
        help="The name of the file containing the galaxy names to be queried")
    parser.add_option("-o", "--output", dest='output',action='store',\
        default="/usr/users/frenchd/galaxy_table_code/GT_update2/",\
        help="Specify the directory in which to save newly generated csv table")
    parser.add_option("-n", "--outname", dest='outname',action='store',\
        help="Specify the output filename (etc.csv) for the newly generated csv table")
    parser.add_option("-v","--verbose",action="store_true",\
        default=False, help="Run verbosely")
    parser.add_option("-t","--test",action="store_true",\
        default=False,help="Run in test mode;do not write results to disk")
        
    opts,args = parser.parse_args()
    
    # check if necessary input exists
    
    if opts.filename is None and opts.outname is None:
        print
        print "NEDcompilerPA5.py  version: 5.0   11/30/2016"
        print 
        print "A program to grab galaxy data from the NED server by means of xml-basedVOTable"
        print "distributions and html parsing, and recompile specified data into a new csv file(s)"
        print
        print "Normal usage example: python NEDcompilerPA5.py -f names.txt -n outFile -o /Users/me/ -v"
        print
        print " -f is a text file of names of objects, one per line, for which data will be retrieved"
        print " -n is the name of the csv file (w/o extension) Nthat will be created and the NED data written to"
        print " -o is the full pathname describing where the file should be saved"
        print " -v tells it to run verbosely"
        print " -t runs the code in test mode, which only prints out the results without writing them to file"
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
    if opts.verbose and not opts.test:
        print
        print "********************** PARAMETERS ******************************"
        print 'filename: ', opts.filename
        print 'output directory: '; print opts.output;
        print 'output filename: ',opts.outname
        print
        print "************************ RESULTS *******************************"
    
    return opts


# def isNumber(s):
#     try:
#         float(s)
#         if math.isnan(s):
#             return False
#         else:
#             return True
#     except Exception:
#         return False
        

class ThreadUrl(threading.Thread):
    """Threaded Url Grab"""
    def __init__(self, queue,outQueue,numberQueue):
        threading.Thread.__init__(self)
        self.queue = queue
        self.outQueue = outQueue
        self.numberQueue = numberQueue

    def run(self):
        while True:
            #grabs host from queue
            host = self.queue.get()
            
            #grabs urls of hosts and reads or parses them according to their type
            try:
#               s1 = time.time()
                url = urllib2.urlopen(host[0])
#               print 'fetch time: ',time.time() - s1
            except Exception, e:
                sys.stderr.write("\n Unable to return url file object. Here is the error message built into the exception: \n %s\n" %e)
                    
            if host[1] ==1:
                warnings.simplefilter("ignore")
                votable1 = 'x'
                try:
                    votable1 = parse(url,pedantic=False)
                except Exception, e:
                    sys.stderr.write("\n Unable to parse voTable. Here is the message built into the exception: \n %s \n" %e)
                
                self.outQueue.put(votable1)
                warnings.resetwarnings()
                
            else:
                html = 'x'
                try:
                    html = url.readlines()
                except Exception, e:
                    sys.stderr.write("\n Unable to read html. Here is the message built into the exception: \n %s \n" %e)
                    
                self.outQueue.put(html)
            self.numberQueue.put(host[1])
            self.queue.task_done()
            url.close()
    
    
class galaxyClass(object):
    # a galaxy object containing it's information and methods of retrieving and 
    # returning it
    
    def __init__(self,name,votable,opts):
        self.name = name
        self.votable = votable
        
        
    def returnDiametersandPA(self,opts):
        # find and return major and minor diameters
        
        chosen = ['x',('x','x'),'x','x']
        completeList = ['x']
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_Diameters_Data")
            freqTargeted = table.array['d_col2'].tolist()
            majorD = table.array['d_col19'].tolist()
            minorD = table.array['d_col26'].tolist()
            dRatio = table.array['d_col21'].tolist()
            pa = table.array['d_col25'].tolist()
            
            # make complete list
            completeList = []
            for i in range(len(freqTargeted)):
                completeList.append([freqTargeted[i],(majorD[i],minorD[i]),dRatio[i],pa[i]])
            
            # now decide which set of values to use
            
            specialSets = []
            completeSets = []
            
            majorDList = []
            for m in range(len(majorD)):
                if isNumber(majorD[m]):
                    majorDList.append(majorD[m])
                    
            maxD = max(majorDList)
            minD = min(majorDList)
            aveD = (maxD+minD)/2
            if opts.verbose:
                print 'max diameter: ',maxD
                print 'average diameter: ',aveD
            
            for m in range(len(freqTargeted)):
                if freqTargeted[m].find('SDSS')!=-1 or freqTargeted[m].find('2MASS')!=-1:
                    specialSets.append(m)
            
            # if SDSS or 2MASS values are present, pick only those out of all possible
            if len(specialSets) !=0:
                usedSet = specialSets
            else:
                usedSet = range(len(freqTargeted))
            
            print 'used set: ',usedSet
            for m in usedSet:
                score = 0
                if isNumber(majorD[m]):
                    if majorD[m] == maxD:
                        if isNumber(minorD[m]):
                            if minorD[m]<majorD[m]:
                                score +=3
                                print 'add 3 points to %s: ' % freqTargeted[m]
                            else:
                                score -=5
                                print 'removing 5 points (minor = major) from %s: ' % freqTargeted[m]
                        else:
                            score +=2
                            print 'adding 2 points (no minor available) to %s: ' % freqTargeted[m]
                    elif majorD[m] >= aveD:
                        if isNumber(minorD[m]):
                            if minorD[m]<majorD[m]:
                                score +=2
                                print 'adding 2 points (< maxD and minor available) to %s: ' % freqTargeted[m]
                            else:
                                score -=5
                                print 'removing 5 points (< maxD, minor = major) from %s: ' % freqTargeted[m]
                        else:
                            score+=1
                            print 'adding 1 point (no minor available) to %s: ' % freqTargeted[m]
                    else:
                        print 'removing 5 points (< aveD) from %s: ' % freqTargeted[m]
                        score-=5
                    if majorD[m] == minorD[m]:
                        score-=5
                        print 'removing 5 points (major = minor final) from %s: ' %freqTargeted[m]
                if isNumber(minorD[m]):
                    score+=.1
                if isNumber(dRatio[m]):
                    score+=.4
                if isNumber(pa[m]):
                    score+=15
                if freqTargeted[m].find('SDSS'):
                    score+=5
                if freqTargeted[m].find('2MASS'):
                    score+=4
                
                completeSets.append((score,m))
                
            completeSets.sort(reverse=True)
                    
            print 'table: ',table
            print 'freqTargeted: ',freqTargeted
            print 'majorD: ',majorD
            print 'minorD: ',minorD
            print 'dRatio: ',dRatio
            print 'pa: ',pa
            print 'completeSets: ',completeSets
            score,m = completeSets[0]
            print 'chosen: ',freqTargeted[m],' = ',majorD[m],', ',minorD[m],', ',dRatio[m],' : ',pa[m]
            print
            chosen = [freqTargeted[m],(majorD[m],minorD[m]),dRatio[m],pa[m]]
            warnings.resetwarnings()
            
        except Exception, e:
            sys.stderr.write("\n Unable to return major or minor diameter. Here is the error message "\
            "built into the exception:\n %s\n" %e)
            chosen = ['x',('x','x'),'x','x']
            
        return chosen,completeList

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


def calcBestDist(rid, distvcorr):
    # calculate a weighted best distance value. This value is used for all further
    # calculations. Weight is 3/4 redshift-independent (rid) to 1/4 corrected redshift velocity
    
    if type(rid) == "<type 'tuple'>":
        rid = rid[0]
    
    if isNumber(rid):
        if isNumber(distvcorr):
            if float(distvcorr)<= 1000:
        
                # if the galaxy is too close, rely only on the real measure
                bestDist = float(rid)
            
            else:
                # if both are available, average them with weights
                bestDist = (float(rid)*3 + float(distvcorr))/4
        else:
            # if no velocity is available, use real measure
            bestDist = rid
            
    # no redshift independent distance measure available        
    elif isNumber(distvcorr):
        bestDist = distvcorr
        
    else:
        # if no distance measures are available, mark it 'x'
        bestDist = 'x'
        
    return bestDist

###########################################################################

    
def main(opts):
    # assuming 'theFile' contains one name per line, read the file
    
    start = time.time()
    maxRetrieve = 10000
    
    csvFile = False
    
    # a list to put just the names in
    if csvFile:
        # assuming the input file is a standard galaxy table:
        theFile = open(opts.filename,'rU')
        reader = csv.DictReader(theFile)
        
        fileLines = []
        for line in reader:
            preferredName = line['preferredName'].replace('\n\r','\n')
            fileLines.append(preferredName)
        
        theFile.close()
        
    else:
        theFile = open(opts.filename,'r')
        fileLines = theFile.readlines()
        theFile.close()
    
    # format the names in the list and return a new list
    newNameList = parseGalaxyNames(fileLines)
    
    fieldnames = ('oldName','freqTargeted','diameters','dRatio','pa','completeList')
    
    
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
            
            # also determine if the last name written to this file shows up in the 
            # names file. If so, assume that we should start from there and ignore
            # all names coming before. Else, start from the top
            oldName = 'x'
            for line in reader:
                oldName = line['oldName']
                
                lineList = [line['oldName'],\
                line['freqTargeted'],\
                line['diameters'],\
                line['dRatio'],\
                line['pa'],\
                line['completeList']]
                
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
    
#     
#     if not opts.test:
#         writerOutFile = open(opts.output+opts.outname+'.csv','wt')
#         writer = createCSVTable(writerOutFile,fieldnames)
    
    queue = Queue.Queue()
    outQueue = Queue.Queue()
    numberQueue = Queue.Queue()
    
    t1 = ThreadUrl(queue,outQueue,numberQueue)
#   t2 = ThreadUrl(queue,outQueue,numberQueue)
#   t3 = ThreadUrl(queue,outQueue,numberQueue)
    t1.setDaemon(True)
#   t2.setDaemon(True)
#   t3.setDaemon(True)
    t1.start()
#   t2.start()
#   t3.start()
        
    total = len(newNameList)
    counter = maxRetrieve
#     counter = total
    for name in newNameList:
        # return information for each object 
        
        totalStart = time.time()
        percentComplete = round(float((total-counter))/total *100,1)
        sys.stdout.write("\n")
        sys.stdout.write("Percent complete: %s %% \n" %percentComplete)
        sys.stdout.write("Galaxies left: %s \n" %counter)
        sys.stdout.write("\n")
        sys.stdout.write("Starting: %s \n" %name)       

        hosts = [["http://ned.ipac.caltech.edu/cgi-bin/datasearch?objname=%s&search_type=Diameters&of=xml_main"\
        %name,1]]
        
        #populate queue with data
        for host in hosts:
            time.sleep(0.9)
            queue.put(host)
    
        #wait on the queue until everything has been processed
        queue.join()
        
        # decide which object is which in the queue
        votable = outQueue.get()
        
        galaxy = galaxyClass(name,votable,opts)
#       galaxy.queryNED()
        sys.stdout.write("1...\n")
        sys.stdout.flush()
        
        diametersandPA,completeList = galaxy.returnDiametersandPA(opts)
        freqTargeted = diametersandPA[0]
        diameters = diametersandPA[1]
        dRatio = diametersandPA[2]
        pa = diametersandPA[3]
        
        if not isNumber(diameters[0]) and str(diameters).find('x')==-1:
            diameters = ('x','x')
        if not isNumber(dRatio) and str(dRatio).find('x')==-1:
            dRatio = 'x'
        if not isNumber(pa) and str(pa).find('x')==-1:
            pa = 'x'
        
        sys.stdout.write("8...")
        sys.stdout.flush()  
        
        oldName = urllib.unquote_plus(name).replace('\n','').strip()
        
        objectInfoList = [oldName.replace(' ',''),freqTargeted,diameters,dRatio,pa,completeList]
        
        # write this data to a row of the csv file
        if not opts.test:
            row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
            if opts.verbose:
                print '\n','row: ',row
            writer.writerow(row)
        
        # update the status bar
        sys.stdout.write("\n")
        print 'Total retrieval time: ',time.time() - totalStart
        counter -=1
        
        if counter == 0:
            print 'Reached maxRetrieve = {0}. Exiting...'.format(maxRetrieve)
            break
        
        
        
    if not opts.test:
        writerOutFile.close()
        
    if opts.verbose:
        print "Done."
        print 
        if not opts.test:
            print "Results can be found in: %s%s" %(opts.output,opts.outname)
            print
        print "Elapsed Time: %s" % (time.time() - start)
        print
    
###############################################################################

if __name__=="__main__":
    # parse commandline
    commandlineOptions = parse_commandline()
    # do the work
    main(commandlineOptions)