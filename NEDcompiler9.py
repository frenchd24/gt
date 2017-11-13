#!/usr/bin/env python

'''A program to grab galaxy data from the NED server by means of xml-based VOTable
distributions and html parsing, and recompile specified data into a new csv file(s)

version 7 and 8: Updates and adding the functionality to continue a search where
    the previous run stopped within a list of names
   
   
v9: major updates to combine with photometry retriever, update the galaxy table for 
    final publication (08/29/16)
    
     
'''

import getpass
import optparse
import sys
import os
import csv
import string
import warnings
import urllib
import numpy
# import atpy
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
__version__="9.0"


def parse_commandline():
    # Parse the options given on the command-line.

    parser = optparse.OptionParser(usage=__doc__,\
                                   version="$Id: NEDcompiler9,v 9.0 08/29/2016")
    parser.add_option("-f", "--filename", dest='filename', action='store', type='string',\
            help="The name of the file containing the galaxy names to be queried")
    parser.add_option("-o", "--output", dest='output',action='store',\
        default="/usr/users/frenchd/Galaxy Table code/NED Update/",\
        help="Specify the directory in which to save newly generated csv table")
    parser.add_option("-n", "--outname", dest='outname',action='store',\
        help="Specify the output filename (etc.csv) for the newly generated csv table")
    parser.add_option("-v","--verbose",action="store_true",\
        default=False, help="Run verbosely")    

    opts,args = parser.parse_args()
    
    # check if necessary input exists
    
    if opts.filename is None and opts.outname is None:
        print
        print "NEDcompiler9.py  version: 9.0   08/29/2016"
        print 
        print "A program to grab galaxy data from the NED server by means of xml-basedVOTable"
        print "distributions and html parsing, and recompile specified data into a new csv file(s)"
        print
        print "Normal usage example: python NEDcompiler9.py -f names.txt -n outFile -o /Users/me/ -v"
        print
        print " -f is a text file of names of objects, one per line, for which data will be retrieved"
        print " -n is the name of the csv file that will be created and the NED data written to"
        print " -o is the full pathname describing where the file should be saved"
        print " -v tells it to run verbosely"
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
    
# def isNumber(x):
#     try:
#         float(x)
#         return True
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
    
    def __init__(self,name,votable,html,fullHtml):
        self.name = name
        self.votable = votable
        self.html = html
        self.fullHtml = fullHtml
        
#   def queryNED(self):
#       # query NED server for each object name, returning the file object representing
#       # the xml VOTable output from NED as well as the redshift-independent distance html file
#       # Note: htmlFileObject is the html page for redshift independent distance, while
#       # fullHtmlFileObject is the whole html page for the object
# #         queue = Queue.Queue()
# #         outQueue = Queue.Queue()
# #         voList = []
# 
#       hosts = [["http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname=%s"\
#       %self.name,1],["http://ned.ipac.caltech.edu/cgi-bin/nDistance?name=%s"\
#       %self.name,2],["http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES"\
#       %self.name,3]]
#           
#       #spawn a pool of threads, and pass them queue instance
# #         for i in range(len(hosts)):
# #             t = ThreadUrl(queue,outQueue,voList)
# #             t.setDaemon(True)
# #             t.start()
# #             time.sleep(0.9)
# 
#       #populate queue with data
#       for host in hosts:
#           queue.put(host)
#   
#       #wait on the queue until everything has been processed
#       queue.join()
#       
#       # decide which object is which in the queue
#       first = outQueue.get()
#       second = outQueue.get()
#       third = outQueue.get()
#       
#       if voList[0] == 1:
#           self.votable = first
#       elif voList[0] == 2:
#           self.html = first
#       elif voList[0] == 3:
#           self.fullHtml = first
#           
#       if voList[1] == 1:
#           self.votable = second
#       elif voList[1] == 2:
#           self.html = second
#       elif voList[1] == 3:
#           self.fullHtml = second
#           
#       if voList[2] == 1:
#           self.votable = third
#       elif voList[2] == 2:
#           self.html = third
#       elif voList[2] == 3:
#           self.fullHtml = third
            
        
    def returnRedshift(self):
        # find and return redshift for each object
        
        redshift = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_MainTable")
            redshift = table.array['Redshift'][0]
            warnings.resetwarnings()
            
        except Exception, e:
            sys.stderr.write("\n Unable to return redshift. Here is the error message "\
            "built into the exception:\n %s\n" %e)
            
        if isNumber(redshift):
            if math.isnan(float(redshift)):
                redshift = 'x'
        return redshift

    def returnJ2000Position(self):
        # find and return equatorial J2000 position RA and Dec for each object
        
        ra = 'x'
        dec = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_MainTable")
            ra = table.array['RA(deg)'][0]
            dec = table.array['DEC(deg)'][0]
            warnings.resetwarnings()
            
        except Exception, e:
            sys.stderr.write("\n Unable to return J2000 position. Here is the error message "\
            "built into the exception:\n %s\n" %e)
            
        if isNumber(ra) or isNumber(dec):
            if math.isnan(float(ra)) or math.isnan(float(dec)):
                ra = 'x'
                dec = 'x'
        return ra,dec

    def returnGalactic(self):
        # find and return galactic position
        
        galacticLong = 'x'
        galacticLat = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_PositionDataTable")
            galacticLong = table.array['pos_lon_gal_d'][0]
            galacticLat = table.array['pos_lat_gal_d'][0]
            warnings.resetwarnings()
            
        except Exception, e:
            sys.stderr.write("\n Unable to return galactic position. Here is the error message "\
            "built into the exception:\n %s\n" %e)
            
        if isNumber(galacticLong) or isNumber(galacticLat):
            if math.isnan(float(galacticLong)) or math.isnan(float(galacticLat)):
                galacticLong ='x'
                galacticLat = 'x'
        return galacticLong,galacticLat

    def returnRedIndependentDist(self):
        # find and return redshift-independent distance with mean, std. dev., min and max
        
        mean = 'x'
        stddev = 'x'
        min = 'x'
        max = 'x'
        
        try:
            for line in self.html:
                index = line.find('Mean')
                if index != -1:
                    indexMagStart = line[index:].find('<td>') + index
                    indexMagEnd = line[indexMagStart:].find('</td>') + indexMagStart
                    indexMpcStart = line[indexMagEnd:].find('<td>') + indexMagEnd
                    indexMpcEnd = line[indexMpcStart:].find('</td>') + indexMpcStart            
                    mean = line[indexMpcStart+4:indexMpcEnd]
                    
                index = line.find('Std. Dev.')
                if index != -1:
                    indexMagStart = line[index:].find('<td>') + index
                    indexMagEnd = line[indexMagStart:].find('</td>') + indexMagStart
                    indexMpcStart = line[indexMagEnd:].find('<td>') + indexMagEnd
                    indexMpcEnd = line[indexMpcStart:].find('</td>') + indexMpcStart            
                    stddev = line[indexMpcStart+4:indexMpcEnd]
                    
                index = line.find('Min.')
                if index != -1:
                    indexMagStart = line[index:].find('<td>') + index
                    indexMagEnd = line[indexMagStart:].find('</td>') + indexMagStart
                    indexMpcStart = line[indexMagEnd:].find('<td>') + indexMagEnd
                    indexMpcEnd = line[indexMpcStart:].find('</td>') + indexMpcStart            
                    min = line[indexMpcStart+4:indexMpcEnd]
                    
                index = line.find('Max.')
                if index != -1:
                    indexMagStart = line[index:].find('<td>') + index
                    indexMagEnd = line[indexMagStart:].find('</td>') + indexMagStart
                    indexMpcStart = line[indexMagEnd:].find('<td>') + indexMagEnd
                    indexMpcEnd = line[indexMpcStart:].find('</td>') + indexMpcStart            
                    max = line[indexMpcStart+4:indexMpcEnd]
                                
        except Exception, e:
            sys.stderr.write("\n Unable to find redshift-independent distance. Here's the error "\
                "message built into the exception: \n %s\n" %e)
        
        if isNumber(mean):
            if math.isnan(float(mean)):
                mean = 'x'
        if isNumber(stddev):
            if math.isnan(float(stddev)):
                stddev = 'x'
        if isNumber(min):
            if math.isnan(float(min)):
                min = 'x'
        if isNumber(max):
            if math.isnan(float(max)):
                max = 'x'
        return mean,stddev,min,max
    
    def returnMorphology(self):
        # find and return NED homogenized galaxy morphology
        
        morphology = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_BasicDataTable")
            morphology = table.array['morph_type'][0]
            warnings.resetwarnings()
            
        except Exception, e:
            sys.stderr.write("\n Unable to return galaxy morphology. Here is the error message "\
            "built into the exception:\n %s\n" %e)
        
        if morphology == '' or morphology == ' ' or morphology == '\n' or morphology == ' \n':
            morphology = 'x'
        elif isNumber(morphology):
            if math.isnan(float(morphology)):
                morphology = 'x'
        return morphology

    def returnDistIndicator(self):
        # find and return distance indicator (under Classifications)
        
        distIndicator = 'x'
        try:
            for line in self.fullHtml:
                index = line.find('Distance Indicator')
                if index != -1:
                    index2 = line[index:].find('<TD>')
                    index3 = line[index2+index:].find('</TD>')
                    distIndicator = line[index2+index+4:index3+index+index2]
                    break
                    
        except Exception, e:
            sys.stderr.write("\n Unable to return distance indicator. Here is the error message "\
            "built into the exception:\n %s\n" %e)
            
        if isNumber(distIndicator):
            if math.isnan(float(distIndicator)):
                distIndicator = 'x'
        return distIndicator

    def returnLuminosityClass(self):
        # find and return luminosity class (under Classifications)
        
        luminosityClass = 'x'
        try:
            for line in self.fullHtml:
                index = line.find('Luminosity Class')
                if index != -1:
                    index2 = line[index:].find('<TD>')
                    index3 = line[index2+index:].find('</TD>')
                    index4 = line[index3+index2+index:].find('<TD>')
                    index5 = line[index4+index3+index2+index:].find('</TD>')
                    luminosityClass = line[index4+index3+index2+index+4:index5+index4+index3+index+index2]
                    break
                    
        except Exception,e:
            sys.stderr.write("\n Unable to return luminosity class. Here is the error message "\
            "built into the exception:\n %s\n" %e)
        
        if isNumber(luminosityClass):
            if math.isnan(float(luminosityClass)):
                luminosityClass = 'x'
        return luminosityClass
        
    def returnEBminusV(self):
        # find and return foreground galactic extinction E(B-V)
        
        EBminusV = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_BasicDataTable")
            EBminusV = table.array['gal_extinc_E(B-V)'][0]
            warnings.resetwarnings()
            
        except Exception, e:
            sys.stderr.write("\n Unable to return EBminusV. Here is the error message "\
            "built into the exception:\n %s\n" %e)
        
        if EBminusV ==' ' or EBminusV == '' or not isNumber(EBminusV):
            EBminusV = 'x'
        elif math.isnan(float(EBminusV)):
            EBminusV = 'x'
        return EBminusV
        
    def returnRadialVelocity(self):
        # find and return radial velocity
        
        radialVelocity = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_MainTable")
            radialVelocity = table.array['main_col6'][0]
            warnings.resetwarnings()
            
        except Exception, e:
            sys.stderr.write("\n Unable to return Radial Velocity. Here is the error message "\
            "built into the exception:\n %s\n" %e)
        
        if radialVelocity == ' ' or radialVelocity == '' or not isNumber(radialVelocity):
            radialVelocity = 'x'
        elif math.isnan(float(radialVelocity)):
            radialVelocity = 'x'
        return radialVelocity
        
    
    def returnDiameters(self):
        # find and return major and minor diameters
        
        majorDiameter = 'x'
        minorDiameter = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_BasicDataTable")
            majorDiameter = table.array['diam_maj'][0]
            minorDiameter = table.array['diam_min'][0]
            warnings.resetwarnings()
            
        except Exception, e:
            sys.stderr.write("\n Unable to return major or minor diameter. Here is the error message "\
            "built into the exception:\n %s\n" %e)
        
        if majorDiameter == ' ' or majorDiameter == '' or not isNumber(majorDiameter):
            majorDiameter = 'x'
        if minorDiameter == ' ' or minorDiameter == '' or not isNumber(minorDiameter):
            minorDiameter == 'x'
        if isNumber(majorDiameter):
            if math.isnan(float(majorDiameter)):
                majorDiameter = 'x'
        if isNumber(minorDiameter):
            if math.isnan(float(minorDiameter)):
                minorDiameter = 'x'
        return majorDiameter,minorDiameter
        
        
    def returnPhotometry(self):
        # find and return quick-look photometry information
        
        photometry = []
        try:
            k = 0
            lineNumber = 0
            found = False
            for line in self.fullHtml:
                k +=1
                index = line.find('Apparent Mag or Flux')
                if index != -1 and found == False:
#                   print 'line: ',line
                    lineNumber = k
                    found = True
                if found and k >= lineNumber +1:
#                   print 'found and k>=lineNumber+1: ',line
                    if line.find('<b>NOTE: </b>The above') != -1:
#                       print 'found end, breaking: ',line
                        break
                    else:
                        indexBand = line.find('<td>')+4
                        band = line[indexBand:indexBand+2]
                        for i in ['U ','B ','V ','R ','H ','NU','U_','B_','V_','R_','H_','b_']:
                            if band == i:
                                indexBandEnd = line[indexBand:].find('</td>')+indexBand
                                fullBand = line[indexBand:indexBandEnd].replace(' ','')
                                
                                indexApparentMag = line[indexBandEnd:].find('<td>')+4 + indexBandEnd
                                indexApparentMagEnd = line[indexApparentMag:].find('</td>')+indexApparentMag
                                fullApparentMag = line[indexApparentMag:indexApparentMagEnd].replace(' ','')
                                
#                               indexRef = line[indexApparentMagEnd:].find('<td>')+4 + indexApparentMagEnd
#                               indexRefEnd = line[indexRef:].find('</td>') + indexRef
                                indexRef = line.find('TARGET="ned_dw">') +16
                                indexRefEnd = line[indexRef:].find('</a>') + indexRef
#                               fullRef = line[indexRef:indexRefEnd].replace(' ','')
                                fullRef = line[indexRef:indexRefEnd]
                                
                                indexAbsoluteMag = line[indexRefEnd:].find('<td>')+4 + indexRefEnd
                                indexAbsoluteMagEnd = line[indexAbsoluteMag:].find('</td>')+indexAbsoluteMag
                                fullAbsoluteMag = line[indexAbsoluteMag:indexAbsoluteMagEnd].replace(' ','')
                                
                                indexBolometric = line[indexAbsoluteMagEnd:].find('<td>')+4 + indexAbsoluteMagEnd
                                indexBolometricEnd = line[indexBolometric:].find('</td>')+indexBolometric
                                fullBolometric = line[indexBolometric:indexBolometricEnd].replace(' ','')
                                
                                photometry.append((fullBand,fullApparentMag,fullAbsoluteMag,\
                                    fullBolometric,fullRef))
                    
        except Exception,e:
            sys.stderr.write("Unable to return photometry data. Here is the error message "\
            "built into the exception:\n %s\n" %e)
        
        if len(photometry) == 0:
            photometry.append('x')
        return photometry
        
        
    def returnNames(self):
        # find and return a list of all the names for a given object
        
        formattedNames = []
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_NamesTable")
            names = table.array['name_col1']
            parsedNames = parseGalaxyNames(names)
            for name in parsedNames:
                formattedName = urllib.unquote_plus(name).replace('\n','').strip()
                formattedNames.append(formattedName)
            warnings.resetwarnings()

        except Exception, e:
            sys.stderr.write("\n Unable to return alternative names. Here is the error message "\
            "built into the exception:\n %s\n" %e)
            return 'x'
        if len(formattedNames) ==0:
            formattedNames.append('x')
        return formattedNames
        

def parseGalaxyNames(nameList):
    # format galaxy names for the url
    
    newNameList = []
    for name in nameList:
        nname = name.strip()
        nname = nname.replace('*','')
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
    

def decidePreferredName(listOfAlternatives,oldName):
    # decides which name to use as the first name, returns that name
    
    print 
    bestName = None
    found = False
    
    # priority:
    first = ['NGC ','IC ','UGC ']
    second = ['MRK ','MCG ','ISO ','SBS ']
    third = ['CGCG ','IRAS ','RXJ ','FGC ','KUG ','PGC ','SDSS ','VCC ']
    fourth = ['2MASS ','2DF ','6DF ','HIPASS ','2MASX ']
    
    for f in first:
        if oldName.find(f) != -1:
            bestName = oldName
            found = True
            break
    
    if not found:
        for name in listOfAlternatives:
            for f in first:
                if name.find(f) != -1 and name.find('[') == -1:
                    bestName = name
                    found = True
                    break
                    
            if not found:
                for s in second:
                    if name.find(s) != -1 and name.find('[') == -1:
                        bestName = name
                        found = True            
                        break
                        
            if not found:
                for t in third:
                    if name.find(t) != -1 and name.find('[') == -1:
                        bestName = name
                        found = True
                        break
    
            if not found:
                for q in fourth:
                    if name.find(q) != -1 and name.find('[') == -1:
                        bestName = name
                        found = True
                        break
            if not found:
                bestName = oldName
    
    print 'bestName: ',bestName
#   if bestName != oldName:
#       listOfAlternatives.remove(bestName)
#       listOfAlternatives.append(oldName)
    return bestName
    
    
def pickPreferredName2(altnames, oldname):
    # oldname is the original name, and possibly not part of altnames
    # altnames is the list of alternative names (probably returned from NED)

    order = {1:'MESSIER',2:'NGC',3:'MRK',4:'UGC',5:'PHL',6:'3C',7:'IC',\
    8:'SBS',9:'MCG',10:'ISO',11:'TON',12:'PGC',13:'PG',14:'PB',15:'FGC',16:'HS',17:'HE',\
    18:'KUG',19:'IRAS',20:'RX',21:'CGCG',22:'FBQS',23:'LBQS',24:'SDSS',25:'VCC',26:'2MASS',\
    27:'2DF',28:'6DF',29:'HIPASS',30:'2MASX'}

    scores = []
    found = False
    for n in altnames:
        for i in range(1,31):
            if bfind(n,order[i]):
                scores.append([i,n])
                found = True
                break

        if not found:
            scores.append([99,n])

    scores.sort()
    bestName = scores[0][1]

    sortedNames = []
    for n in scores:
        sortedNames.append(n[0])

    if scores[0][0] == 99:
        bestName = oldname

    if bfind(bestName,'SDSS'):
        bestName = oldname

    return bestName, sortedNames
    
    
def pickPreferredName3(altNames, oldName):
    # oldname is the original name, and possibly not part of altnames
    # altnames is the list of alternative names (probably returned from NED)

    order = ['MESSIER','NGC','MRK','UGC','PHL','3C','IC','SBS','MCG','ISO','TON','PGC',\
    'PG','PB','FGC','HS','HE','KUG','IRAS','RX','CGCG','FBQS','LBQS','SDSS','VCC',\
    '2MASS','2DF','6DF','HIPASS','2MASX']
    
    # add oldName to the list of alternate names if it is not already there
    if len(altNames) >=1:
        if not bfind(str(altNames),str(oldName)):
            try:
                altNames.append(oldName)
            except Exception,e:
                sys.stdout.write("Issue with picking preferred name: {0}".format(e))
                sys.stdout.write("altNames list was: {0}".format(altNames))
    else:
        altNames = [oldName]
    
    found = False
    for i in order:
        for n in altNames:
            # n is input string, i is what i'm looking for
            if bfind(n,i) and not found:
                finalName = n
                found = True
                break

    if not found:
        finalName = oldName

    return finalName
    
    
    
    
# def returnRAandDec(ra,dec):
#     # this function converts ra and dec in degrees to HH:MM:SS and DD:MM:SS
#     
#     raHours,raMinutes = str(((float(ra)*24)/360)).split('.')
#     raMinutes,raSeconds = str(float('0.'+raMinutes) *60).split('.')
#     raSeconds = float('0.'+raSeconds) *60
# 
#     decDegrees,decMinutes = str(dec).split('.')
#     decMinutes,decSeconds = str(float('0.'+decMinutes)*60).split('.')
#     decSeconds = float('0.'+decSeconds)*60
#     
#     return (raHours,raMinutes,round(raSeconds,2)),(decDegrees,decMinutes,round(decSeconds,2))
# 
# 
# def calculatevcorr(ra,dec,velocity):
#     rav = 186.7833
#     decv = 12.9333
#     vcorr = velocity + 300*(math.sin(math.radians(dec)) * math.sin(math.radians(decv)) + \
#     math.cos(math.radians(dec))*math.cos(math.radians(decv)) * math.cos(math.radians(ra-rav)))
#     return vcorr
    
    
# def returnLinDiameters(major,minor,distance):
#     # input major and minor in arcsec, distance in Mpc
#     # outputs major and minor in kpc
#     newMajor = math.tan(math.radians(float(major)))*(distance*1000)
#     newMinor = math.tan(math.radians(float(minor)))*(distance*1000)
#     return (newMajor,newMinor)
    
    
# def returnInclination(major,minor):
#     # outputs inclination in degress
#     inclination = math.acos(float(minor)/float(major))*180/math.pi
#     return inclination
    

# def returnAngDiameters(major,minor,distance):
#     # input distances in mpc, major and minor is in kpc
#     # outputs angular diameters in arcsec
#     newMajor = math.atan((float(major)/1000)/float(distance))*(1/3600)
#     newMinor = math.atan((float(minor)/1000)/float(distance))*(1/3600)
#     return (newMajor,newMinor)
    
    
    
def returnGalaxyName(ra,dec,radius):
    # INPUT: ra and dec in sexagesimal for a galaxy, radius of search around that point
    #
    # OUTPUT: returns the NED preferred name for the first (closest) object in the 
    # search radius results
    
    # format RA and Dec
    newra = ra.replace('+','%2B')
    newra = newra.replace('-','%2D')
    newdec = dec.replace('+','%2B')
    newdec = newdec.replace('-','%2D')
    
    mainhost = 'http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon={0}&lat={1}&radius=0.5&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=xml_main&zv_breaker=30000.0&list_limit=5&img_stamp=YES'.format(newra,newdec)
    try:
        mainurl = urllib2.urlopen(mainhost)
        print 'opened mainurl',mainurl
    except Exception, e:
        sys.stderr.write("\n Unable to return url file object. Here is the error message built into the exception: \n %s\n" %e)
    
    warnings.simplefilter("ignore")
    mainvotable = 'x'
    try:
        mainvotable = parse(mainurl,pedantic=False)
        print 'parsed mainvotable',mainvotable
    except Exception, e:
        sys.stderr.write("\n Unable to parse voTable. Here is the message built into the exception: \n %s \n" %e)
    
    warnings.resetwarnings()
    mainurl.close()
    
    name = False
    try:
        warnings.simplefilter("ignore")
        maintable = mainvotable.get_table_by_id("NED_MainTable")
        mainName = maintable.array['main_col2'][0]
        mainName = str(mainName).strip()
        name = urllib.quote_plus(mainName)
        name = name.replace('\n','')
        print 'found name: ',name
        warnings.resetwarnings()
    except Exception, e:
        sys.stderr.write("\n Unable to return alternative names. Here is the error message "\
        "built into the exception:\n %s\n" %e)
        name = False
        
    return name
    
    
###########################################################################

    
def main(opts):
    # assuming 'theFile' contains one name per line, read the file
    
    hubbleConstant = 71
    
    # set a limit to how many queries to make?
    maxRetrieve = 10000
    
    start = time.time()

    theFile = open(opts.filename,'rU')
#     fileLines = csv.DictReader(theFile,delimiter='|')
    fileLines = csv.DictReader(theFile)

    
    nameList = []
    coordList = []
    
#     for i in fileLines:
#         name = i['name']
#         ra = i['ra']
#         dec = i['dec']

    for i in fileLines:
        name = i['Object Name']
        ra = i['RA(deg)']
        dec = i['DEC(deg)']
                
        nameList.append(name)
        coordList.append([ra,dec])
        
    theFile.close()
    
    # format the names in the list and return a new list
    newNameList = parseGalaxyNames(nameList)

    
    fieldnames = ('preferredName',\
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
                lineList = [line['preferredName'],\
                line['oldName'],\
                line['redshift'],\
                line['degreesJ2000RA_Dec'],\
                line['J2000RA_Dec'],\
                line['galacticLong_Lat'],\
                line['rIndependentDistMean_sd_min_max (Mpc)'],\
                line['morphology'],\
                line['distanceIndicator'],\
                line['luminosityClass'],\
                line['EBminusV'],\
                line['radialVelocity (km/s)'],\
                line['vcorr (km/s)'],\
                line['angDiameters (arcsec)'],\
                line['linDiameters (kpc)'],\
                line['distvcorr (Mpc)'],\
                line['inclination (deg)'],\
                line['photometry'],\
                line['alternativeNames']]
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
                    coordList = coordList[nameListIndex+1:]
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
                # create a new file based on the newly entered file name
                opts.outname = aTwo
                writerOutFile = open(opts.output+opts.outname+'.csv','wt')
                writer = createCSVTable(writerOutFile,fieldnames)

    else:
        # It does not already exist. Create it as a dictionary csv file
        writerOutFile = open(opts.output+opts.outname+'.csv','wt')
        writer = createCSVTable(writerOutFile,fieldnames)
    
    
    queue = Queue.Queue()
    outQueue = Queue.Queue()
    numberQueue = Queue.Queue()
    
    t1 = ThreadUrl(queue,outQueue,numberQueue)
    t2 = ThreadUrl(queue,outQueue,numberQueue)
    t3 = ThreadUrl(queue,outQueue,numberQueue)
    t1.setDaemon(True)
    t2.setDaemon(True)
    t3.setDaemon(True)
    t1.start()
    t2.start()
    t3.start()
        
    total = len(newNameList)
    counter = maxRetrieve
    for coord,name in zip(coordList,newNameList):
        
#         if len(n) >= 14:
#             # grab the name for this galaxy
#             ra, dec = coord
#             searchRadius = 0.1
#             name = returnGalaxyName(ra,dec,searchRadius)
#         else:
#             name = n
    
        # return basic information for each object
        
        totalStart = time.time()
        percentComplete = round(float((maxRetrieve-counter))/total *100,1)
        sys.stdout.write("\n")
        sys.stdout.write("Percent complete: %s %% \n" %percentComplete)
        sys.stdout.write("Galaxies left: %s \n" %counter)
        sys.stdout.write("\n")
        sys.stdout.write("Starting: %s \n" %name)

        start2 = time.time()
        hosts = [["http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname=%s"\
        %name,1],["http://ned.ipac.caltech.edu/cgi-bin/nDistance?name=%s"\
        %name,2],["http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES"\
        %name,3]]

        #populate queue with data
        for host in hosts:
            time.sleep(0.9)
            queue.put(host)
    
        #wait on the queue until everything has been processed
        queue.join()
        
        # decide which object is which in the queue
        first = outQueue.get()
        second = outQueue.get()
        third = outQueue.get()
        
        number1 = numberQueue.get()
        number2 = numberQueue.get()
        number3 = numberQueue.get()
        
        if number1 == 1:
            votable = first
        elif number1 == 2:
            html = first
        elif number1 == 3:
            fullHtml = first
            
        if number2 == 1:
            votable = second
        elif number2 == 2:
            html = second
        elif number2 == 3:
            fullHtml = second
            
        if number3 == 1:
            votable = third
        elif number3 == 2:
            html = third
        elif number3 == 3:
            fullHtml = third
        
        queryTime = time.time() - start2

        print 'queue is empty: ',queue.empty()
        print 'outqueue is empty: ',outQueue.empty()
        print 'numberqueue is empty: ',numberQueue.empty()

        galaxy = galaxyClass(name,votable,html,fullHtml)
#       galaxy.queryNED()
        sys.stdout.write("1...")
        sys.stdout.flush()
        
        start3 = time.time()
        redshift = galaxy.returnRedshift()
        sys.stdout.write("2...")
        sys.stdout.flush()
        
        degreesJ2000RA,degreesJ2000Dec = galaxy.returnJ2000Position()
        if isNumber(degreesJ2000RA) and isNumber(degreesJ2000Dec):
            J2000RA_Dec = convertRAandDec(degreesJ2000RA,degreesJ2000Dec,'sexagesimal')
#             J2000RA_Dec = returnRAandDec(degreesJ2000RA,degreesJ2000Dec)
        else:
            J2000RA_Dec = ('x','x')
            degreesJ2000RA,degreesJ2000Dec = coord
        sys.stdout.write("3...")
        sys.stdout.flush()      
        
        galacticLong_Lat = galaxy.returnGalactic()
        sys.stdout.write("4...")
        sys.stdout.flush()          
        
        rIndependentDistMean_sd_min_max = galaxy.returnRedIndependentDist()
        sys.stdout.write("5...")
        sys.stdout.flush()  
        
        morphology = galaxy.returnMorphology()
        sys.stdout.write("6...")
        sys.stdout.flush()  
        
        distanceIndicator = galaxy.returnDistIndicator()
        sys.stdout.write("7...")
        sys.stdout.flush()
        
        luminosityClass = galaxy.returnLuminosityClass()
        
        EBminusV = galaxy.returnEBminusV()
        
        radialVelocity = galaxy.returnRadialVelocity()
        if isNumber(radialVelocity):
            if float(radialVelocity) >=0:
                vcorr = calculatevcorr(degreesJ2000RA,degreesJ2000Dec,radialVelocity)
                if vcorr >0:
                    distvcorr = vcorr / hubbleConstant
                else:
                    distvcorr ='x'
            else:
                vcorr = 'x'
                distvcorr = 'x'
        else:
            vcorr = 'x'
            distvcorr = 'x'
        
        
        # I am assuming that diameters are found in arcmin, and I'm converting to arcsec
        diameters = galaxy.returnDiameters()
        angMaj = diameters[0]
        if isNumber(angMaj):
            angMaj = float(angMaj)*60
        angMin = diameters[1]
        if isNumber(angMin):
            angMin = float(angMin)*60
            # switch them if the major diameter is larger than the minor
            if isNumber(angMaj):
                if angMaj < angMin:
                    angMaj,angMin = angMin, angMaj
        else:
            angMin = 'x'
            
        if isNumber(angMaj) and isNumber(angMin) and isNumber(distvcorr):
#             linDiameters = returnLinDiameters(distvcorr,angMaj,angMin)
            linDiameters = calculateLinearDiameters(angMaj,angMin,distvcorr)
        else:
            linDiameters = ('x','x')
        sys.stdout.write("8...")
        sys.stdout.flush()
        
        if isNumber(angMaj) and isNumber(angMin):
#             inclination = returnInclination(angMaj,angMin)
            inclination = calculateInclination(angMaj,angMin)

        else:
            inclination = 'x'
        
#       photometry = galaxy.returnPhotometry()
        photometry = ['x']
        sys.stdout.write("9...")
        sys.stdout.flush()
        
        alternativeNames = galaxy.returnNames()
        sys.stdout.write("10")
        sys.stdout.flush()
        
        oldName = urllib.unquote_plus(name).replace('\n','').strip()
#         preferredName = decidePreferredName(alternativeNames,oldName)
        preferredName = pickPreferredName3(alternativeNames,oldName)

        
        strippedAlternativeNames = []
        for alternate in alternativeNames:
            stripped = alternate.replace(' ','')
            strippedAlternativeNames.append(stripped)
            
        print 'query time: ',queryTime
        print 'time for other stuff: ',time.time() - start3
        
        objectInfoList = [preferredName.replace(' ',''),\
        oldName.replace(' ',''),\
        redshift,\
        (float(degreesJ2000RA),float(degreesJ2000Dec)),\
        J2000RA_Dec,galacticLong_Lat,\
        rIndependentDistMean_sd_min_max,\
        morphology,\
        distanceIndicator,\
        luminosityClass,\
        EBminusV,\
        radialVelocity,\
        vcorr,\
        (angMaj,angMin),\
        linDiameters,\
        distvcorr,\
        inclination,\
        photometry,\
        strippedAlternativeNames]
                
        # write info to file
        row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
        print 'type, row: ',type(row),', ',row
        if opts.verbose:
            print '\n','row: ',row
            
        writer.writerow(row)
        sys.stdout.write("\n")
        print 'Total retrieval time: ',time.time() - totalStart
        counter -=1
        if counter == 0:
            print 'Reached maxRetrieve = {0}. Exiting...'.format(maxRetrieve)
            break

        
    writerOutFile.close()
    if opts.verbose:
        print "Done."
        print 
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