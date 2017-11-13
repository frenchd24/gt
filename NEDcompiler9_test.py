#!/usr/bin/env python

'''
NEDcompiler9_test.py 09/06/17

Just stripping everything but the core away to see why some galaxies won't return right

--- from NEDcompiler9.py ---

A program to grab galaxy data from the NED server by means of xml-based VOTable
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
    
    
class galaxyClass(object):
    # a galaxy object containing it's information and methods of retrieving and 
    # returning it
    
    def __init__(self,name,votable,html,fullHtml):
        self.name = name
        self.votable = votable
        self.html = html
        self.fullHtml = fullHtml

        
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

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    hubbleConstant = 71
    
    # set a limit to how many queries to make?
    maxRetrieve = 10000
    
    galaxyName = '2MASX+J18363718-3218138'
    ra = '279.15488d'
    dec = '-32.30381d'
    radius = 0.1
    
    print 'galaxyName: ',returnGalaxyName(ra,dec,radius)
    print
        
    host1 = "http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname={0}".format(galaxyName)
    host2 ="http://ned.ipac.caltech.edu/cgi-bin/nDistance?name={0}".format(galaxyName)
    host3 = "http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?objname={0}&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES".format(galaxyName)


    #grabs urls of hosts and reads or parses them according to their type
    url = 'x'
    try:
        url = urllib2.urlopen(host1)
        print 'url: ',url
        print
    except Exception, e:
        sys.stderr.write("\n Unable to return url file object. Here is the error message built into the exception: \n %s\n" %e)
        
    warnings.simplefilter("ignore")
    votable = 'x'
    try:
        votable = parse(url,pedantic=False)
        print 'votable1: ',votable1
        print
    except Exception, e:
        sys.stderr.write("\n Unable to parse voTable. Here is the message built into the exception: \n %s \n" %e)
    
    warnings.resetwarnings()

    html = 'x'
    try:
        html = url.readlines()
        print 'html: ',html
        print
    except Exception, e:
        sys.stderr.write("\n Unable to read html. Here is the message built into the exception: \n %s \n" %e)
        

    redshift = 'x'
    try:
        warnings.simplefilter("ignore")
        table = votable.get_table_by_id("NED_MainTable")
        redshift = table.array['Redshift'][0]
        warnings.resetwarnings()
        
        print 'table: ',table
        print
        print 'redshift: ',redshift
        
    except Exception, e:
        sys.stderr.write("\n Unable to return redshift. Here is the error message "\
        "built into the exception:\n %s\n" %e)

    url.close()

###############################################################################

if __name__=="__main__":
    # parse commandline
    # do the work
    main()