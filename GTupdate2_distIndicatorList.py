#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_distIndicatorList.py, v 1 09/28/2017

List all the distIndicators in FinalGalaxyTable7.csv

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


def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':
        filename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable7.csv'
    
    elif user =='David':
        filename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable7.csv'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the files
    csv.field_size_limit(sys.maxsize)
    
    # basic info file
    basicFile = open(filename,'rU')
    basicReader = csv.DictReader(basicFile)
  
    d_rc3 = {}
    d_morph = {}
    d_lum = {}
    d_dist = {}
    
    excludeList = [\
    'M-star',\
    'M_star',\
    'Opt.var.',\
    'K4-K5;Candidate_WD',\
    'F6-F8;Candidate_WD',\
    'A',\
    'E',\
    'Candidate_AGN',\
    'M1',\
    'star??',\
    'O',\
    'K_Star',\
    'PN?',\
    'K1',\
    'M0',\
    'M0V',\
    'A0',\
    'DA-star',\
    'High_vel._cloud',\
    'O',\
    'Carbon',\
    'Point_Src_[SDSS]',\
    'Possible_star',\
    'Planetary_nebula',\
    'M3-M4',\
    'F2',\
    'A-star'\
    'PN:',\
    'Cand._glob._cluster',\
    'Candidate_PN',\
    'F']
    
    okList = ['pec','SBc','Sab','SA0','SB(','SA(','SB',\
    'SA',\
    'Sb',\
    'Sc',\
    'S0',\
    'Irr',\
    'dSph',\
    'LINER',\
    'Dwarf',\
    'E5',\
    'Im',\
    'E0',\
    'E+',\
    'BCD',\
    'Sbrst',\
    'Elliptical',\
    'Spiral',\
    'Im_HII',\
    'dwarf',\
    'BLLAC',\
    'Extended_Src_[SDSS]',\
    'dE4',\
    'dE2',\
    'Sd',\
    'IBm',\
    'ELG',\
    'Ir',\
    'dE7',\
    'Sm',\
    'cD',\
    'Sa',\
    'S...',\
    'AGN']
      
    count = 0
    # loop through and do the work
    for line in basicReader:
        # update the counter
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
        
        # stuff in basicReader:
        Name = line['Name']
        distIndicator = str(line['distIndicator'])
        lumClass = str(line['lumClass'])
        MType = str(line['MType'])
        RC3_type = str(line['RC3_type'])
            
        if d_rc3.has_key(RC3_type):
            d_rc3[RC3_type] +=1
        else:
            d_rc3[RC3_type]=1
            
#         if bfind(MType,'WR'):
#             print 'MType = ',MType
#             print 'Name = ',Name
        
        findi = False
        for i in okList:
            if MType.find(i) != -1:
                findi = True
                
        if not findi and MType != 'x':
#             print 'MType = ',MType
#             print 'Name = ',Name
            
#         if len(MType) < 4:
#             if d_morph.has_key(MType):
#                 d_morph[MType] +=1 
#             else:
#                 d_morph[MType]=1      
            
            if d_morph.has_key(MType):
                d_morph[MType] +=1 
            else:
                d_morph[MType]=1
            
            
        if d_dist.has_key(distIndicator):
            d_dist[distIndicator] +=1 
        else:
            d_dist[distIndicator]=1

        if d_lum.has_key(lumClass):
            d_lum[lumClass] +=1 
        else:
            d_lum[lumClass]=1

            
    print 'Results - RC3: ',
    for i in d_rc3:
        print i
   
    print
    print 'Results - morph: ',
    for i in d_morph:
        print i

    print
    print 'Results - distIndicator: ',
    for i in d_dist:
        print i
        
    print
    print 'Results - lumClass: ',
    for i in d_lum:
        print i
        

    basicFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
