#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: pickleGT.py v 1.0 01/03/2018

Pickle the full galaxy table


Comes from: buildDataLists3.py, v 5.0 12/02/2015

Makes a pickle file with all the info from:
salt_sightlines_all_results_include.csv



'''

import sys
import os
import csv

from pylab import *
# import atpy
from math import *
from utilities import *
import getpass
import pickle

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from matplotlib import rc


###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    if getpass.getuser() == 'frenchd':
        filename = '/Users/frenchd/Research/gt/FinalGalaxyTable13_filtered.csv'
        pickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()

    
    theFile = open(filename,'rU')
    pickleFile = open(pickleFilename,'wt')
    reader = csv.DictReader(theFile)

    nullFloat = -99.99
    nullStr = 'x'
    nullInt = -99
    
    # overall structure: gtDict is a dictionary containing a bunch of shit from the full
    # galaxy table
    gtDict = {}

    
    # lists for stuff in galaxy table
    majorAxisL = []
    incL = []
    adjustedIncL = []
    paL = []
    BmagL = []
    Bmag_sdssL = []
    RID_medianL = []
    RID_meanL = []
    RID_stdL = []
    VhelL = []
    RAdegL = []
    DEdegL = []
    NameL = []
    LstarL = []
    bestDistL = []
    vcorrL = []
    
    for l in reader:
        majorAxis = l['MajDiam']
        inc = l['inc']
        adjustedInc = l['adjustedInc']
        pa = l['PA']
        Bmag = float(l['Bmag'])
        Bmag_sdss = float(l['Bmag_sdss'])
        RID_median = l['RID_median']
        RID_mean = l['RID_mean']
        RID_std = l['RID_std']
        Vhel = l['Vhel']
        RAdeg = l['RAdeg']
        DEdeg = l['DEdeg']
        Name = l['Name']
        Lstar_med = float(l['Lstar_med'])
        Lstar_sdss = float(l['Lstar_sdss'])
        bestDist = l['bestDist']
        vcorr = l['vcorr']
        flag = int(l['flag'])
        
#         if flag !=1:
        if flag ==0:

            # use Lstar_sdss if needed
            lstar = -99.99
            if Lstar_med > 0:
                lstar = Lstar_med
            elif Lstar_med < 0 and Lstar_sdss > 0:
                lstar = Lstar_sdss
            else:
                lstar = Lstar_med
            
            # use Bmag_sdss if needed
            b = -99.99
            if Bmag > 0:
                b = Bmag
            elif Bmag < 0 and Bmag_sdss > 0:
                b = Bmag_sdss
            else:
                b = Bmag

                
            if b > 0 and lstar < 0:
                print 'Name: ',Name, ' : ',b,', ',lstar
                print

            majorAxisL.append(float(majorAxis))
            incL.append(int(inc))
            adjustedIncL.append(int(adjustedInc))
            paL.append(int(pa))
            BmagL.append(float(b))
            RID_medianL.append(float(RID_median))
            RID_meanL.append(float(RID_mean))
            RID_stdL.append(float(RID_std))
            VhelL.append(int(Vhel))
            RAdegL.append(float(RAdeg))
            DEdegL.append(float(DEdeg))
            NameL.append(Name)
            LstarL.append(float(lstar))
            bestDistL.append(float(bestDist))
            vcorrL.append(int(vcorr))
        
    
    print 'len Lstar: ',len(LstarL)
    print
    gtDict['majorAxis'] = majorAxisL
    gtDict['inc'] = incL
    gtDict['adjustedInc'] = adjustedIncL
    gtDict['PA'] = paL
    gtDict['Bmag'] = BmagL
    gtDict['RID_median'] = RID_medianL
    gtDict['RID_mean'] = RID_meanL
    gtDict['RID_std'] = RID_stdL
    gtDict['Vhel'] = VhelL
    gtDict['RAdeg'] = RAdegL
    gtDict['DEdeg'] = DEdegL
    gtDict['Name'] = NameL
    gtDict['Lstar'] = LstarL
    gtDict['bestDist'] = bestDistL
    gtDict['vcorr'] = vcorrL
    
    pickle.dump(gtDict,pickleFile)
    pickleFile.close()
    theFile.close()
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    