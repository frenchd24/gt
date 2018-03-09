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
        pickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT.p'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()

    
    theFile = open(filename,'rU')
    pickleFile = open(pickleFilename,'wt')
    reader = csv.DictReader(theFile)

    
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
    
    for l in reader:
        majorAxis = l['MajDiam']
        inc = l['inc']
        adjustedInc = l['adjustedInc']
        pa = l['PA']
        Bmag = l['Bmag']
        Bmag_sdss = l['Bmag_sdss']
        RID_median = l['RID_median']
        RID_mean = l['RID_mean']
        RID_std = l['RID_std']
        Vhel = l['Vhel']
        RAdeg = l['RAdeg']
        DEdeg = l['DEdeg']
        Name = l['Name']

        majorAxisL.append(majorAxis)
        incL.append(inc)
        adjustedIncL.append(adjustedInc)
        paL.append(pa)
        BmagL.append(Bmag)
        Bmag_sdssL.append(Bmag_sdss)
        RID_medianL.append(RID_median)
        RID_meanL.append(RID_mean)
        RID_stdL.append(RID_std)
        VhelL.append(Vhel)
        RAdegL.append(RAdeg)
        DEdegL.append(DEdeg)
        NameL.append(Name)
    
    gtDict['majorAxis'] = majorAxisL
    gtDict['inc'] = incL
    gtDict['adjustedInc'] = adjustedIncL
    gtDict['PA'] = paL
    gtDict['Bmag'] = BmagL
    gtDict['Bmag_sdss'] = Bmag_sdssL
    gtDict['RID_median'] = RID_medianL
    gtDict['RID_mean'] = RID_meanL
    gtDict['RID_std'] = RID_stdL
    gtDict['Vhel'] = VhelL
    gtDict['RAdeg'] = RAdegL
    gtDict['DEdeg'] = DEdegL
    gtDict['Name'] = NameL
    
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
    