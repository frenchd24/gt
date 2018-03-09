#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_mag_v_diam.py, v1 01/08/18

Plot magnitude vs diameter for the whole galaxy table


Based on: GTupdate2_allSkyPlot.py, v1 10/31/17

Makes all sky plots of a bunch of stuff 

MUST BE RUN INSIDE THE PYSALT CONDA ENVIRONMENT:
>>> source activate pysalt
>>> python GTupdate2_allSkyPlot.py
>>> ...sweet success...

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
from pylab import *

import math

import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

from astropy import wcs
from astropy.io import fits

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################

def calculate_absoluteMag(m,dm,d,dd,e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    #
    # dm is the error in apparent magnitude, dd is the error in distance
    
    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)
    
    # now do the error
    # error is dominated by dd, so don't even bother with an extinction error

    dM = math.sqrt(float(dm)**2 + ((5*float(dd))**2 / (float(d) * math.log10(10))**2))
    
    return M, dM



def median_low(l):
    # returns the closest element in a list to the median, rounding down

    # E.g.,
    # list = [1, 2, 3, 4, 5]
    # median_low(list) = 3
    #
    # list = [1, 2, 3, 4]
    # median_low(list) = 2
    
#     l.sort()
    l = np.array(l)
    med = np.median(l)
    
    diff = abs(np.subtract(med,l))
    
    diff = list(diff)
    l = list(l)
    indexMin = diff.index(min(diff))
    
    return l[indexMin]

    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':        
        inFilename = '/Users/frenchd/Research/gt/FinalGalaxyTable12_filtered.csv'

    elif user =='David':
        pass
#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_filtered.csv'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
        
    # open the files
        # open the files
    csv.field_size_limit(sys.maxsize)
    
    # read in the tables
    inFile = open(inFilename,'rU')
    reader = csv.DictReader(inFile)
    
    # what are the null values equal to?
    nullFloat = -99.99
    nullInt = -99
    nullStr = 'x'
    
    # the header/column names
    fieldnames = (\
    'Name',\
    'NEDname',\
    'z',\
    'RAdeg',\
    'DEdeg',\
    'RAh',\
    'RAm',\
    'RAs',\
    'DE-',\
    'DEd',\
    'DEm',\
    'DEs',\
    'GLON',\
    'GLAT',\
    'Vhel',\
    'vcorr',\
    'distvcorr',\
    'RID_mean',\
    'RID_median',\
    'RID_std',\
    'RID_min',\
    'RID_max',\
    'bestDist',\
    'e_bestDist',\
    'MajDiam_ang',\
    'MinDiam_ang',\
    'e_MajDiam_ang',\
    'e_MinDiam_ang',\
    'MajDiam',\
    'MinDiam',\
    'e_MajDiam',\
    'e_MinDiam',\
    'R_vir',\
    'inc',\
    'adjustedInc',\
    'e_inc',\
    'PA',\
    'diameterKey',\
    'ratioKey',\
    'paKey',\
    'RC3_type',\
    'RC3_d25',\
    'RC3_r25',\
    'RC3_pa',\
    'group_num',\
    'group_mem',\
    'group_dist',\
    'MType',\
    'flag',\
    'distIndicator',\
    'lumClass',\
    'E(B-V)',\
    'Bmag',\
    'Bmag_key',\
    'Bmag_max',\
    'Bmag_max_key',\
    'Bmag_min',\
    'Bmag_min_key',\
    'Bmag_sdss',\
    'gmag_sdss',\
    'rmag_sdss',
    'zmag_sdss',\
    'Lstar_med',\
    'e_Lstar_med',\
    'Lstar_max',\
    'e_Lstar_max',\
    'Lstar_min',\
    'e_Lstar_min',\
    'Lstar_sdss',\
    'e_Lstar_sdss',\
    'altNames')
    
    
##########################################################################################
##########################################################################################    

    hubbleC = 71.0
    
    vcorrs = []
    vhels = []
    ras = []
    decs = []

    BmagL = []
    MajDiamL = []
        
    # do the work
    count = 0
    for l in reader:
        count +=1
        
        Name = l['Name']
        NEDname = l['NEDname']
        vcorr = eval(l['vcorr'])
        RID_mean = eval(l['RID_mean'])
        RID_median = eval(l['RID_median'])
        bestDist = l['bestDist']
        ebminusv = l['E(B-V)']
        Vhel = eval(l['Vhel'])
        flag = l['flag']
        Bmag = l['Bmag']
        MajDiam = l['MajDiam']

        ra = eval(l['RAdeg'])
        dec = eval(l['DEdeg'])
        
        if str(Bmag) != str(nullFloat) and str(MajDiam) != str(nullFloat):
            
            if str(bestDist) != str(nullFloat):
                Bmag_abs = calculate_absoluteMag(float(Bmag),0.01,bestDist,1.,ebminusv)[0]
            
                vcorrs.append(vcorr)
                vhels.append(Vhel)
                ras.append(ra)
                decs.append(dec)
                BmagL.append(float(Bmag_abs))
                MajDiamL.append(float(MajDiam))

        # update the counter
        percentComplete = round((float(count)/130759)*100,2)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        sys.stdout.flush()
            
##########################################################################################
    # plot it
    colmap = cm.RdBu_r
    
    ras = numpy.array(ras)
    decs = numpy.array(decs)
#     colors = numpy.array(vhels)
#     
#     vmaxVal = max(colors)
#     vminVal = min(colors)
# 
#     norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
#     m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)


    # Import all required packages.
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import matplotlib.pyplot as plt
    import matplotlib.pylab as lab
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable


    # Generate random data, for RA between 0 and 360 degrees, for DEC between
    # -90 and +90 degrees.
#     ra_random = ras
#     dec_random = decs

    # Transform the data into a SkyCoord object.
#     c = SkyCoord(ra=ra_random*u.degree, dec=dec_random*u.degree, frame='icrs')


    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
#     ra_rad = c.ra.radian
#     dec_rad = c.dec.radian
#     ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
#     ax = lab.subplot(111,projection="aitoff")
    ax = lab.subplot(111)

    # lab.subplot(111)
    
    plot1 = plt.scatter(MajDiamL,BmagL,marker='o',lw=0,s=4,alpha=0.5)
    xlabel(r'$\rm MajDiam ~ [kpc]$')
    ylabel(r'$\rm M_B$')
    ylim(-6,-26)
    ax.set_xscale('log')
    xlim(1,90)


#     plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
#     vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=0.5)
#     
#     cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
#     cbar.set_label(r'$\rm Vhel ~[km ~s^{-1}]$')
#     
#     xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
#     ax.set_xticklabels(xlab, weight=510)
#     ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    saveDirectory = '/Users/frenchd/Research/GT_update2/'
    plt.savefig('{0}mag_v_diam.pdf'.format(saveDirectory),format='pdf')


##########################################################################################
##########################################################################################
    
    # close the files
    inFile.close()
    
    print "Done."
    print
    
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
