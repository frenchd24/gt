#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_RIDplot.py, v1 10/31/17

Makes a plot of RID*HubbleC - vcorr across the sky. 

MUST BE RUN INSIDE THE PYSALT CONDA ENVIRONMENT:
>>> source activate pysalt
>>> python GTupdate2_RIDplot.py
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
        inFilename = '/usr/data/moosejaw/frenchd/GT_update2/FinalGalaxyTable190_filtered.csv'

    elif user =='David':
        inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_filtered.csv'

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
    
    RID_dif_means = []
    RID_dif_medians = []
    ras = []
    decs = []

    # these are for the outlier cases
    RID_dif_means_2 = []
    RID_dif_medians_2 = []
    ras_2 = []
    decs_2 = []
    
    
    cutoff = 3000
    
    # do the work
    count = 0
    for l in reader:
        count +=1
        
        Name = l['Name']
        NEDname = l['NEDname']
        vcorr = eval(l['vcorr'])
        RID_mean = eval(l['RID_mean'])
        RID_median = eval(l['RID_median'])
        flag = l['flag']

        ra = eval(l['RAdeg'])
        dec = eval(l['DEdeg'])

        
        if RID_mean != nullFloat:
            rmeanv = float(RID_mean) * hubbleC
            rmedianv = float(RID_median) * hubbleC
        
            RID_dif_mean = rmeanv-vcorr
            RID_dif_median = rmedianv-vcorr
            
#             if abs(RID_dif_mean) >= 10000:
#                 RID_dif_means.append(10000)
#                 print 'Name = {0}, RID_mean = {1}, vcorr = {2}'.format(Name,RID_mean,vcorr}
            
#             elif abs(RID_dif_mean) >= 7500:
#                 RID_dif_means.append(7500)
# 
#             elif abs(RID_dif_mean) >= 5000:
#                 RID_dif_means.append(5000)
#                 
#             elif abs(RID_dif_mean) >= 2500:
#                 RID_dif_means.append(2500)
#                 
#             elif abs(RID_dif_mean) >= 1500:
#                 RID_dif_means.append(1500)
#             
#             else:
#                 RID_dif_means.append(RID_dif_mean)

            if abs(RID_dif_mean) <= cutoff:
#                 RID_dif_means.append(RID_dif_mean)
#                 RID_dif_medians.append(RID_dif_median)
                RID_dif_means.append(abs(RID_dif_mean))
                RID_dif_medians.append(abs(RID_dif_median))
                ras.append(ra)
                decs.append(dec)

            else:
#                 RID_dif_means_2.append(RID_dif_mean)
#                 RID_dif_medians_2.append(RID_dif_median)
                RID_dif_means_2.append(abs(RID_dif_mean))
                RID_dif_medians_2.append(abs(RID_dif_median))
                ras_2.append(ra)
                decs_2.append(dec)
                
#             if abs(RID_dif_mean) >= 5000 and abs(RID_dif_mean) <= 10000:
#                 RID_dif_means.append(RID_dif_mean)
#                 RID_dif_medians.append(RID_dif_median)
#                 ras.append(ra)
#                 decs.append(dec)
        
        # update the counter
        percentComplete = round((float(count)/130759)*100,2)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        sys.stdout.flush()
            
##########################################################################################
    # plot it
    colmap = cm.RdBu_r
    
    vmaxVal = max(RID_dif_means)
    vminVal = min(RID_dif_means)
    
    ras = numpy.array(ras)
    decs = numpy.array(decs)
    colors = numpy.array(RID_dif_means)
    
    ras_2 = numpy.array(ras_2)
    decs_2 = numpy.array(decs_2)
    colors_2 = numpy.array(RID_dif_means_2)
    
    
#     print
#     print 'ra, dec : ',thing.transpose()
#     print

    
    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
#     norm = matplotlib.colors.Normalize(vmin = vmaxVal, vmax = vminVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)


#     fig = figure(figsize=(9,7))
#     ax = fig.add_subplot(111,projection='aitoff')
    
#     plt.subplot(111, projection="aitoff")
#     fig, ax1 = plt.subplot(111, projection="aitoff")
#     fig, (ax1, ax2) = plt.subplots(figsize=(8, 3), ncols=2)
#     plot1 = ax1.scatter(ras,decs,s=10,c=RID_dif_means,vmin=vminVal,vmax=vmaxVal,marker='.',cmap=colmap)
#     plot1 = scatter(px,py,s=30,c=RID_dif_means,vmin=vminVal,vmax=vmaxVal,marker='.',lw=0,cmap=colmap)

#     plot1 = ax.scatter(x,y,s=30,c=colors,vmin=vminVal,vmax=vmaxVal,marker='.',lw=0,cmap=colmap)
#     plt.grid(True)
#     cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical')
#     cbar.set_label(r'$\rm RID_{mean} -vcorr ~[km ~s^{-1}]$')
#     plt.title("Aitoff - RID_dif_means")
#     plt.show()


    # Import all required packages.
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import matplotlib.pyplot as plt
    import matplotlib.pylab as lab
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable


    # Generate random data, for RA between 0 and 360 degrees, for DEC between
    # -90 and +90 degrees.
#     ra_random = np.random.rand(100)*360.0
#     dec_random = np.random.rand(100)*180.0-90.0
    ra_random = ras
    dec_random = decs
    
    ra_random2 = ras_2
    dec_random2 = decs_2

    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ra_random*u.degree, dec=dec_random*u.degree, frame='icrs')
    c2 = SkyCoord(ra=ra_random2*u.degree, dec=dec_random2*u.degree, frame='icrs')


    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi
    
    ra_rad2 = c2.ra.radian
    dec_rad2 = c2.dec.radian
    ra_rad2[ra_rad2 > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

#     lab.title("Aitoff projection of our random data")
#     lab.grid(True)
#     plot1 = plt.plot(ra_rad, dec_rad,'o',c=colors,lw=0,cmap=colmap)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,vmax=vmaxVal,lw=0,cmap=colmap,s=10)
#     plot2 = plt.scatter(ra_rad2, dec_rad2,marker='o',c='green',lw=0,s=10)
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.046, pad=0.04)
    cbar.set_label(r'$\rm RID_{mean} -vcorr ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=500)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    
#     plt.show()
    saveDirectory = '/Users/David/Research_Documents/GT_update2/'
    plt.savefig('{0}RID-vcorr_{1}cut_nogreen_abs.pdf'.format(saveDirectory,cutoff),format='pdf')


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
