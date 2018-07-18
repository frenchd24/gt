#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_allSkyPlot2.py, v2 06/14/18

v2: Make 10 separate plots, one for each 1000 km/s velocity step

Makes all sky plots of a bunch of stuff 

MUST BE RUN INSIDE THE PYSALT CONDA ENVIRONMENT:
>>> source activate pysalt
>>> python GTupdate2_allSkyPlot.py
>>> ...sweet success...

'''

# Import all required packages.
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.pylab as lab
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable



import sys
import os
import csv
# import string
import warnings
import numpy as np
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
        inFilename = '/Users/frenchd/Research/gt/FinalGalaxyTable13_filtered.csv'
        saveDirectory = '/Users/frenchd/Research/GT_update2/galaxyTable_paper/figures/all_sky_12h/'

    elif user =='David':
        pass
#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_filtered.csv'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    
    # plot all sky map in 10 velocity bins centered on 0h RA
    plot_all_sky = False
    
    # plot all sky map in 10 velocity bins centered on 12h RA
    plot_all_sky_12h = True
        
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
    
    all_vcorrs = []
    all_vhels = []
    all_ras = []
    all_decs = []
    
    
    v1 = []
    ra1 = []
    dec1 = []

    v2 = []
    ra2 = []
    dec2 = []

    v3 = []
    ra3 = []
    dec3 = []
        
    v4 = []
    ra4 = []
    dec4 = []

    v5 = []
    ra5 = []
    dec5 = []
        
    v6 = []
    ra6 = []
    dec6 = []

    v7 = []
    ra7 = []
    dec7 = []

    v8 = []
    ra8 = []
    dec8 = []
    
    v9 = []
    ra9 = []
    dec9 = []

    v10 = []
    ra10 = []
    dec10 = []
    
            
    # do the work
    count = 0
    for l in reader:
        count +=1
        
        Name = l['Name']
        NEDname = l['NEDname']
        vcorr = eval(l['vcorr'])
        RID_mean = eval(l['RID_mean'])
        RID_median = eval(l['RID_median'])
        Vhel = eval(l['Vhel'])
        flag = eval(l['flag'])

        ra = eval(l['RAdeg'])
        dec = eval(l['DEdeg'])
        
#         if float(ra) <10. and float(ra) > -10:

        if flag == 0:
            go = True
        else:
            go = False
            
        if Vhel <= 450:
            go = False
            
        if go:
            all_vcorrs.append(vcorr)
            all_vhels.append(Vhel)
            all_ras.append(ra)
            all_decs.append(dec)
        
            if Vhel <= 1000.:
                v1.append(Vhel)
                ra1.append(ra)
                dec1.append(dec)
        
            if Vhel <= 2000. and Vhel >1000:
                v2.append(Vhel)
                ra2.append(ra)
                dec2.append(dec)

            if Vhel <= 3000. and Vhel >2000:
                v3.append(Vhel)
                ra3.append(ra)
                dec3.append(dec)

            if Vhel <= 4000. and Vhel >3000:
                v4.append(Vhel)
                ra4.append(ra)
                dec4.append(dec)
            
            if Vhel <= 5000. and Vhel >4000:
                v5.append(Vhel)
                ra5.append(ra)
                dec5.append(dec)
            
            if Vhel <= 6000. and Vhel >5000:
                v6.append(Vhel)
                ra6.append(ra)
                dec6.append(dec)
            
            if Vhel <= 7000. and Vhel >6000:
                v7.append(Vhel)
                ra7.append(ra)
                dec7.append(dec)
            
            if Vhel <= 8000. and Vhel >7000:
                v8.append(Vhel)
                ra8.append(ra)
                dec8.append(dec)
            
            if Vhel <= 9000. and Vhel >8000:
                v9.append(Vhel)
                ra9.append(ra)
                dec9.append(dec)

            if Vhel <= 10000. and Vhel >9000:
                v10.append(Vhel)
                ra10.append(ra)
                dec10.append(dec)
            

        # update the counter
        percentComplete = round((float(count)/130759)*100,2)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        sys.stdout.flush()
            
##########################################################################################
##########################################################################################

    alpha_galaxy = 0.6
    size_galaxy = 2

    # plot it
#     colmap = cm.RdBu_r
#     colmap = cm.cool
    colmap = cm.plasma

    colors = np.array(all_vhels)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)

##########################################################################################
##########################################################################################

    if plot_all_sky:

    ##########################################################################################
        # 0 - 1000 km/s
    
        ras = np.array(ra1)
        decs = np.array(dec1)
        colors = np.array(v1)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}1.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 1000 - 2000 km/s
    
        ras = np.array(ra2)
        decs = np.array(dec2)
        colors = np.array(v2)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=4, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}2.pdf'.format(saveDirectory),format='pdf')
    
    ##########################################################################################
        # 2000 - 3000 km/s
    
        ras = np.array(ra3)
        decs = np.array(dec3)
        colors = np.array(v3)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}3.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 3000 - 4000 km/s
    
        ras = np.array(ra4)
        decs = np.array(dec4)
        colors = np.array(v4)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}4.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 4000 - 5000 km/s
    
        ras = np.array(ra5)
        decs = np.array(dec5)
        colors = np.array(v5)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}5.pdf'.format(saveDirectory),format='pdf')



    ##########################################################################################
        # 5000 - 6000 km/s
    
        ras = np.array(ra6)
        decs = np.array(dec6)
        colors = np.array(v6)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}6.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 6000 - 7000 km/s
    
        ras = np.array(ra7)
        decs = np.array(dec7)
        colors = np.array(v7)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}7.pdf'.format(saveDirectory),format='pdf')
    
    
    
    ##########################################################################################
        # 7000 - 8000 km/s
    
        ras = np.array(ra8)
        decs = np.array(dec8)
        colors = np.array(v8)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}8.pdf'.format(saveDirectory),format='pdf')
    
    
    ##########################################################################################
        # 8000 - 9000 km/s
    
        ras = np.array(ra9)
        decs = np.array(dec9)
        colors = np.array(v9)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}9.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 9000 - 10000 km/s
    
        ras = np.array(ra10)
        decs = np.array(dec10)
        colors = np.array(v10)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}10.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################

    ##########################################################################################
        # Full table
    
        ras = np.array(all_ras)
        decs = np.array(all_decs)
        colors = np.array(all_vhels)
    
        print
        print 'len(ras) : ',len(ras)
        print 'len(colors) : ',len(colors)
        print 'colors[:10] : ', colors[:10]
        print
    
        vmaxVal = max(colors)
        vminVal = min(colors)

        norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
        m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad, marker='o', c=colors, vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}all_sky.pdf'.format(saveDirectory),format='pdf')
         
        
##########################################################################################
##########################################################################################
##########################################################################################

    ra_shift = -180.0
#     xlabels = ['20h','22h','0h','2h','4h','6h','8h','10h','14h','16h','18h']

    xlabels = ['2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h']
#     xlabels = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']

    if plot_all_sky_12h:

    ##########################################################################################
        # 0 - 1000 km/s
        
        print 'ra1[0:5] : ',ra1[0:5]
        
        ras = np.array(ra1) + ra_shift
        decs = np.array(dec1)
        colors = np.array(v1)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

        print 'ras[0:5] : ',ras[0:5]

        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}1.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 1000 - 2000 km/s
    
        ras = np.array(ra2) + ra_shift
        decs = np.array(dec2)
        colors = np.array(v2)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=4, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}2.pdf'.format(saveDirectory),format='pdf')
    
    ##########################################################################################
        # 2000 - 3000 km/s
    
        ras = np.array(ra3) + ra_shift
        decs = np.array(dec3)
        colors = np.array(v3)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}3.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 3000 - 4000 km/s
    
        ras = np.array(ra4) + ra_shift
        decs = np.array(dec4)
        colors = np.array(v4)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}4.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 4000 - 5000 km/s
    
        ras = np.array(ra5) + ra_shift
        decs = np.array(dec5)
        colors = np.array(v5)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}5.pdf'.format(saveDirectory),format='pdf')



    ##########################################################################################
        # 5000 - 6000 km/s
    
        ras = np.array(ra6) + ra_shift
        decs = np.array(dec6)
        colors = np.array(v6)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}6.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 6000 - 7000 km/s
    
        ras = np.array(ra7) +ra_shift
        decs = np.array(dec7)
        colors = np.array(v7)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}7.pdf'.format(saveDirectory),format='pdf')
    
    
    
    ##########################################################################################
        # 7000 - 8000 km/s
    
        ras = np.array(ra8) +ra_shift
        decs = np.array(dec8)
        colors = np.array(v8)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}8.pdf'.format(saveDirectory),format='pdf')
    
    
    ##########################################################################################
        # 8000 - 9000 km/s
    
        ras = np.array(ra9) + ra_shift
        decs = np.array(dec9)
        colors = np.array(v9)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}9.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################
        # 9000 - 10000 km/s
    
        ras = np.array(ra10) + ra_shift
        decs = np.array(dec10)
        colors = np.array(v10)
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}10.pdf'.format(saveDirectory),format='pdf')

    ##########################################################################################

    ##########################################################################################
        # Full table
        import astropy.coordinates as coord
    
        ras = np.array(all_ras) + ra_shift
        decs = np.array(all_decs)
        colors = np.array(all_vhels)
        
#         colors = np.array([5000, 9000])
#         size_galaxy = 20
    
        print
        print 'len(ras) : ',len(ras)
        print 'len(colors) : ',len(colors)
        print 'colors[:10] : ', colors[:10]
        print
    
        vmaxVal = max(colors)
        vminVal = min(colors)

        norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
        m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
        # Transform the data into a SkyCoord object.
        c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


        # Matplotlib needs the coordinates in radians, so we have to convert them.
        # Furtermore matplotlib needs the RA coordinate in the range between -pi
        # and pi, not 0 and 2pi.
        ra_rad = c.ra.radian
        dec_rad = c.dec.radian
        ra_rad[ra_rad > np.pi] -= 2. * np.pi
        
#         ras = np.array([0, 10]) + 2.
#         decs = np.array([10, 20])
#         
#         ra = coord.Angle(ras*u.degree)
#         ra = ra.wrap_at(180*u.degree)
#         dec = coord.Angle(decs*u.degree)
        

        # Now plot the data in Aitoff projection with a grid.
        fig = plt.figure(figsize=(12,6))
        ax = lab.subplot(111,projection="aitoff")
        # lab.subplot(111)
        
        
#         plot1 = ax.scatter(ra.radian, dec.radian, marker='o', c=colors, vmin=vminVal,\
#         vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
        
        
        plot1 = plt.scatter(ra_rad, dec_rad, marker='o', c=colors, vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
        cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
        cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
        xlab = xlabels
        ax.set_xticklabels(xlab, weight=530)
        ax.grid(color='k', linestyle='solid', linewidth=0.5)
        tight_layout()
    
    #     plt.show()
        plt.savefig('{0}all_sky.pdf'.format(saveDirectory),format='pdf')

##########################################################################################
    
    # close the files
    inFile.close()
    
    print "Done."
    print
    
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
