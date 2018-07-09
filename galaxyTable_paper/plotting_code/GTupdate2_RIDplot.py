#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_RIDplot.py, v1.1 07/07/18

Update for thesis.


Makes a plot of RID*HubbleC - vcorr across the sky. (10/31/17)

MUST BE RUN INSIDE THE PYSALT CONDA ENVIRONMENT:
>>> source activate pysalt
>>> python GTupdate2_RIDplot.py
>>> ...sweet success...

'''

import sys
import os
import csv
# import string
# import warnings
# import numpy
# import atpy
import getpass
from utilities import *
from pylab import *

import math

# import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

# from astropy import wcs
# from astropy.io import fits
from scipy import stats

import numpy as np



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
        inFilename = '/Users/frenchd/Research/GT_update2_files/FinalGalaxyTable13_filtered.csv'
        saveDirectory = '/Users/frenchd/Research/GT_update2/'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
        
    all_sky_plot = False
    average_plot = True
        
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
    vhels = []
    vcorrs = []

    # these are for the outlier cases
    RID_vheldif_means = []
    RID_vheldif_medians = []
    ras_vhel = []
    decs_vhel = []
    
    
    cutoff = 3000
    
    # do the work
    count = 0
    for l in reader:
        count +=1
        
        Name = l['Name']
        NEDname = l['NEDname']
        vcorr = eval(l['vcorr'])
        vhel = eval(l['Vhel'])
        RID_mean = eval(l['RID_mean'])
        RID_median = eval(l['RID_median'])
        flag = eval(l['flag'])

        ra = eval(l['RAdeg'])
        dec = eval(l['DEdeg'])

        
        if RID_mean >= 0 and RID_median >= 0 and vhel >0:
            vhel = float(vhel)
            vcorr = float(vcorr)
            RID_mean = float(RID_mean)
            RID_median = float(RID_median)
        
            RID_dif_mean = vcorr/RID_mean
            RID_dif_median = vcorr/RID_median

            RID_dif_means.append(RID_dif_mean)
            RID_dif_medians.append(RID_dif_median)
            ras.append(ra)
            decs.append(dec)
            
            
            RID_vheldif_mean = vhel/RID_mean
            RID_vheldif_median = vhel/RID_median
            
            RID_vheldif_means.append(RID_vheldif_mean)
            RID_vheldif_medians.append(RID_vheldif_median)
            ras_vhel.append(ra)
            decs_vhel.append(dec)
            
            vhels.append(vhel)
            vcorrs.append(vcorr)
            
        
        # update the counter
        percentComplete = round((float(count)/130759)*100,2)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        sys.stdout.flush()
            
##########################################################################################
    if all_sky_plot:
        # plot it
        colmap = cm.RdBu_r
    
        vmaxVal = max(RID_dif_medians)
        vminVal = min(RID_dif_medians)
    
        ras = numpy.array(ras)
        decs = numpy.array(decs)
        colors = numpy.array(RID_dif_medians)
    
        ras_2 = numpy.array(ras_vhel)
        decs_2 = numpy.array(decs_vhel)
        colors_2 = numpy.array(RID_vheldif_medians)
    
    
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
        from mpl_toolkits.axes_grid1 import make_axes_locatable


        # Generate random data, for RA between 0 and 360 degrees, for DEC between
        # -90 and +90 degrees.
    #     ra_random = np.random.rand(100)*360.0
    #     dec_random = np.random.rand(100)*180.0-90.0
        ra_random = ras
        dec_random = decs
    
        ra_random2 = ras_vhel
        dec_random2 = decs_vhel

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
        plt.savefig('{0}RID_median-vcorr.pdf'.format(saveDirectory),format='pdf')


##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
    if average_plot:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        color_purple = '#7570b3'
        color_purple2 = '#984ea3'
        color_blue = '#436bad'      # french blue
        color_red = '#ec2d01'     # tomato red
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        

        alpha_points = 0.2
        alpha_hist = 0.9
        markerSize = 1
        
        binSize = 250
        bins = arange(0, 10250, binSize)
        
        label_points = r'$\rm Vhel/RID\_median$'
#         label_points = r'$\rm RID\_median/Vcorrs$'


#         color_isolated = color_green
#         color_assoc = color_orange
#         color_two = color_purple3
#         color_group = color_yellow

        color_points = 'grey'
        color_hist = 'black'
        color_hist_mean = color_blue

##########################################################################################
        # do the plotting 
             
        isolated_ys_mean = np.array(RID_vheldif_means)
        
        isolated_xs = vhels
        isolated_ys = np.array(RID_vheldif_medians)

#         isolated_xs = np.array(vcorrs)
#         isolated_ys = np.array(RID_dif_medians)

        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker='o',
                        c=color_points,
                        s=markerSize,
                        edgecolor=color_points,
                        alpha=alpha_points,
                        label=label_points)
                        
                        
        # median histogram
        bin_means, edges, binNumber = stats.binned_statistic(isolated_xs,
                                                            isolated_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='solid',
            color=color_hist,
            lw=2.0,
            alpha=alpha_hist,
            label=r'$\rm Median ~H_0$')
            
            
        # mean histograms
#         bin_means, edges, binNumber = stats.binned_statistic(isolated_xs,
#                                                             isolated_ys_mean,
#                                                             statistic='mean',
#                                                             bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,
#             Y,
#             ls='dashed',
#             color=color_hist_mean,
#             lw=2.0,
#             alpha=alpha_hist,
#             label=r'$\rm Mean ~H_0$')
            
            
            
        # associated
#         if plot_number >=2:
#             plot1 = scatter(associated_xs,
#                             associated_ys,
#                             marker=symbol_assoc,
#                             c=color_assoc,
#                             s=markerSize,
#                             edgecolor='black',
#                             alpha=alpha_assoc,
#                             label=label_assoc)
#         
#             # histogram associated
#             bin_means, edges, binNumber = stats.binned_statistic(associated_xs,
#                                                                 associated_ys,
#                                                                 statistic='mean',
#                                                                 bins=bins)
#             left,right = edges[:-1],edges[1:]        
#             X = array([left,right]).T.flatten()
#             Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#             plot(X,
#                 Y,
#                 ls='solid',
#                 color=color_assoc,
#                 lw=2.0,
#                 alpha=alpha_bins,
#                 label=r'$\rm Assoc. ~Mean ~EW$')
        
    
        
        # x-axis
        majorLocator   = MultipleLocator(2000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1000)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        xlabel(r'$\rm v_{{hel}}~[km s^{{-1}}]$')
#         xlabel(r'$\rm Vcorr~[km s^{{-1}}]$')
        ylabel(r'$H_{{\rm 0}}~\rm [km s^{{-1}}~Mpc^{{-1}}]$')
        leg = ax.legend(scatterpoints=5,prop={'size':12},loc='upper right',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        # grid
        ax.grid(b=None,which='major',axis='both')
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed')
        
        # axis limits
        ylim(30,100)
        xlim(0, 10000)

        if average_plot:
            savefig('{0}/H0(vhels)_median.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()


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
