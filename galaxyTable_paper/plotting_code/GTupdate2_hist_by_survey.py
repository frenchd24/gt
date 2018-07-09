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


Make a violin plot

'''

import sys
# import os
import csv
import numpy as np
import getpass
from utilities import *
# from pylab import *

import math
import matplotlib.pyplot as plt


from matplotlib import rc
fontScale = 16
rc('text', usetex=True)
rc('font', size=16, family='serif', weight='normal')
rc('xtick.major',size=8,width=0.6)
rc('xtick.minor',size=5,width=0.6)
rc('ytick.major',size=8,width=0.6)
rc('ytick.minor',size=5,width=0.6)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
rc('axes',labelsize = fontScale)
rc('xtick', labelsize = fontScale)
rc('ytick',labelsize = fontScale)
# rc('font', weight = 450)
# rc('axes',labelweight = 'bold')
rc('axes',linewidth = 1)
rc('axes',labelsize=fontScale+4)


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



# def median_low(l):
#     # returns the closest element in a list to the median, rounding down
# 
#     # E.g.,
#     # list = [1, 2, 3, 4, 5]
#     # median_low(list) = 3
#     #
#     # list = [1, 2, 3, 4]
#     # median_low(list) = 2
#     
# #     l.sort()
#     l = np.array(l)
#     med = np.median(l)
#     
#     diff = abs(np.subtract(med,l))
#     
#     diff = list(diff)
#     l = list(l)
#     indexMin = diff.index(min(diff))
#     
#     return l[indexMin]

    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':        
        inFilename = '/Users/frenchd/Research/gt/FinalGalaxyTable13_filtered.csv'
        namesFilename = '/Users/frenchd/Research/gt/FinalGalaxyTable13_altNames.csv'

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
    namesFile = open(namesFilename,'rU')
    names_reader = csv.DictReader(namesFile)
    
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
    altNames = []
    
    sdss_vhel = []
    otherThing = []
    two_mass_vhel = []
    two_df_vhel = []
    six_df_vhel = []
    rc3_vhel = []
    other_vhel = []
        
    # do the work
    count = 0
    for l,n in zip(reader,names_reader):
        count +=1
        
        altNames = str(n['altNames'])
        
        Name = l['Name']
        NEDname = l['NEDname']
        vcorr = eval(l['vcorr'])
        RID_mean = eval(l['RID_mean'])
        RID_median = eval(l['RID_median'])
        bestDist = l['bestDist']
        ebminusv = l['E(B-V)']
        vhel = eval(l['Vhel'])
        flag = eval(l['flag'])
        Bmag = l['Bmag']
        MajDiam = l['MajDiam']
        RC3_type = l['RC3_type']

        ra = eval(l['RAdeg'])
        dec = eval(l['DEdeg'])
        
        found = False
#         if flag != 1:
        if flag == 0:

            if bfind(altNames,'2MAS'):
                two_mass_vhel.append(vhel)
                found = True
            
            if bfind(altNames,'SDSS'):
                sdss_vhel.append(vhel)
                found = True
            
            if bfind(altNames,'2dF'):
                two_df_vhel.append(vhel)
                found = True

            if bfind(altNames,'6dF'):
                six_df_vhel.append(vhel)
                found = True

            if RC3_type !='x':
                rc3_vhel.append(vhel)
                found = True

            if not found:
                other_vhel.append(vhel)
        

        # update the counter
        percentComplete = round((float(count)/130759)*100,2)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        sys.stdout.flush()
            
##########################################################################################

    
    # plot it
    fig = plt.figure(figsize=(12,8))
    ax = plt.subplot(111)
    
    print 'len(two_mass_vhel): ',len(two_mass_vhel)
    print
    print 'len(sdss_vhel): ',len(sdss_vhel)
    print
    print 'len(two_df_vhel): ',len(two_df_vhel)
    print
    print 'len(six_df_vhel): ',len(six_df_vhel)
    print
    print 'len(rc3_vhel): ',len(rc3_vhel)
    print
    print 'len(other_vhel): ',len(other_vhel)
    print
    
    
#     prop_cycle = plt.rcParams['axes.prop_cycle']
#     colors = prop_cycle.by_key()['color']
#     
#     for i in colors:
#         print 'i: ',i

    
    n_bins = 35
    stack = False
#     plot1 = plt.scatter(MajDiamL,BmagL,marker='o',lw=0,s=4,alpha=0.8,c='grey')

    color_blue = '#436bad'      # french blue
    color_red = '#ec2d01'     # tomato red

    color_green = '#1b9e77'
    color_purple = '#7570b3'
    color_orange = '#d95f02'
    color_purple2 = '#984ea3'
    
    color_green2 = '#4daf4a'
    color_orange = '#d95f02'
    color_purple3 = '#7570b3'
    color_pink = '#e7298a'
    color_lime = '#66a61e'
    color_yellow = '#e6ab02'
    color_brown = '#a6761d'
    color_coal = '#666666'
    
    color_blue2 = '#377eb8'


    style_2mass = 'solid'
    marker_2mass = ''
    color_2mass = 'black'
    lw_2mass = 2.5
    
    style_sdss = '--'
    marker_sdss = ''
#     color_sdss = 'red'
    color_sdss = color_red
    lw_sdss = 2
    
    style_2df = 'solid'
    marker_2df = 'x'
#     color_2df = 'goldenrod'
    color_2df = color_yellow
    lw_2df = 1.5
    
    style_6df = '-.'
    marker_6df = ''
#     color_6df = 'darkgreen'
#     color_6df = color_purple2
    color_6df = 'darkgreen'
    lw_6df = 2
    
    style_rc3 = ':'
    marker_rc3 = 'd'
#     color_rc3 = 'blue'
    color_rc3 = color_blue2
    lw_rc3 = 2
    
    style_other = 'solid'
    marker_other = ''
#     color_other = 'grey'
    color_other = 'grey'
    lw_other = 1.5

##########################################################################################
    # 2MAS
#     ax.hist(two_mass_vhel, n_bins, histtype='step', stacked=stack, fill=False,
#     label = r'$\rm 2MASS$',color='black',lw=3,ls='solid')
    
    hist, bin_edges = np.histogram(two_mass_vhel, bins=n_bins)
#     print 'hist: ',hist
#     print 'bin_edges: ',bin_edges
#     print

    ax.plot(bin_edges[:-1], hist, color=color_2mass, lw=lw_2mass, label= r'$\rm 2MASS$', ls=style_2mass, marker=marker_2mass)
##########################################################################################
    # SDSS
#     ax.hist(sdss_vhel, n_bins, histtype='step', stacked=stack, fill=False,
#     label = r'$\rm SDSS$',color='red',lw=2,ls='dashed')
    
    hist, bin_edges = np.histogram(sdss_vhel, bins=n_bins)
    ax.plot(bin_edges[:-1], hist, color=color_sdss, lw=lw_sdss, label= r'$\rm SDSS$', ls=style_sdss, marker=marker_sdss)
##########################################################################################
    # 2df
#     ax.hist(two_df_vhel, n_bins, histtype='step', stacked=stack, fill=False,
#     label = r'$\rm 2dF$',color='goldenrod',lw=2,ls='dotted')
    
    hist, bin_edges = np.histogram(two_df_vhel, bins=n_bins)
    ax.plot(bin_edges[:-1], hist, color=color_2df, lw=lw_2df, label= r'$\rm 2dF$', ls=style_2df, marker=marker_2df,ms=10)
##########################################################################################
    # 6df
#     ax.hist(six_df_vhel, n_bins, histtype='step', stacked=stack, fill=False,
#     label = r'$\rm 6dF$',color='darkgreen',lw=2,marker='-.')
    
    hist, bin_edges = np.histogram(six_df_vhel, bins=n_bins)
    ax.plot(bin_edges[:-1], hist, color=color_6df, lw=lw_6df, label= r'$\rm 6dF$', ls=style_6df, marker=marker_6df)
##########################################################################################
    # RC3
#     ax.hist(rc3_vhel, n_bins, histtype='step', stacked=stack, fill=False,
#     label = r'$\rm RC3$',lw=2,color='blue',ls='dotted')
    
    hist, bin_edges = np.histogram(rc3_vhel, bins=n_bins)
    ax.plot(bin_edges[:-1], hist, color=color_rc3, lw=lw_rc3, label= r'$\rm RC3$', ls=style_rc3, marker=marker_rc3)
##########################################################################################
    # all else
#     ax.hist(other_vhel, n_bins, histtype='step', stacked=stack, fill=False,
#     label = r'$\rm Other$',color='grey',lw=1,ls='solid')
    
    hist, bin_edges = np.histogram(other_vhel, bins=n_bins)
    ax.plot(bin_edges[:-1], hist, color=color_other, lw=lw_other, label= r'$\rm Other$', ls=style_other, marker=marker_other)
##########################################################################################

    plt.xlabel(r'$\rm V_{hel}~ [km s^{{-1}}]$')
    plt.ylabel(r'$\rm Number $')

#     plt.ylim(-6,-26)
#     ax.set_xscale('log')
#     plt.xlim(1,90)

    # x-axis
    majorLocator   = plt.MultipleLocator(2000)
    majorFormatter = plt.FormatStrFormatter(r'$\rm %d$')
    minorLocator   = plt.MultipleLocator(500)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.xaxis.set_minor_locator(minorLocator)
#     ax.set_xticks([1, 2, 5, 10, 20, 50, 100])


    # y axis
    majorLocator   = plt.MultipleLocator(1000)
    majorFormatter = plt.FormatStrFormatter(r'$\rm %d$')
    minorLocator   = plt.MultipleLocator(200)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_minor_locator(minorLocator)


    # make the legend
#     import matplotlib.patches as mpatches
#     import matplotlib.lines as mlines
# 
#     blue_line = mlines.Line2D([], [], color='blue', alpha = 0.8,
#     markersize=10, label=r'$\rm M = {0} ~ Log_{{10}}(D) {1}$'.format(a,b))
#                               
#     plt.legend(handles=[blue_line],loc='upper left')
    
    plt.xlim(0,10000)
    plt.legend(loc='upper left')
    plt.tight_layout()
    
    saveDirectory = '/Users/frenchd/Research/GT_update2/galaxyTable_paper/figures/'
    plt.savefig('{0}hist_by_survey8_filtered_0only.pdf'.format(saveDirectory),format='pdf')


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
