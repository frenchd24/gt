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
# import os
import csv
import numpy as np
import getpass
from utilities import *
# from pylab import *

import math
import matplotlib.pyplot as plt


from matplotlib import rc
fontScale = 18
rc('text', usetex=True)
rc('font', size=18, family='serif', weight='normal')
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
    
    def magvdiam(d,a,b):
        '''
            functional form (from Wakker & Savage 2009):
                
                M = -5.78 log R - 12
        '''
        try:
            M = a * np.log10(d) + b
        except Exception,e:
            print 'failed on d = {0}'.format(d)
            M = 0
            
        return M
        
    
    # plot it
    fig = plt.figure(figsize=(12,6))
    ax = plt.subplot(111)
    
    plot1 = plt.scatter(MajDiamL,BmagL,marker='o',lw=0,s=4,alpha=0.8,c='grey')
    plt.xlabel(r'$\rm Diameter ~ [kpc]$')
    plt.ylabel(r'$\rm M_B$')
    plt.ylim(-6,-26)
    ax.set_xscale('log')
    plt.xlim(1,90)
    
    
    from scipy.interpolate import interp1d
    from scipy import interpolate, optimize
    
    a = -5.8
    b = -12.
#     y = magvdiam(MajDiamL, a, b)
    
    # form is optimize.curve_fit(function_name, x_data, y_data)
    
    xFit = []
    yFit = []
    for m,b in zip(MajDiamL,BmagL):
        if m >0:
            xFit.append(m)
            yFit.append(b)
    
    xFit = np.array(xFit)
    yFit = np.array(yFit)
    
    popt, pcov = optimize.curve_fit(magvdiam, xFit, yFit)
    perr = np.sqrt(np.diag(pcov))
    print
    print 'popt: ',popt
    print
    print 'pcov: ',pcov
    print
    print 'perr: ',perr
    print
    
    
    a = round(popt[0],3)
    b = round(popt[1],3)
    
    x2 = np.linspace(0.1,200,500)
    
    plt.plot(x2, magvdiam(x2, *popt), 'r-',color='blue',alpha = 0.8)
# ,label=r'$\rm M = {0}, b={1}'.format(a,b)
    plt.legend()

    # x-axis
#     ax.xaxis.set_major_locator(plt.LogLocator(base=10.0, numticks=20))
#     majors = [0, 1, 2, 5, 10, 20, 50, 100]
#     ax.xaxis.set_major_locator(plt.FixedLocator(majors))
#     minors = np.append(np.linspace(0, 1, 11)[1:-1], np.array([30, 40, 60, 70, 80, 90]))
#     ax.xaxis.set_minor_locator(plt.FixedLocator(minors))


    majorLocator   = plt.MultipleLocator(10)
    majorFormatter = plt.FormatStrFormatter(r'$\rm %d$')
#     minorLocator   = plt.MultipleLocator(10)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(majorFormatter)
#     ax.xaxis.set_minor_locator(minorLocator)
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100])


    # y axis
    majorLocator   = plt.MultipleLocator(2)
    majorFormatter = plt.FormatStrFormatter(r'$\rm %d$')
    minorLocator   = plt.MultipleLocator(0.5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_minor_locator(minorLocator)


    # make the legend
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    blue_line = mlines.Line2D([], [], color='blue', alpha = 0.8,
    markersize=10, label=r'$\rm M = {0} ~ Log_{{10}}(D) {1}$'.format(a,b))
                              
    plt.legend(handles=[blue_line],loc='upper left')


    plt.tight_layout()
    
    saveDirectory = '/Users/frenchd/Research/GT_update2/galaxyTable_paper/figures/'
    plt.savefig('{0}mag_v_diam_fit.pdf'.format(saveDirectory),format='pdf')


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
