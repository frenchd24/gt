 #!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: make_Lstar_histograms3_vert.py v 4 03/09/18

Comes from: filament_make_histograms.py, v 1.0 08/11/2015

but updated for the entire data set for the pilot paper. Makes photometry plots for the 
full galaxy table. (11/17/2015)

v1.1 Make nicer looking Texify'd plots (4/27/16)

v2: change name from 'pilot_make_histograms.py' to 'pilot_make_Lstar_histograms2.py'
    Make a panel of plots for different v_max values
    (6/3/16)

v2.1: formatting updates for the paper (8/11/16)

v3: remove the functions that are not used anymore. This ONLY plots the single,
    multi-panel plot used in the pilot paper (9/8/16)
    
v3.1: slight update to remove space between the plots (i.e. share axes)
    -> (09/27/16)
    
v3.2: slight update for the first referee report (12/12/16)
    - make labels bigger, try a vertical stack instead of 2x2 panels
    
v4: update for galaxy table paper (03/09/18)

'''

import sys
# import os
import csv
from pylab import *
# import atpy
import math
import getpass
import scipy.optimize as optimization
import pickle
import itertools
from utilities import *

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


def schechter(m,phi,mstar,alpha):
    # construct a Schechter luminosity function and return the associated density function
    
#     s = 0.4*log(10)*phi*(10**(0.4*(mstar-m)))**(alpha +1) * exp(-10**(0.4*(mstar-m)))
#     s = 0.4*log(10)*phi*(10**(0.4*(m-mstar)*(alpha +1))) * exp(-10**(0.4*(m-mstar)))

    s = 0.4 * math.log10(10) * phi * (10**(0.4*(mstar-m)))**(alpha +1) * exp(-10**(0.4*(mstar-m)))

    return s
    
    
def apparentMag(M,d):
    # convert M to apparent magnitude at distance d (input in Mpc)
    
    return 5*log10(d*10**6)-5+ M
    
    
    
def absoluteMag_noExtinc(m,d):
    # m is apparent magnitude, d is distance in Mpc
    M = float(m) - 5*math.log10((float(d)*10**6)/10)
    return M


def absoluteMag(m,d,e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)
    return M
    

# def extinction(m,e):
#     # m is apparent magnitude, e is extinction E(B-V)
#     
#     M = float(m) - 3.1*float(e)
#     return M
    

def lstarValue(mstar,m):
    # calculate and return L/Lstar
    lratio = 10**(-0.4*(m-mstar))
    return lratio



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


def make_histogram_lstar_split(d1,vlo,vhi,saveDirectory,save):
    # unpack data
    # d_gt = {'Bmedian':gt_BmedianList,'BLstar':gt_BLstarList,'ra':gt_raList,'dec':gt_decList,'dist':gt_distList,'vcorr':gt_vcorrList}
    #
    # make a 4 panel Lstar histogram, each panel with progressively higher v_max values
    #
    # make a vertical stack instead of 2x2
    
#     phot = d1['Bmedian']
#     BLstar = d1['BLstar']
#     ra = d1['ra']
#     dec = d1['dec']
#     dist = d1['dist']
#     vcorr = d1['vcorr']

    phot = d1['Bmag']
#     BLstar = d1['Lstar']
    BLstar = d1['Lstar']

    ra = d1['RAdeg']
    dec = d1['DEdeg']
    dist = d1['bestDist']
#     vcorr = d1['vcorr']
    vhel = d1['Vhel']
    flag = d1['flag']
    name = d1['Name']
    
    MajDiam = d1['majorAxis']
    bestDist = d1['bestDist']
    
#     BLstar = []
#     for m, d, n in zip(MajDiam, bestDist, name):
#         if m >0 and d>0:
#             M = -4.696 * math.log10(float(m)) - 12.963
#            
#             mstar = -19.57
#             lstar = lstarValue(mstar, M)
#             BLstar.append(lstar)
#            
#             if n == 'IC0210':
#                 print '---- IC0210 ----'
#                 print 'MajDiam = ',m
#                 print 'M = ',M
#                 print 'lstar = ',lstar
#                 print
#                 print
#         else:
#             BLstar.append(0)

    
    # breakdown of velocity splits
    vlo_1 = 0
    vlo_2 = 2500
    vlo_3 = 6000
    vlo_4 = 8000

    vhi_1 = 2500
    vhi_2 = 6000
    vhi_3 = 8000
    vhi_4 = 10000
    
    label1 = r'$\rm {0} \leq cz < {1}~ km\, s^{{-1}}$'.format(vlo_1,vhi_1)
    label2 = r'$\rm {0} \leq cz < {1}~ km\, s^{{-1}}$'.format(vlo_2,vhi_2)
    label3 = r'$\rm {0} \leq cz < {1}~ km\, s^{{-1}}$'.format(vlo_3,vhi_3)
    label4 = r'$\rm {0} \leq cz \leq {1}~ km\, s^{{-1}}$'.format(vlo_4,vhi_4)
    
    binsize = 0.2
    
#     fig = figure(figsize=(14,9))
    fig = figure(figsize=(10,14))
    subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.001,hspace=0.001)

    # first
    ax = fig.add_subplot(4,1,1)
#     ax = gca()
    
    cutPhot = []
    cutBLstar = []
    cutRA = []
    cutDec = []
    cutDist = []
    cutVhel = []
    
    for p, b, r, d, dis, v, f, n in zip(phot, BLstar, ra, dec, dist, vhel, flag, name):
        if float(v) >= vlo_1 and float(v) < vhi_1:
            if p > 0 and b !=0:
                logLstar = log10(float(b))
                if not isNumber(logLstar):
                    print 'logLstar, Lstar: ',logLstar,', ',b
                    print 'phot: ',p
                cutPhot.append(float(p))
                cutBLstar.append(logLstar)
                cutRA.append(float(r))
                cutDec.append(float(d))
                cutDist.append(float(dis))
                cutVhel.append(float(v))
                
                if logLstar < -3.9:
                    print '{0} - Log(Lstar) = {1}, Lstar = {2}, B = {3}, {4}'.format(n, logLstar, b, p, v)
            
#     print 'cutBLstar: ',cutBLstar
    print
    print 'max: ',max(cutBLstar)
    print 'min: ',min(cutBLstar)
    print 'med: ',median(cutBLstar)
    print 'len: ',len(cutBLstar)
    print

#     setbins = arange(-3,2,binsize)
    setbins = arange(-4,2,binsize)
    counts,bins = histogram(cutBLstar,bins=setbins)
    
    ax.hist(cutBLstar,setbins,histtype='bar',color='grey',label=label1)
    maxHeight = max(counts)+(max(counts)*0.25)
#     ylim(0,maxHeight)
    
#     ylabel(r'$\rm Number$')
#     xlabel(r'$\rm log_{10} (L/L_*)$')
#     ax.annotate(label1,xy=(log10(1.2),maxHeight*0.7))
#     title(label1)

#     these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', alpha=1, facecolor='none')

#     place a text box in upper right in axes coords
#     ax.text(0.702, 0.95, label1, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax.text(0.02, 0.95, label1, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    
    # x coordinate adjustment for annotations
    xAn = -0.2
    
    # label color
    lcolor = 'red'
    lstyle = 'dashed'
    lalpha = 0.7

    
    # label size
    lsize = 15
    
    # draw and annotate Lstar = 1 line
    axvline(x=log10(1),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
    l1 = r'$\rm L^*=1.0$'
    annotate(s=l1,xy=(log10(1)+xAn,maxHeight*0.96),size=lsize)
    
    # draw Lstar = 0.5 line
    axvline(x=log10(0.5),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
    l_5 = r'$\rm L^*=0.5$'
    annotate(s=l_5,xy=(log10(0.5)+xAn+xAn/1.3,maxHeight*0.85),size=lsize)
  
    # draw Lstar = 0.1 line
    axvline(x=log10(0.1),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
    l_1 = r'$\rm L^*=0.1$'
    annotate(s=l_1,xy=(log10(0.1)+xAn,maxHeight*0.96),size=lsize)
    
    # draw Lstar = 0.05 line
    axvline(x=log10(0.05),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
    l_05 = r'$\rm L^*=0.05$'
    annotate(s=l_05,xy=(log10(0.05)+xAn+xAn/1.3,maxHeight*0.85),size=lsize)
    
    # draw Lstar = 0.01 line
    axvline(x=log10(0.01),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
    l_01 = r'$\rm L^*=0.01$'
    annotate(s=l_01,xy=(log10(0.01)+xAn,maxHeight*0.96),size=lsize)
    
    # format the axes
    #
    # x-axis
    majorLocator   = MultipleLocator(1)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(0.5)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.set_xticks([-4,-3,-2,-1,0,1,2])
    plt.setp(ax.get_xticklabels(), visible=False)


    # y axis
    majorLocator   = MultipleLocator(200)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(100)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_minor_locator(minorLocator)
    
    ax.set_ylim([0,1000])
    ax.set_xlim([-4,2])

#     ax.yaxis.tick_right()

#     tight_layout()
    
##########################################################################################

    # second
#     ax2 = fig.add_subplot(4,1,2)
    ax2 = fig.add_subplot(4,1,2,sharex=ax)

#     ax = gca()
    
    cutPhot2 = []
    cutBLstar2 = []
    cutRA2 = []
    cutDec2 = []
    cutDist2 = []
    cutVhel2 = []
    
    for p, b, r, d, dis, v, f in zip(phot, BLstar, ra, dec, dist, vhel, flag):
        if float(v) >= vlo_2 and float(v) < vhi_2:
#             print 'v: ',v
            if p > 0 and b !=0:
                cutPhot2.append(float(p))
                cutBLstar2.append(log10(float(b)))
                cutRA2.append(float(r))
                cutDec2.append(float(d))
                cutDist2.append(float(dis))
                cutVhel2.append(float(v))
            
#     print 'cutBLstar: ',cutBLstar2
    print
    print 'max: ',max(cutBLstar2)
    print 'min: ',min(cutBLstar2)
    print 'med: ',median(cutBLstar2)
    print 'len2: ',len(cutBLstar2)
    print


#     setbins = arange(-3,2,binsize)
    setbins = arange(-4,2,binsize)
    counts2,bins2 = histogram(cutBLstar2,bins=setbins)
    
    ax2.hist(cutBLstar2,setbins,histtype='bar',color='grey',label=label2)
    maxHeight = max(counts2)+(max(counts2)*0.25)
#     ylim(0,maxHeight)
    
#     ylabel(r'$\rm Number$')
#     xlabel(r'$\rm log_{10} (L/L_*)$')
#     ax.annotate(label2,xy=(log10(1.2),maxHeight*0.7))
#     title(label2)

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', alpha=1, facecolor='none')

    # place a text box in upper right in axes coords
#     ax2.text(0.662, 0.95, label2, transform=ax2.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax2.text(0.02, 0.95, label2, transform=ax2.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    # x coordinate adjustment for annotations
    xAn = -0.2
    
    # label color
    lcolor = 'red'
    lstyle = 'dashed'
    
    # label size
    lsize = 14
    
#     draw and annotate Lstar = 1 line
    axvline(x=log10(1),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l1 = r'$\rm L_*=1.0$'
#     annotate(s=l1,xy=(log10(1)+xAn,maxHeight*0.92),size=lsize)
    
#     draw Lstar = 0.5 line
    axvline(x=log10(0.5),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_5 = r'$\rm L_*=0.5$'
#     annotate(s=l_5,xy=(log10(0.5)+xAn+xAn,maxHeight*0.85),size=lsize)
  
#     draw Lstar = 0.1 line
    axvline(x=log10(0.1),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_1 = r'$\rm L_*=0.1$'
#     annotate(s=l_1,xy=(log10(0.1)+xAn,maxHeight*0.92),size=lsize)
    
#     draw Lstar = 0.05 line
    axvline(x=log10(0.05),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_05 = r'$\rm L_*=0.05$'
#     annotate(s=l_05,xy=(log10(0.05)+xAn+xAn,maxHeight*0.85),size=lsize)
    
#     draw Lstar = 0.01 line
    axvline(x=log10(0.01),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_01 = r'$\rm L_*=0.01$'
#     annotate(s=l_01,xy=(log10(0.01)+xAn,maxHeight*0.92),size=lsize)
    
    # format the axes
    #
    # x-axis
    majorLocator   = MultipleLocator(1)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(0.5)
    ax2.xaxis.set_major_locator(majorLocator)
    ax2.xaxis.set_major_formatter(majorFormatter)
    ax2.xaxis.set_minor_locator(minorLocator)
    ax2.set_xticks([-4,-3,-2,-1,0,1,2])
    plt.setp(ax2.get_xticklabels(), visible=False)


    # y axis
    majorLocator   = MultipleLocator(1000)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(500)
    ax2.yaxis.set_major_locator(majorLocator)
    ax2.yaxis.set_major_formatter(majorFormatter)
    ax2.yaxis.set_minor_locator(minorLocator)
#     plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=True)
    ax2.set_yticks([0,1000,2000])


    ax2.set_ylim([0,3000])
    ax2.set_xlim([-4,2])

#     ax2.yaxis.tick_right()
#     tight_layout()

##########################################################################################

    # third
    ax3 = fig.add_subplot(4,1,3,sharex=ax)
#     ax = gca()
    
    cutPhot3 = []
    cutBLstar3 = []
    cutRA3 = []
    cutDec3 = []
    cutDist3 = []
    cutVhel3 = []
    
    for p, b, r, d, dis, v, f in zip(phot, BLstar, ra, dec, dist, vhel, flag):
        if float(v) >= vlo_3 and float(v) < vhi_3:
#             print 'v: ',v
            if p > 0 and b !=0:
                cutPhot3.append(float(p))
                cutBLstar3.append(log10(float(b)))
                cutRA3.append(float(r))
                cutDec3.append(float(d))
                cutDist3.append(float(dis))
                cutVhel3.append(float(v))
            
#     print 'cutBLstar: ',cutBLstar3
    print
    print 'max: ',max(cutBLstar3)
    print 'min: ',min(cutBLstar3)
    print 'med: ',median(cutBLstar3)
    print 'len3: ',len(cutBLstar3)
    print


#     setbins = arange(-3,2,binsize)
    setbins = arange(-4,2,binsize)
    counts3,bins3 = histogram(cutBLstar3,bins=setbins)
    
    ax3.hist(cutBLstar3,setbins,histtype='bar',color='grey',label=label3)
    maxHeight = max(counts3)+(max(counts3)*0.25)
#     ylim(0,maxHeight)
    
#     ylabel(r'$\rm Number$')
    xlabel(r'$\rm log_{10} (L/L^*)$')
#     ax.annotate(label3,xy=(log10(1.2),maxHeight*0.7))
#     title(label3)

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', alpha=1, facecolor='none')

    # place a text box in upper right in axes coords
#     ax3.text(0.663, 0.95, label3, transform=ax3.transAxes, fontsize=14,verticalalignment='top', bbox=props)
    ax3.text(0.02, 0.95, label3, transform=ax3.transAxes, fontsize=14,verticalalignment='top', bbox=props)

    # x coordinate adjustment for annotations
    xAn = -0.2
    
    # label color
    lcolor = 'red'
    lstyle = 'dashed'
    
    # label size
    lsize = 14
    
#     draw and annotate Lstar = 1 line
    axvline(x=log10(1),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l1 = r'$\rm L_*=1.0$'
#     annotate(s=l1,xy=(log10(1)+xAn,maxHeight*0.92),size=lsize)
    
#     draw Lstar = 0.5 line
    axvline(x=log10(0.5),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_5 = r'$\rm L_*=0.5$'
#     annotate(s=l_5,xy=(log10(0.5)+xAn+xAn,maxHeight*0.85),size=lsize)
  
#     draw Lstar = 0.1 line
    axvline(x=log10(0.1),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_1 = r'$\rm L_*=0.1$'
#     annotate(s=l_1,xy=(log10(0.1)+xAn,maxHeight*0.92),size=lsize)
    
#     draw Lstar = 0.05 line
    axvline(x=log10(0.05),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_05 = r'$\rm L_*=0.05$'
#     annotate(s=l_05,xy=(log10(0.05)+xAn+xAn,maxHeight*0.85),size=lsize)
    
#     draw Lstar = 0.01 line
    axvline(x=log10(0.01),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_01 = r'$\rm L_*=0.01$'
#     annotate(s=l_01,xy=(log10(0.01)+xAn,maxHeight*0.92),size=lsize)
    
    # format the axes
    #
    # x-axis
    majorLocator   = MultipleLocator(1)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(0.5)
    ax3.xaxis.set_major_locator(majorLocator)
    ax3.xaxis.set_major_formatter(majorFormatter)
    ax3.xaxis.set_minor_locator(minorLocator)
    ax3.set_xticks([-4,-3,-2,-1,0,1,2])
    plt.setp(ax3.get_xticklabels(), visible=False)

    # y axis
    majorLocator   = MultipleLocator(1000)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(500)
    ax3.yaxis.set_major_locator(majorLocator)
    ax3.yaxis.set_major_formatter(majorFormatter)
    ax3.yaxis.set_minor_locator(minorLocator)
    ax3.set_yticks([0,1000,2000,3000])

    ax3.set_ylim([0,4000])
    ax3.set_xlim([-4,2])

#     ax3.yaxis.tick_right()

#     tight_layout()

##########################################################################################

    # fourth
#     ax4 = fig.add_subplot(4,1,4,sharex=ax2)
    ax4 = fig.add_subplot(4,1,4,sharex=ax)

#     ax = gca()
    
    cutPhot4 = []
    cutBLstar4 = []
    cutRA4 = []
    cutDec4 = []
    cutDist4 = []
    cutVhel4 = []
    
    for p, b, r, d, dis, v, f in zip(phot, BLstar, ra, dec, dist, vhel, flag):
        if float(v) >= vlo_4 and float(v) <= vhi_4:
#             print 'v: ',v
            if p > 0 and b !=0:
                cutPhot4.append(float(p))
                cutBLstar4.append(log10(float(b)))
                cutRA4.append(float(r))
                cutDec4.append(float(d))
                cutDist4.append(float(dis))
                cutVhel4.append(float(v))
            
#     print 'cutBLstar: ',cutBLstar4
    print
    print 'max: ',max(cutBLstar4)
    print 'min: ',min(cutBLstar4)
    print 'med: ',median(cutBLstar4)
    print 'len4: ',len(cutBLstar4)
    print


#     setbins = arange(-3,2,binsize)
    setbins = arange(-4,2,binsize)
    counts4,bins4 = histogram(cutBLstar4,bins=setbins)
        
    ax4.hist(cutBLstar4,setbins,histtype='bar',color='grey',label=label4)
    maxHeight = max(counts4)+(max(counts4)*0.25)
#     ylim(0,maxHeight)
    
#     ylabel(r'$\rm Number$')
    xlabel(r'$\rm log_{10} (L/L^*)$')
#     ax.annotate(label4,xy=(log10(1.2),maxHeight*0.7))
#     title(label4)

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', alpha=1, facecolor='none')

    # place a text box in upper right in axes coords
#     ax4.text(0.650, 0.95, label4, transform=ax4.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax4.text(0.02, 0.95, label4, transform=ax4.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    # x coordinate adjustment for annotations
    xAn = -0.2
    
    # label color
    lcolor = 'red'
    lstyle = 'dashed'
    
    # label size
    lsize = 14
    
#     draw and annotate Lstar = 1 line
    axvline(x=log10(1),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l1 = r'$\rm L_*=1.0$'
#     annotate(s=l1,xy=(log10(1)+xAn,maxHeight*0.92),size=lsize)
#     
#     draw Lstar = 0.5 line
    axvline(x=log10(0.5),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_5 = r'$\rm L_*=0.5$'
#     annotate(s=l_5,xy=(log10(0.5)+xAn+xAn,maxHeight*0.85),size=lsize)
#   
#     draw Lstar = 0.1 line
    axvline(x=log10(0.1),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_1 = r'$\rm L_*=0.1$'
#     annotate(s=l_1,xy=(log10(0.1)+xAn,maxHeight*0.92),size=lsize)
#     
#     draw Lstar = 0.05 line
    axvline(x=log10(0.05),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_05 = r'$\rm L_*=0.05$'
#     annotate(s=l_05,xy=(log10(0.05)+xAn+xAn,maxHeight*0.85),size=lsize)
#     
#     draw Lstar = 0.01 line
    axvline(x=log10(0.01),linewidth=1, color=lcolor,linestyle=lstyle,alpha=lalpha)
#     l_01 = r'$\rm L_*=0.01$'
#     annotate(s=l_01,xy=(log10(0.01)+xAn,maxHeight*0.92),size=lsize)
    
    # format the axes
    #
    # x-axis
    majorLocator   = MultipleLocator(1)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(0.5)
    ax4.xaxis.set_major_locator(majorLocator)
    ax4.xaxis.set_major_formatter(majorFormatter)
    ax4.xaxis.set_minor_locator(minorLocator)
    ax4.set_xticks([-4,-3,-2,-1,0,1,2])

    # y axis
    majorLocator   = MultipleLocator(1000)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(500)
    ax4.yaxis.set_major_locator(majorLocator)
    ax4.yaxis.set_major_formatter(majorFormatter)
    ax4.yaxis.set_minor_locator(minorLocator)
    ax4.set_yticks([0,1000,2000,3000,4000])


#     ax4.set_ylim([0,4050])
    ax4.set_ylim([0,5000])
    ax4.set_xlim([-4,2])

#     ax4.yaxis.tick_right()


    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.ylabel(r'$\rm Number$',labelpad=20)
#     plt.tight_layout()


    if save:
        savefig('{0}Lstar_histogram_4bins_final_{1}-{2}_v3_vert_flag0.pdf'.format(saveDirectory, vlo_1, vhi_4),format='pdf',bbox_inches='tight')
    else:
        show()
       
        
        
def main():
        
    '''
    The data set looks like this:
    
    d_b_j = {'phot':b_j_phot,'err_up':b_j_up_phot,'err_lo':b_j_lo_phot,'ra':b_j_raList,'dec':b_j_decList,'dist':b_j_distList,'vcorr':b_j_vcorrList}
    d_b_all = {'phot':b_all_phot,'err_up':b_all_up_phot,'err_lo':b_j_lo_phot,'ra':b_all_raList,'dec':b_all_decList,'dist':b_all_distList,'vcorr':b_all_vcorrList}
    d_j = {'phot':j_phot,'err_up':j_up_phot,'err_lo':j_lo_phot,'ra':j_raList,'dec':j_decList,'dist':j_distList,'vcorr':j_vcorrList}
    d_k = {'phot':k_phot,'err_up':k_up_phot,'err_lo':k_lo_phot,'ra':k_raList,'dec':k_decList,'dist':k_distList,'vcorr':k_vcorrList}
    d_h = {'phot':h_phot,'err_up':h_up_phot,'err_lo':h_lo_phot,'ra':h_raList,'dec':h_decList,'dist':h_distList,'vcorr':h_vcorrList}
    d_g = {'phot':g_phot,'err_up':g_up_phot,'err_lo':g_lo_phot,'ra':g_raList,'dec':g_decList,'dist':g_distList,'vcorr':g_vcorrList}
    d_i = {'phot':i_phot,'err_up':i_up_phot,'err_lo':i_lo_phot,'ra':i_raList,'dec':i_decList,'dist':i_distList,'vcorr':i_vcorrList}
    d_r = {'phot':r_phot,'err_up':r_up_phot,'err_lo':r_lo_phot,'ra':r_raList,'dec':r_decList,'dist':r_distList,'vcorr':r_vcorrList}
    d_v = {'phot':v_phot,'err_up':v_up_phot,'err_lo':v_lo_phot,'ra':v_raList,'dec':v_decList,'dist':v_distList,'vcorr':v_vcorrList}
    d_gt = {'Bmedian':gt_BmedianList,'BLstar':gt_BLstarList,'ra':gt_raList,'dec':gt_decList,'dist':gt_distList,'vcorr':gt_vcorrList}
    d_all = {'ra':all_raList,'dec':all_decList,'dist':all_distList,'vcorr':all_vcorrList}
    d_full = {'phot':full_phot,'err_up':full_errUp,'err_lo':full_errLo,'name':full_name,'source':full_source,'ra':full_raList,'dec':full_decList,'dist':full_distList,'vcorr':full_vcorrList}
                
    # add them all to the final dictionary
    d = {'d_b_j':d_b_j,\
    'd_b_all':d_b_all,\
    'd_j':d_j,\
    'd_k':d_k,\
    'd_h':d_h,\
    'd_g':d_g,\
    'd_i':d_i,\
    'd_r':d_r,\
    'd_v':d_v,\
    'd_gt':d_gt,\
    'd_all':d_all,\
    'd_full':d_full}  

    '''

    
    # open and read the galaxy table    
    if getpass.getuser() == 'frenchd':
#         photFilename = '/Users/frenchd/Research/GT_update2/pickleGT.p'
#         photFilename = '/Users/frenchd/Research/inclination/fullPhotPickle.p'

#         photFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filtered.p'
        photFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'
        galaxyFilename = '/Users/frenchd/Research/gt/FinalGalaxyTable13_filtered.csv'
        saveDirectory = '/Users/frenchd/Research/GT_update2/galaxyTable_paper/figures/'
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # open the galaxy file
    theFile = open(galaxyFilename,'rU')
    reader = csv.DictReader(theFile)
    
    # open the pickle file to get at the compiled photometry file
    photFile = open(photFilename,'rU')
    d = pickle.load(photFile)
    photFile.close()
    

#     d1 = d['d_gt']
#     print d1.keys()
#     print
# 
#     lstar = d1['BLstar']
#     count = 0
#     for l in lstar:
#         if l !='x':
#             count +=1
#     print 'count: ',count
#     print
    
    
    lstar = d['Lstar']
    count = 0
    for l in lstar:
        if l != -99.99:
            count +=1
    print 'count: ',count
    print
    
    # what to make?
    makeHistogramLstar_split = True
    
    
##########################################################################################
    if makeHistogramLstar_split:
        # plot a CDF for Lstar values
#         dataset = 'd_gt'
#         d1 = d[dataset]
        d1 = d
        
        # velocity limits
        vlo = 0
        vhi = 10000
        
        save = True

        # make histogram plot of Lstar values
        make_histogram_lstar_split(d1,vlo,vhi,saveDirectory,save)
    

    # close the galaxy table
    theFile.close()



if __name__=="__main__":
    main()
    
    