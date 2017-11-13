#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: GTupdate2_diameterCompare3.py, v 3.0 03/17/2017

Determine how distance measurements compare in order to convert between measurements

v3: remove duplicates, only selecting the largest value when there are multiple measurements
from the same telescope available. This is assuming that the different values indicate
different percentages of total light included.


updated for GTupdate2 (3/17/2017)

"""

import sys
import os
import csv
import string
from math import *
import numpy
from pylab import *
import getpass
import scipy.optimize as optimization
from utilities import *

    
################################################################

        
def isNumber_andGreater(s,n):
    try:
        s=float(s)
        if s >= n:
            return True
        else:
            return False
    except Exception,e:
        return False

def round_to_sig(x,sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)
    
    
def isOdd(val,med,avg,dif,num):
    # returns true if:
    #
    #       abs(val-med) > val*dif  AND  abs(val - avg) > val*dif
    
    if isNumber(val) and isNumber(med) and isNumber(avg) and isNumber(dif):
        if abs(float(val) - float(med)) > float(val) * float(dif):
            if abs(float(val) - float(avg)) > float(val) * float(dif):
                return True
            else:
                return False
        else:
            return False
    else:
        return False
    
    
def isOdd_old(val,med,avg,dif,num):
    # sigma clipping.
    # returns true if:
    #
    #       abs(val-med) > val*dif  AND  abs(val - avg) > val*dif
    
    if isNumber(val) and isNumber(med) and isNumber(avg) and isNumber(dif):
        if sqrt((1.0/num) *(abs(float(val) - float(med)))**2) > float(dif):
            if sqrt((1.0/num) *(abs(float(val) - float(avg)))**2) > float(dif):
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def lineFunc(x,m,b):
    # return y = m*x+b
    return m*x + b
    
    
    
def main():
    # This function reformats Bart's file
    
    # hubble constant used throughout
    hubbleC = 71
    
    # open the files
    if getpass.getuser() == 'frenchd':
        filename = '/usr/data/moosejaw/frenchd/GT_update2/returnPA11111111111111.csv'
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    diameterFile = open(filename,'rU')        
    
    # read in the csv files as dictionary csv files
    reader = csv.DictReader(diameterFile)

    count = 0
    dict = {}
    
#     surveys = [\
#     'r (SDSS Petrosian)',\
#     'r (SDSS Exponential)',\
#     'r (SDSS Isophotal)',\
#     'K_s (2MASS isophotal)',\
#     'K_s (2MASS "total")',\
#     'ESO-LV "Quick Blue" IIa-O',\
#     'POSS1 103a-O',\
#     'r (SDSS deVaucouleurs)',\
#     'r (SDSS de Vaucouleurs)',\
#     'RC3 D_25, R_25 (blue)',\
#     'RC3 D_0 (blue)',\
#     'ESO-Uppsala "Quick Blue" IIa-O',\
#     'ESO-LV IIIa-F',\
#     'RC3 A_e (Johnson B)',\
#     'POSS1 103a-E',\
#     'K_s (LGA/2MASS isophotal)',\
#     'K_s (LGA/2MASS "total")',\
#     'B (Johnson)',\
#     'R (Kron-Cousins)',\
#     '408 & 1407 MHz',\
#     '1407 MHz']
    
    
    # survey 2 is always the same, we're comparing survey1 to 'K_s (2MASS "total")' each time
    survey2 = 'K_s (2MASS "total")'
#     survey2 = 'r (SDSS Petrosian)'

    s2List = []
    
#     survey1 = 'ESO-LV "Quick Blue" IIa-O'
#     survey1 = 'r (SDSS Petrosian)'
#     survey1 = 'r (SDSS Exponential)'
#     survey1 = 'r (SDSS Isophotal)'
#     survey1 = 'K_s (2MASS isophotal)'
#     survey1 = 'K_s (2MASS "total")'
    survey1 = 'ESO-LV "Quick Blue" IIa-O'
#     survey1 = 'POSS1 103a-O'
#     survey1 = 'r (SDSS deVaucouleurs)'
#     survey1 = 'r (SDSS de Vaucouleurs)'
#     survey1 = 'RC3 D_25, R_25 (blue)'
#     survey1 = 'RC3 D_0 (blue)'
#     survey1 = 'ESO-Uppsala "Quick Blue" IIa-O'
#     survey1 = 'ESO-Uppsala "Quick Blue" IIa-O' - none w/ r Petro
#     survey1 = 'ESO-LV IIIa-F' - none w/ r Petro
#     survey1 = 'ESO-LV IIIa-F'
#     survey1 = 'RC3 A_e (Johnson B)'
#     survey1 = 'POSS1 103a-E' - none w/ r Petro
#     survey1 = 'POSS1 103a-E'
#     survey1 = 'K_s (LGA/2MASS isophotal)' - only 1 w/ K total
#     survey1 = 'K_s (LGA/2MASS "total")' - only 1 w/ K total
#     survey1 = 'B (Johnson)'
#     survey1 = 'R (Kron-Cousins)'
#     survey1 = '408 & 1407 MHz' - none w/ r Petro
#     survey1 = '408 & 1407 MHz' - none w/ K total
#     survey1 = '1407 MHz' - only 1 w/ r Petro
#     survey1 = '1407 MHz' - only 1 w/ K total


    s1List = []
    
    # true for major axis comparison, false for minor axis comparison
    major = True
    dif = 0.4
    cutoff = 2000000
    for line in reader:
        count+=1
        if count <455000000:
            oldName = line['oldName']
            oldCompleteList = eval(line['completeList'])
            
            # switch major and minor if they are mixed up.
            completeList = []
            if len(str(oldCompleteList))>5:
#                 print 'completeList: ',oldCompleteList
                for m in oldCompleteList:
                    diameters = m[1]
                    ratio = m[2]
                    maj,min = diameters[0],diameters[1]
                
                    if isNumber(maj):
                        if float(maj) < cutoff:
                            # switch them if the minor axis is larger than the major
                            if isNumber(maj) and isNumber(min):
                                if maj < min:
                                    min = float(maj)*float(ratio)
                        
                            if isNumber(min) and not isNumber(maj):
                                maj,min = min, maj
                
                            completeList.append([m[0],(maj,min),m[2],m[3]])

            if len(completeList)>1 and str(completeList).find(survey1)!=-1 and\
            str(completeList).find(survey2)!=-1:
            
                # first choose only one of the measurements from each telescope if there
                # are multiple available (choose the largest)
                innerDict = {}
                for m in completeList:
                    survey = m[0]
                    diameters = m[1]
                    maj,min = diameters[0],diameters[1]
                    
                    if major:
                        axis = maj
                    else:
                        axis = min
                
                    # check if the survey has been added previously
                    if isNumber(axis):
                        if float(axis) >0:
                            if innerDict.has_key(survey):
                                i = innerDict[survey]
                                prevMaj,prevMin = i[1][0],i[1][1]
                        
                                if major:
                                    prevAxis = prevMaj
                                else:
                                    prevAxis = prevMin
                        
                                if isNumber(prevAxis) and isNumber(axis):
                                    if prevAxis < axis:
                                        innerDict[survey] = m
                    
                            else:
                                innerDict[survey] = m
            
                allSurveys = innerDict.keys()
                uniqueList = []
                for survey in allSurveys:
                    m = innerDict[survey]
                    uniqueList.append(m)         
            
#                 print 'uniqueList: ',uniqueList
                
                len1before = len(s1List)
                len2before = len(s2List)
                allList = []         
                for m in uniqueList:
                    if major:
                        allList.append(m[1][0])
                    else:
                        allList.append(m[1][1])
                
#                 print 'allList: ',allList
                med = median(allList)
                val1 = 'x'
                val2 = 'x'
                avg = average(allList)
                num = len(allList)
                for m in uniqueList:
                    survey = m[0]
                    diameters = m[1]
                    ratio = m[2]
                    pa = m[3]
                    
                    if survey == survey1 and len(s1List) == len1before:
#                         if float(diameters[0]) > 53 and float(diameters[0]) <55:
#                             print 'Oldname: ',oldName
#                             print completeList
#                             print
                        if major:
                            if isOdd(diameters[0],med,avg,dif,num):
                                print 'odd: ',m,'  --  ',med
                                print 'for: ',oldName
                                print
                                break
                            else:
                                val1 = diameters[0]
                        else:
                            if isOdd(diameters[1],med,avg,dif,num):
                                print 'odd: ',m,'  --  ',med
                                print 'for: ',oldName
                                print
                                break
                            else:
                                val1 = diameters[1]
                    
                    # add this to the list if it has not already been added
                    if survey == survey2 and len(s2List) == len2before:
                        if major:
                            if isOdd(diameters[0],med,avg,dif,num):
                                print 'odd: ',m,'  --  ',med
                                print 'for: ',oldName
                                print
                                break
                            else:
                                val2 = diameters[0]
                        else:
                            if isOdd(diameters[1],med,avg,dif,num):
                                print 'odd: ',m,'  --  ',med
                                print 'for: ',oldName
                                print
                                break
                            else:
                                val2 = diameters[1]
                                
                if isNumber(val1) and isNumber(val2):
                    s1List.append(val1)
                    s2List.append(val2)


    print 'len(1): ',len(s1List)
    print 'len(2): ',len(s2List)
    
    max1 = max(s1List)
    max2 = max(s2List)
    maxAll = max(max1,max2)
    print 'maxAll: ', maxAll
    
    fig = figure()
    ax = fig.add_subplot(111)
#     ax.scatter(s1List,s2List)
    f = np.polyfit(s1List, s2List,1,full=True)
    fit, residuals, rank, singular_values, rcond = f
    print
    print 'fit: ',fit
    print 'residuals: ',residuals
    print 'rank: ',rank
    print 'singular_values: ',singular_values
    print 'rcond: ',rcond
    print
    z = np.poly1d(fit)
    xp = np.linspace(0, max1, 100)
    ax.plot(s1List,s2List,'.',label='Data with dif={0}'.format(dif))
    ax.plot(xp,z(xp),'-',label='1D fit: y = {0}*x + {1}'.format(str(round(fit[0],3)),str(round(fit[1],3))))
#     ax.plot(xp, lineFunc(xp,3.95257,-56.5652),'--',c='black',label='1D fit: y = {0}*x + {1}'.format(3.95257,-56.5652))

    # round stuff for nice formatting
    res = round(residuals[0],0)
    singVals = []
    for i in singular_values:
        i2 = round(i,4)
        singVals.append(i2)
    
    rcond2 = round_to_sig(rcond,sig=4)
    
    xLabelPos = 0.6
    yLabelPos = 0.7
    yShift = 0.05
    lsize = 12
    fitLabel = '1D fit: y = {0}*x + {1}'.format(str(round(fit[0],3)),str(round(fit[1],3)))
    residualsLabel = 'Res = {0}'.format(res)
    rankLabel = 'rank = {0}'.format(rank)
    singularValueLabel = 'Singular vals = {0}'.format(singVals)
    rcondLabel = 'rcond = {0}'.format(rcond2)
    
    ax.annotate(fitLabel,xy=(xLabelPos,yLabelPos),xycoords='axes fraction',size=lsize)
    ax.annotate(residualsLabel,xy=(xLabelPos,yLabelPos-yShift),xycoords='axes fraction',size=lsize)
    ax.annotate(rankLabel,xy=(xLabelPos,yLabelPos-yShift-yShift),xycoords='axes fraction',size=lsize)
    ax.annotate(singularValueLabel,xy=(xLabelPos,yLabelPos-yShift-yShift-yShift),xycoords='axes fraction',size=lsize)   
    ax.annotate(rcondLabel,xy=(xLabelPos,yLabelPos-yShift-yShift-yShift-yShift),xycoords='axes fraction',size=lsize)   

#     legend()
    
    ax.set_xlabel('{0} diameters'.format(survey1))
    ax.set_ylabel('{0} diameters'.format(survey2))
    ax.set_ylim(0,max2+max2*0.3)
    ax.set_xlim(0,max1+max1*0.3)
    
#     ans = raw_input("Save (y/n)? ")
#     while ans != 'y' and ans != 'n':
#         ans = raw_input("Save (y/n)? ")
    
    survey1 = survey1.replace('/','_')
    survey2 = survey2.replace('/','_')
    ans = 'y'
    if ans == 'y':
        fig.savefig('/usr/data/moosejaw/frenchd/GT_update2/diameterFits/{0}-{1}_dif{2}.pdf'.format(survey1,survey2,dif),format='pdf')
#         fig.savefig('/usr/users/frenchd/Galaxy Table code/diameterFix/finalFits/ESO-LV_"Quick Blue"_IIa-O_vs_K_s(LGA_2MASS)_total.pdf',format='pdf')
    else:
        show()


    diameterFile.close()

if __name__=="__main__":
    main()
