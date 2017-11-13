#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: GTupdate2_diameterCompare7.py, v 7.0 09/25/17

Determine how distance measurements compare in order to convert between measurements

v3: remove duplicates, only selecting the largest value when there are multiple measurements
from the same telescope available. This is assuming that the different values indicate
different percentages of total light included.

updated for GTupdate2 (3/17/2017)

v4: change fit to use scipy.stats.linregress instead, so i get errors (3/31/17)

v5: change the fitting again, try some shit I found on the internet here:
https://micropore.wordpress.com/2017/02/07/python-fit-with-error-on-both-axis/

(4/11/17)

v6: progressing on v5, use the othoregress.py package to fit using orthogonal distance
regression. (4/13/17)

v7: make errors equal on all axes - doesn't really help (09/25/17)


"""

import sys
import os
import csv
import string
from math import *
from pylab import *
import getpass
from scipy import optimize
from utilities import *
from scipy.stats import linregress

import numpy as np
from scipy.optimize import curve_fit

# from scipy import odr
import orthoregress



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
    try:
        return round(x, sig-int(floor(log10(abs(x))))-1)
    except Exception,e:
        sys.stdout.write('Error: {0}'.format(e))
        return x
        
    
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


def fit_bootstrap(p0, datax, datay, function, yerr_systematic=0.0):

    errfunc = lambda p, x, y: function(x,p) - y

    # Fit first time
    pfit, perr = optimize.leastsq(errfunc, p0, args=(datax, datay), full_output=0)
    print 'first: ',pfit,perr

    # Get the stdev of the residuals
    residuals = errfunc(pfit, datax, datay)
    sigma_res = np.std(residuals)

    sigma_err_total = np.sqrt(sigma_res**2 + yerr_systematic**2)

    # 100 random data sets are generated and fitted
    ps = []
    for i in range(100):

        randomDelta = np.random.normal(0., sigma_err_total, len(datay))
        randomdataY = datay + randomDelta

        randomfit, randomcov = \
            optimize.leastsq(errfunc, p0, args=(datax, randomdataY),\
                             full_output=0)

        ps.append(randomfit) 

    ps = np.array(ps)
    mean_pfit = np.mean(ps,0)

    # You can choose the confidence interval that you want for your
    # parameter estimates: 
    Nsigma = 2. # 1sigma gets approximately the same as methods above
                # 1sigma corresponds to 68.3% confidence interval
                # 2sigma corresponds to 95.44% confidence interval
    err_pfit = Nsigma * np.std(ps,0)

    pfit_bootstrap = mean_pfit
    perr_bootstrap = err_pfit
    print 'pfit_bootstrap, perr_bootstrap : ',pfit_bootstrap, perr_bootstrap
    return pfit_bootstrap, perr_bootstrap




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
    survey1 = 'r (SDSS Petrosian)'
#     survey1 = 'r (SDSS Exponential)'
#     survey1 = 'r (SDSS Isophotal)'
#     survey1 = 'K_s (2MASS isophotal)'
#     survey1 = 'K_s (2MASS "total")'
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
    dif = 10000
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
                    maj,minor = diameters[0],diameters[1]
                
                    if isNumber(maj):
                        if float(maj) < cutoff:
                            # switch them if the minor axis is larger than the major
                            if isNumber(maj) and isNumber(minor):
                                if maj < minor:
                                    minor = float(maj)*float(ratio)
                        
                            if isNumber(minor) and not isNumber(maj):
                                maj,minor = minor, maj
                
                            completeList.append([m[0],(maj,minor),m[2],m[3]])

            if len(completeList)>1 and str(completeList).find(survey1)!=-1 and\
            str(completeList).find(survey2)!=-1:
            
                # first choose only one of the measurements from each telescope if there
                # are multiple available (choose the largest)
                innerDict = {}
                for m in completeList:
                    survey = m[0]
                    diameters = m[1]
                    maj,minor = diameters[0],diameters[1]
                    
                    if major:
                        axis = maj
                    else:
                        axis = minor
                
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
    

##########################################################################################

    xdata = s1List
    ydata = s2List
    
##########################################################################################

    def func(p, x):
        b, c = p
        return b*x + c

    # Run the regression
#     xerr = np.array(xdata)*0.01
#     yerr = np.array(ydata)*0.01
    
    xerr = np.ones(len(xdata))
    yerr = np.ones(len(ydata))
    
     
    out = orthoregress.orthoregress(xdata,ydata,xerr,yerr)
    m_lr, b_lr, r_value, p_value, std_err = linregress(s1List, s2List)

    print
    print 'out.beta: ',out.beta
    print 'out.sd_beta: ',out.sd_beta

    #print fit parameters and 1-sigma estimates
    popt = out.beta
    perr = out.sd_beta
    print 
    print('fit parameter 1-sigma error')
    print
    print('---------')
    for i in range(len(popt)):
        print(str(popt[i])+' +- '+str(perr[i]))

    # prepare confidence level curves
    # to draw 5-sigma intervals
    nstd = 5.
    popt_up5 = popt + nstd * perr
    popt_dw5 = popt - nstd * perr

    # and 1-sigma intervals
    popt_up1 = popt + perr
    popt_dw1 = popt - perr
    
    x_fit = np.linspace(min(xdata), max(xdata), 100)
    fit = func(popt, x_fit)
    fit_up5 = func(popt_up5, x_fit)
    fit_dw5 = func(popt_dw5, x_fit)
    
    fit_up1 = func(popt_up1, x_fit)
    fit_dw1 = func(popt_dw1, x_fit)

    # linregress result
    fit_lr = func([m_lr,b_lr],x_fit)


##########################################################################################
    

    # round stuff for nice formatting
    m_lr2 = round_to_sig(m_lr,sig=3)
    b_lr2 = round_to_sig(b_lr,sig=3)
    r_value2 = round_to_sig(r_value,sig=3)
    p_value2 = round_to_sig(p_value,sig=3)
    std_err2 = round_to_sig(std_err,sig=3)

    m_odr = round_to_sig(popt[0],sig=3)
    b_odr = round_to_sig(popt[1],sig=3)
    merr_odr = round_to_sig(perr[0],sig=3)
    berr_odr = round_to_sig(perr[1],sig=3)
    
    # plot the fit
#     xp = np.linspace(0, max1, 100)
#     ax.plot(xp, m2 * xp + b2, '-', label='1D fit: y = {0}*x + {1}'.format(m2,b2))
#     ax.plot(xp, (m2 + std_err) * xp + b2, '-',color='red', label='1D fit: y = {0}*x + {1}'.format(m2,b2))
#     ax.plot(xp, (m2 - std_err) * xp + b2, '-',color='red', label='1D fit: y = {0}*x + {1}'.format(m2,b2))
    
    
    fitLabel = '1D fit: y = {0}*x + {1}'.format(m_lr2,b_lr2)
    r_valueLabel = 'r_value = {0}'.format(r_value2)
    p_valueLabel = 'p_value = {0}'.format(p_value2)
    std_errLabel = 'std_err = {0}'.format(std_err2)
    
    #plot
    fig = figure(figsize=(7.7,5.7))
    ax = fig.add_subplot(111)
#     fig, ax = plt.subplots(1)
    fsize = 13
    rcParams['font.size']= fsize
    errorbar(xdata, ydata, yerr=yerr, xerr=xerr, hold=True, ecolor='k', fmt='none', label='data')
    xlabel('{0} diameters'.format(survey1), fontsize=fsize)
    ylabel('{0} diameters'.format(survey2), fontsize=fsize)
    title('fit with error on both axis', fontsize=fsize)
    plot(x_fit, fit, 'r', lw=2, label='ODR = {0}*x + {1}'.format(m_odr,b_odr))
    plot(x_fit, fit_lr, color='green', lw=2, label='Linregress fit curve')
    ax.fill_between(x_fit, fit_up1, fit_dw1, alpha=.45, label='1-sigma: m+/- {0}, b+/- {1}'.format(merr_odr,berr_odr))
    ax.fill_between(x_fit, fit_up5, fit_dw5, alpha=.20, label='5-sigma interval')

    legend(loc='upper left',fontsize=fsize)
#     xlim(0,max(x))
#     ylim(0,max(y))

    print
    print 'm_lr2: ',m_lr2
    print 'b_lr2: ',b_lr2
    print 'r_value2: ',r_value2
    print 'p_value2: ',p_value2
    print 'std_err2: ',std_err2
    print
    
    ax.set_ylim(0,max2+max2*0.35)
    ax.set_xlim(0,max1+max1*0.1)
    
#     show()
#     sys.exit()
    
    survey1 = survey1.replace('/','_')
    survey2 = survey2.replace('/','_')
    ans = 'y'
    if ans == 'y':
        directory = '/usr/data/moosejaw/frenchd/GT_update2/diameterFits_odr3/'
        
        # write the figure
        fig.savefig('{0}/{1}-{2}_odrreg.pdf'.format(directory,survey1,survey2),format='pdf')
        
        # write the information file
        fit_result = open('{0}/{1}-{2}_odrreg_fits.txt'.format(directory,survey1,survey2),'wt')
        
        fit_result.write('----Orthogonal Distance Regression results----\n')
        fit_result.write('[m,b] = {0}\n'.format(out.beta))
        fit_result.write('[m_err,b_err] = {0}\n'.format(out.sd_beta))
        fit_result.write('Covariance = {0}\n'.format(out.cov_beta))
        fit_result.write('Residual Variance = {0}\n'.format(out.res_var))
        fit_result.write('Relative Error = {0}\n'.format(out.rel_error))
        fit_result.write('Info = {0}\n'.format(out.info))
        fit_result.write('Stop Reason = {0}\n'.format(out.stopreason))
        fit_result.write('\n')
        fit_result.write('\n')
        fit_result.write('----Standard Linear Regression results----\n')
        fit_result.write('{0}\n'.format(fitLabel))
        fit_result.write('r_value = {0}\n'.format(r_valueLabel))
        fit_result.write('p_value = {0}\n'.format(p_valueLabel))
        fit_result.write('std_err = {0}\n'.format(std_err))

        fit_result.close()
        
#         fig.savefig('/usr/users/frenchd/Galaxy Table code/diameterFix/finalFits/ESO-LV_"Quick Blue"_IIa-O_vs_K_s(LGA_2MASS)_total.pdf',format='pdf')
    else:
        show()

    diameterFile.close()

if __name__=="__main__":
    main()
