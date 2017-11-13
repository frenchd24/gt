#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: GTupdate2_diameterCompare5.py, v 5.0 4/11/17

Determine how distance measurements compare in order to convert between measurements

v3: remove duplicates, only selecting the largest value when there are multiple measurements
from the same telescope available. This is assuming that the different values indicate
different percentages of total light included.

updated for GTupdate2 (3/17/2017)

v4: change fit to use scipy.stats.linregress instead, so i get errors (3/31/17)

v5: change the fitting again, try some shit I found on the internet here:
https://micropore.wordpress.com/2017/02/07/python-fit-with-error-on-both-axis/

(4/11/17)

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

from scipy import odr


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
    dif = 10
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
    
    fig = figure()
    ax = fig.add_subplot(111)
    
    # do the fit
#     f = np.polyfit(s1List, s2List,1,full=True)
#     fit, residuals, rank, singular_values, rcond = f
#     print
#     print 'fit: ',fit
#     print 'residuals: ',residuals
#     print 'rank: ',rank
#     print 'singular_values: ',singular_values
#     print 'rcond: ',rcond
#     print
    
#     z = np.poly1d(fit)
#     xp = np.linspace(0, max1, 100)
#     ax.plot(xp,z(xp),'-',label='1D fit: y = {0}*x + {1}'.format(str(round(fit[0],3)),str(round(fit[1],3))))

##########################################################################################

#     m, b, r_value, p_value, std_err = linregress(s1List, s2List)

    xdata = s1List
    ydata = s2List
    
#     pfit, pcov = optimize.curve_fit(lineFunc, xdata, ydata)

##########################################################################################

    # These are initial guesses for fits:
    pstart = [1.0,0.0]
    
    # bootstrap method

#     def ff(x, p):
#         return lineFunc(x, *p)

    # These are the true parameters
#     p0 = 1.0
#     p1 = 40
#     p2 = 2.0

    # These are initial guesses for fits:
#     pstart = [
#         p0 + random.random(),
#         p1 + 5.*random.random(), 
#         p2 + random.random()
#     ]

#     %matplotlib inline
#     import matplotlib.pyplot as plt
#     xvals = np.linspace(0., 1, 120)
#     yvals = f(xvals, p0, p1, p2)

    # Generate data with a bit of randomness
#     xdata = np.array(xvals)
#     np.random.seed(42)
#     err_stdev = 0.2
#     yvals_err =  np.random.normal(0., err_stdev, len(xdata))
#     ydata = f(xdata, p0, p1, p2) + yvals_err

#     plt.plot(xvals, yvals)
#     plt.plot(xdata, ydata, 'o', mfc='None')
    
#     pfit, perr = fit_bootstrap(pstart, xdata, ydata, ff)
# 
#     print 'pfit: ',pfit
#     print
#     print 'perr: ',perr

##########################################################################################

#     pfit, perr = fit_curvefit(pstart, xdata, ydata, ff)
# 
#     pfit, pcov = optimize.curve_fit(lineFunc,xdata,ydata,p0=pstart)
#     error = [] 
#     for i in range(len(pfit)):
#         try:
#           error.append(np.absolute(pcov[i][i])**0.5)
#         except:
#           error.append( 0.00 )
#     pfit_curvefit = pfit
#     perr_curvefit = np.array(error)
#     print "curve_fit: ", pfit_curvefit, perr_curvefit 
# 
#     print("\nFit parameters and parameter errors from curve_fit method :")
#     print("pfit = ", pfit)
#     print("perr = ", perr)


#     error = [] 
#     for i in range(len(pfit)):
#         try:
#           error.append(absolute(pcov[i][i])**0.5)
#         except:
#           error.append( 0.00 )
#     pfit_curvefit = pfit
#     perr_curvefit = array(error)

#     print("\nFit parameters and parameter errors from curve_fit method :")
#     print("pfit = ", pfit_curvefit)
#     print("perr = ", perr_curvefit)
    
    
#     pfit_curvefit = pfit
#     perr_curvefit = array(perr)
#     
#     m = pfit_curvefit[0]
#     b = pfit_curvefit[1]
#     
#     std_err = perr_curvefit[0]
#     r_value = 0
#     p_value = 0

##########################################################################################


#     def func(p, x):
#         a, b, c = p
#         return a*x*x + b*x + c

    def func(p, x):
        b, c = p
        return b*x + c

    # Model object
    quad_model = odr.Model(func)
    
    x = np.array(xdata)
    y = np.array(ydata)
    noise_x = x*0.05
    noise_y = y*0.05

    # Create a RealData object
    data = odr.RealData(x, y, sx=noise_x, sy=noise_y)
#     data = odr.Data(x, y)


    # Set up ODR with the model and data
#     odr2 = odr.ODR(data, quad_model, beta0=[0., 1., 1.])
    odr2 = odr.ODR(data, quad_model, beta0=[1., 1.])


    # Run the regression
    out = odr2.run()

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
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr

    x_fit = np.linspace(min(x), max(x), 100)
    fit = func(popt, x_fit)
    fit_up = func(popt_up, x_fit)
    fit_dw= func(popt_dw, x_fit)

    #plot
    fig, ax = plt.subplots(1)
#     rcParams['font.size']= 15
#     errorbar(x, y, yerr=noise_y, xerr=noise_x, hold=True, ecolor='k', fmt='none', label='data')
    xlabel('{0} diameters'.format(survey1), fontsize=15)
    ylabel('{0} diameters'.format(survey2), fontsize=15)
    title('fit with error on both axis', fontsize=15)
    plot(x_fit, fit, 'r', lw=2, label='best fit curve')
#     plot(x0, y0, 'k-', lw=2, label='True curve')
#     ax.fill_between(x_fit, fit_up, fit_dw, alpha=.25, label='5-sigma interval')
    legend(loc='lower right',fontsize=15)
    xlim(0,max(x))
    ylim(0,max(y))
    show()



##########################################################################################

    # try orthoregress modules
    import orthoregress
    
    xerr = np.array(xdata)*0.05
    yerr = np.array(ydata)*0.05
    
    out = orthoregress.orthoregress(xdata,ydata,xerr,yerr)

    print
    print 'out.beta: ',out.beta
    print 'out.sd_beta: ',out.sd_beta
    print 'out.cov_beta: ',out.cov_beta
    print 'out.delta: ',out.delta
    print 'out.eps: ',out.eps
    print 'out.res_var: ',out.res_var
    print 'out.sum_sqare: ',out.sum_sqare
    print 'out.sum_square_delta: ',out.sum_square_delta
    print 'out.sum_square_eps: ',out.sum_square_eps
    print
    sys.exit()

    # round stuff for nice formatting
    m2 = round_to_sig(m,sig=3)
    b2 = round_to_sig(b,sig=3)
    r_value2 = round_to_sig(r_value,sig=3)
    p_value2 = round_to_sig(p_value,sig=3)
    std_err2 = round_to_sig(std_err,sig=3)

    
    # plot the fit
    xp = np.linspace(0, max1, 100)
    ax.plot(xp, m2 * xp + b2, '-', label='1D fit: y = {0}*x + {1}'.format(m2,b2))
    ax.plot(xp, (m2 + std_err) * xp + b2, '-',color='red', label='1D fit: y = {0}*x + {1}'.format(m2,b2))
    ax.plot(xp, (m2 - std_err) * xp + b2, '-',color='red', label='1D fit: y = {0}*x + {1}'.format(m2,b2))

    
    # plot the data
    ax.plot(s1List,s2List,'.',label='Data with dif={0}'.format(dif))

    print
    print 'm2: ',m2
    print 'b2: ',b2
    print 'r_value2: ',r_value2
    print 'p_value2: ',p_value2
    print 'std_err2: ',std_err2
    print

    
    xLabelPos = 0.6
    yLabelPos = 0.7
    yShift = 0.05
    lsize = 12
    fitLabel = '1D fit: y = {0}*x + {1}'.format(m2,b2)
    r_valueLabel = 'r_value = {0}'.format(r_value2)
    p_valueLabel = 'p_value = {0}'.format(p_value2)
    std_errLabel = 'std_err = {0}'.format(std_err2)
    
    ax.annotate(fitLabel,xy=(xLabelPos,yLabelPos),xycoords='axes fraction',size=lsize)
    ax.annotate(r_valueLabel,xy=(xLabelPos,yLabelPos-yShift),xycoords='axes fraction',size=lsize)
    ax.annotate(p_valueLabel,xy=(xLabelPos,yLabelPos-yShift-yShift),xycoords='axes fraction',size=lsize)
    ax.annotate(std_errLabel,xy=(xLabelPos,yLabelPos-yShift-yShift-yShift),xycoords='axes fraction',size=lsize)   

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
        fig.savefig('/usr/data/moosejaw/frenchd/GT_update2/diameterFits/{0}-{1}_dif{2}_linreg3.pdf'.format(survey1,survey2,dif),format='pdf')
#         fig.savefig('/usr/users/frenchd/Galaxy Table code/diameterFix/finalFits/ESO-LV_"Quick Blue"_IIa-O_vs_K_s(LGA_2MASS)_total.pdf',format='pdf')
    else:
        show()


    diameterFile.close()

if __name__=="__main__":
    main()
