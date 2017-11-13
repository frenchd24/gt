#!/usr/bin/env python


'''

Diameter chooser for GTupdate2

Tester is under main()

'''

from utilities import *
from scipy import stats
import numpy as np

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


def is_outlier(points, thresh=2.5):
    """
    
    Taken from: https://github.com/joferkington/oost_paper_code/blob/master/utilities.py
    
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    points = np.array(points)
    
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    
    # deal with scenarios like diff = [11, 0 ,0] -> median is 0
    if med_abs_deviation == 0.0:
        med_abs_deviation = np.mean(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def is_val_outlier(val,points,thresh = 3.5):
    # applies is_outlier(points,thresh=3.5) and returns only the outlier result of 'val'
    
    outliers = is_outlier(points,thresh=thresh)
    pIndex = list(points).index(val)
    val_outlier = outliers[pIndex]
    
    return val_outlier


def choose_diameter_measurement(options, piecemeal=True,verbose=False):
    # choose which measurement of diameter to use from the list 'options'
    # returns the list [key,(major,minor),dRatio,pa] for the chosen measurement
    # 
    # if piecemeal == True, the returned diameter array is constructed from
    # the highest ranking each of major, ratio, and pa. If False, instead is returned
    # the highest ranking, MOST COMPLETE diameter set
        
    order = ['K_s (2MASS "total")',\
    'K_s (LGA/2MASS "total")',\
    'K_s (2MASS isophotal)',\
    'K_s (LGA/2MASS "isophotal")',\
    'POSS1 103a-O',\
    'POSS1 103a-E',\
    'ESO-LV "Quick Blue" IIa-O',\
    'r (SDSS Isophotal)',\
    'RC3 D_0 (blue)',\
    'RC3 D_25, R_25 (blue)',\
    'r (SDSS Petrosian)',\
    'r (SDSS deVaucouleurs)',\
    'r (SDSS de Vaucouleurs)',\
    'RC3 A_e (Johnson B)',\
    'r (SDSS Exponential)',\
    'B (Johnson)',\
    'R (Kron-Cousins)',\
    'ESO-Uppsala "Quick Blue" IIa-O',\
    'ESO-LV IIIa-F',\
    ]
    
    # try evaluating options
    try:
        options = eval(options)
    except Exception, e:
        sys.stdout.write('\n eval(options) failed on {0}, continuing...\n'.format(options))

    corrected_options = options
    secondary = options
    
    # the stye is nested dictionaries, the top level is the survey
    #
    # style: d_all = {'K_s (2MASS "total")':{'major':[51],'ratio':[0.5],'pa':[35]},
    #                 ''POSS1 103a-O':{'major':[32.5,28.5,50.5],'ratio':[0.5,0.5,0.4],'pa':[35,24,30]}}
    #
    
    pa_all = []
    ratio_all = []
    major_all = []
    d_all = {}
    majorThresh = 5.5
    ratioThresh = 2.0
    paThresh = 3.5
    
    for c in corrected_options:
        if len(c)>1:
            survey = c[0]
            major,minor = c[1]
            ratio = c[2]
            pa = c[3]
        
            if isNumber(major):
                major_all.append(major)
        
            if isNumber(ratio):
                if float(ratio) >1.0:
                    ratio = 1/float(ratio)
                    
                ratio_all.append(ratio)
        
            if isNumber(pa):
                pa_all.append(pa)
        
            if d_all.has_key(survey):
                prev = d_all[survey]
            
                prev['major'].append(major)
                prev['ratio'].append(ratio)
                prev['pa'].append(pa)
            
            else:
                d_all[survey] = {'major':[major],'ratio':[ratio],'pa':[pa]}  
    
    
    majorSurvey = None
    ratioSurvey = None
    paSurvey = None
    majorValue = None
    ratioValue = None
    paValue = None
    
    for survey in order:
        # only proceed if any of these things have not yet been set
        if not majorSurvey or not ratioSurvey or not paSurvey:
            if d_all.has_key(survey):
                majors = d_all[survey]['major']
                ratios = d_all[survey]['ratio']
                pas = d_all[survey]['pa']
            
                majorsNum = []
                ratiosNum = []
                pasNum = []
                for m in majors:
                    if isNumber(m):
                        majorsNum.append(m)
            
                for r in ratios:
                    if isNumber(r):
                        ratiosNum.append(r)
            
                for p in pas:
                    if isNumber(p):
                        pasNum.append(p)
            
                if len(majorsNum)>0:
                    # choose the largest major diameter
                    majorsNum.sort(reverse=True)
                    largestMajor = majorsNum.pop(0)
                                  
                    # check against all for outliers - major
                    isOutlier = is_val_outlier(largestMajor, major_all, thresh=majorThresh)
                    while isOutlier and len(majorsNum)>0:
                        if len(majorsNum)>0:
                            largestMajor = majorsNum.pop(0)
                            isOutlier = is_val_outlier(largestMajor, major_all, thresh=majorThresh)
                        else:
                            isOutlier = True
                
                    # set this as the one if it is not an outlier and has not already been
                    # set
                    if not isOutlier and not majorSurvey:
                        majorSurvey = survey
                        majorValue = largestMajor
                    
                    
                if len(ratiosNum)>0:
                    # choose the median value in the list (rounding low for even lists)
                    ratiosNum.sort(reverse=True)
                    medianRatio = median_low(ratiosNum)
            
                    # check against all for outliers - ratio
                    isOutlier = is_val_outlier(medianRatio, ratio_all, thresh=ratioThresh)
                    while isOutlier and len(ratiosNum)>0:
                        _ = ratiosNum.pop(ratiosNum.index(medianRatio))
                        if len(ratiosNum)>0:
                            medianRatio = median_low(ratiosNum)
                            isOutlier = is_val_outlier(medianRatio, ratio_all, thresh=ratioThresh)

                        else:
                            isOutlier = True

                    # set this as the one if it is not an outlier and has not already been
                    # set
                    if not isOutlier and not ratioSurvey:
                        ratioSurvey = survey
                        ratioValue = medianRatio
            
                if len(pasNum)>0:
                    # choose the median value in the list (rounding low for even lists)
                    pasNum.sort(reverse=True)
                    medianPA = median_low(pasNum)
    
                    # check against all for outliers - pa
                    isOutlier = is_val_outlier(medianPA, pa_all, thresh=paThresh)
                    while isOutlier and len(pasNum)>0:
                        _ = pasNum.pop(pasNum.index(medianPA))
                        if len(pasNum)>0:
                            medianPA = median_low(pasNum)
                            isOutlier = is_val_outlier(medianPA, pa_all, thresh=paThresh)
                            
                        else:
                            isOutlier = True
                            
                    # set this as the one if it is not an outlier and has not already been
                    # set
                    if not isOutlier and not paSurvey:
                        paSurvey = survey
                        paValue = medianPA
    
    
    # now calculate a minor axis
    if isNumber(majorValue) and isNumber(ratioValue):
        minorValue = float(majorValue) * float(ratioValue)
    else:
        minorValue = None

    toReturn = [majorSurvey, (majorValue,minorValue), ratioValue, paValue]
    keys = [majorSurvey, ratioSurvey, paSurvey]
    
    return toReturn, keys
    
    
def main():
    
#     options = [['r (SDSS Isophotal)', (5.35, 4.91), 0.92, 152.0], ['r (SDSS Petrosian)', (1.18, None), None, None], ['r (SDSS deVaucouleurs)', (0.01, 0.01), 1.0, 163.0], ['r (SDSS Exponential)', (0.0, 0.0), 0.81, 163.0]]
#     options = [['r (SDSS Exponential)', (0.0, 0.0), 0.81, 163.0], ['r (SDSS Petrosian)', (1.18, None), None, None], ['r (SDSS deVaucouleurs)', (0.01, 0.01), 1.0, 163.0]]
#     options = [['POSS1 103a-O', (270.0, 33.0), 0.12, None], ['POSS1 103a-O', (1350.0, 372.0), 0.28, None], ['RC3 D_25, R_25 (blue)', (1312.7, 434.5), 0.33, 162.0], ['RC3 A_e (Johnson B)', (476.6, None), None, None], ['RC3 D_0 (blue)', (1312.7, None), None, None], ['K_s (LGA/2MASS isophotal)', (146.4, 146.4), 1.0, 16.0], ['K_s (LGA/2MASS "total")', (250.2, 250.2), 1.0, 16.0], ['POSS1 103a-O', (1380.0, 480.0), 0.35, 162.0]]
#     options = ['x']
#     options = [['POSS1 103a-O', (54.0, 96.0), 1.78, None], ['RC3 D_25, R_25 (blue)', (101.9, 60.02), 0.59, None], ['RC3 A_e (Johnson B)', (42.5, None), None, None], ['RC3 D_0 (blue)', (106.7, None), None, None]]
#     options = [['ESO-LV "Quick Blue" IIa-O', (54.95, 104.87), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (120.23, 229.45), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (142.89, 272.69), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (94.41, 180.17), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (107.15, 204.48), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (124.45, 237.5), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (71.61, 136.66), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (60.26, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (28.43, 108.93), 0.26, 104.0], ['ESO-LV IIIa-F', (25.35, None), None, None]]
#     options = [['ESO-Uppsala "Quick Blue" IIa-O', (96.0, 96.0), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (97.72, 97.72), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (117.49, 117.49), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (134.9, 134.9), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (74.13, 74.13), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (86.1, 86.1), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (103.51, 103.51), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (51.29, 51.29), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (51.29, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (12.15, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (5.34, 8.59), 0.62, 173.0], ['ESO-LV IIIa-F', (5.82, None), None, None], ['RC3 D_25, R_25 (blue)', (95.1, 90.82), 0.95, None], ['RC3 D_0 (blue)', (95.1, None), None, None]]
    options = [['ESO-Uppsala "Quick Blue" IIa-O', (300.0, 60.0), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (197.24, 986.2), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (323.59, 1617.95), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (506.99, 2534.95), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (367.28, 1836.4), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (489.78, 2448.9), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (668.34, 3341.7), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (213.8, 1069.0), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (49.55, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (34.62, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (46.46, 227.75), 0.2, 94.0], ['ESO-LV IIIa-F', (16.03, None), None, None], ['RC3 D_25, R_25 (blue)', (198.7, 41.53), 0.21, 31.0], ['RC3 D_0 (blue)', (198.7, None), None, None]]

    best,best_keys = choose_diameter_measurement(options)
    
    
    print 'best: ',best
    print' best_keys: ',best_keys
    print
    print
    print 'Options were: '
    for i in options:
        print i
    
    print
    print 'done'
    
    
if __name__ == '__main__':
    main()
    