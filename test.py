#!/usr/bin/env python

'''

Test median_low() method

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
    print 'median: ',median
    diff = np.sum((points - median)**2, axis=-1)
    print 'diff: ',diff
    diff = np.sqrt(diff)
    print 'diff2 = ',diff
    med_abs_deviation = np.median(diff)
    if med_abs_deviation == 0.0:
        med_abs_deviation = np.mean(diff)
    print 'med_abs_deviation: ',med_abs_deviation

    modified_z_score = 0.6745 * diff / med_abs_deviation
    
    print 'modified_z_score = ',modified_z_score
    
    return modified_z_score > thresh


def is_val_outlier(val,points,thresh = 3.5):
    # applies is_outlier(points,thresh=3.5) and returns only the outlier result of 'val'
    
    outliers = is_outlier(points,thresh=thresh)
    
    pIndex = list(points).index(val)
    
    val_outlier = outliers[pIndex]
    
    return val_outlier
        

def main():
#     a = [0.12, 0.28, 0.33, 1.0, 1.0, 0.35]
#     a = [0.12, 0.28, 0.33, 1.0, 1.0, 0.35, None]
    a = [152.0, 163.0, 163.0]


    val = 152.0
    thresh = 4.5
    
    print median_low(a)
    print
    print 'a = ',a
    print 'outliers: ',is_outlier(a,thresh=thresh)
    print 'val: ',val
    print 'is_val_outlier(val,points,thresh = {0}): '.format(thresh),is_val_outlier(val,a,thresh=thresh)
    

if __name__ == '__main__':
    main()