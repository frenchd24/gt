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

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def choose_diameter_measurement(options, piecemeal=True,verbose=False):
    # choose which measurement of diameter to use from the list 'options'
    # returns the list [key,(major,minor),dRatio,pa] for the chosen measurement
    # 
    # if piecemeal == True, the returned diameter array is constructed from
    # the highest ranking each of major, ratio, and pa. If False, instead is returned
    # the highest ranking, MOST COMPLETE diameter set
    
    manual_choose = True
    
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
    
    # correct weirdness
#     corrected_options = []
#     for n in options:
#         if len(n)>1:
#             major,minor = n[1]
#             ratio = n[2]
#             if isNumber(major) and isNumber(minor):
#                 if float(minor) > float(major):
#                     major, minor = minor, major
#                 
#             if isNumber(ratio):
#                 if ratio >1.0:
#                     ratio = 1.0/ratio
#                 
#             corrected_options.append([n[0],(major,minor),ratio,n[3]])
#         else:
#             corrected_options.append(n)

#     corrected_options = []
#     secondary = []
#     for n in options:
#         if len(n)>1:
#             major,minor = n[1]
#             ratio = n[2]
#             if isNumber(major) and isNumber(minor):
#                 if float(major) > float(minor):
#                     corrected_options.append(n)
#                 else:
#                     secondary.append(n)
#             else:
#                 secondary.append(n)
#         else:
#             corrected_options.append(n)
#     
#     if len(corrected_options) == 0:
#         corrected_options = secondary
#     
#     hasMajor = False
#     hasMinor = False
#     hasRatio = False
#     hasPA = False
#     for c in corrected_options:
#         major,minor = c[1]
#         ratio = c[2]
#         pa = c[3]
#         if isNumber(major) and not hasMajor:
#             hasMajor = True
#         if isNumber(minor) and not hasMinor:
#             hasMinor = True
#         if isNumber(ratio) and not hasRatio:
#             hasRatio = True
#         if isNumber(pa) and not hasPA:
#             hasPA = True
#             
#     if not hasMajor:
#         for s in secondary:
#             major,minor = s[1]
#             ratio = s[2]
#             pa = s[3]
#             if isNumber(major):
#                 corrected_options.append(s)
#                 
#     if not hasMinor:
#         for s in secondary:
#             major,minor = s[1]
#             ratio = s[2]
#             pa = s[3]
#             if isNumber(minor):
#                 corrected_options.append(s)
# 
#     if not hasRatio:
#         for s in secondary:
#             major,minor = s[1]
#             ratio = s[2]
#             pa = s[3]
#             if isNumber(ratio):
#                 corrected_options.append(s)
#                 
#     if not hasPA:
#         for s in secondary:
#             major,minor = s[1]
#             ratio = s[2]
#             pa = s[3]
#             if isNumber(pa):
#                 corrected_options.append(s)

    corrected_options = options
    secondary = options
    
    ratio_orderFull = []
    ratio_onlyFull = []
    pa_orderFull = []
    pa_onlyFull = []
    d_ratio = {}
    d_pa = {}
    for c in corrected_options:
        survey = c[0]
        ratio = c[2]
        pa = c[3]
        
        if isNumber(ratio):
            ratio_orderFull.append([ratio,survey])
            ratio_onlyFull.append(ratio)
            
            # add to the dictionary
            if d_ratio.has_key(survey): 
                n = d_ratio[survey]
                n.append(ratio)
                d_ratio[survey] = n
            else:
                d_ratio[survey] = [ratio]
    
        if isNumber(pa):
            pa_orderFull.append([pa,survey])
            pa_onlyFull.append(pa)
            
            # add to the dictionary
            if d_pa.has_key(survey): 
                n = d_pa[survey]
                n.append(pa)
                d_pa[survey] = n
            else:
                d_pa[survey] = [pa]

    allowUnityRatio = True
    
    print 'd_ratio: ',d_ratio
    print 'd_pa: ', d_pa

#     print 'corrected_options: ',corrected_options
    
    # first prepare the list. If there are multiple surveys under the same name,
    # keep only the best
    d = {}
    for n in corrected_options:
        key = n[0]        
        
        # if you encounter a duplicate survey, take the best pieces
        if d.has_key(key):
            other = d[key]
            
            otherMajor,otherMinor = other[1]
            otherRatio = other[2]
            otherPA = other[3]
            
            major,minor = n[1]
            ratio = n[2]
            pa = n[3]
                    
            # now the main thing
            if isNumber(major) and isNumber(otherMajor):
                if float(major) > float(otherMajor):
                    if not isNumber(pa) and isNumber(otherPA):
                        sys.stdout.write("\n Largest not complete: {0} vs {1}\n".format(n, other))
                elif float(otherMajor) >float(major):
                    if not isNumber(otherPA) and isNumber(pa):
                        sys.stdout.write("\n Largest not complete: {0} vs {1}\n".format(n, other))

            
            # best case, just choose the larger measurement
            if isNumber(major) and isNumber(minor) and isNumber(pa):
                if isNumber(otherMajor) and isNumber(otherMinor) and isNumber(otherPA):
                    if major > minor or allowUnityRatio:
                        if float(major) > float(otherMajor):
                            # if this one is better than the last one, replace it
                            d[key] = n
                        elif float(otherMajor) == float(otherMinor):
                            # if the current major axis is not the larger, but is 
                            # different from the minor axis value whereas the other
                            # survey's is not, replace it
                            d[key] = n
                    
                    # if they're both equal, choose the bigger one
                    elif float(otherMajor) == float(otherMinor):
                        if float(major) > float(otherMajor):
                            d[key] = n
                
                # if not isNumber(*) for all the things, choose the new one
                else:
                    # if piecemeal take all the biggest pieces
                    if piecemeal:
                        if isNumber(major) and isNumber(otherMajor):
                            if float(major) > float(otherMajor):
                                new_n = d[key]
                                new_n[1] = (major,minor)
                                d[key] = new_n
                        # if the other one sucks
                        elif isNumber(major):
                            new_n = d[key]
                            new_n[1] = (major,minor)
                            d[key] = new_n
                                                            
                        if isNumber(ratio) and isNumber(otherRatio):
                            if float(otherRatio) >=1.0 or float(otherRatio) <=0.0:
                                new_n = d[key]
                                new_n[2] = ratio
                                d[key] = new_n
                                
                        elif isNumber(ratio):
                                new_n = d[key]
                                new_n[2] = ratio
                                d[key] = new_n
                                
                        if isNumber(pa) and not isNumber(otherPA):
                            new_n = d[key]
                            new_n[3] = pa
                            d[key] = new_n
                            
                    # if not piecemeal, take the complete one
                    else:
                        d[key] = n
            
            # if the current one sucks, see how much the old one sucks
            elif isNumber(major) and isNumber(otherMajor):
                # both have major
                if isNumber(minor) and isNumber(otherMinor):
                    # both have minor - hence not PA
                    if float(major) > float(otherMajor):
                        new_n = d[key]
                        new_n[1] = (major,minor)
                        d[key] = new_n
                    
                    # now check PA
                    if isNumber(pa) and not isNumber(otherPA):
                        new_n = d[key]
                        new_n[3] = pa
                        d[key] = new_n
                    
                    elif not isNumber(pa) and isNumber(otherPA):
                        new_n = d[key]
                        new_n[3] = otherPA
                        d[key] = new_n
                    
                elif isNumber(minor) and not isNumber(otherMinor):
                    # OTHER does not have minor
                    if isNumber(pa) and not isNumber(otherPA):
                        # if the other one just sucks, choose this one
                        d[key] = n
                    else:
                        # if OTHER has a PA and THIS doesn't, but otherwise
                        # THIS is better, just grab the PA from OTHER
                        new_n = n
                        new_n[3] = otherPA
                        d[key] = new_n
                    
                # if THIS one has a PA, but not a minor, take the PA
                elif isNumber(otherMinor) and not isNumber(minor):
                    if isNumber(pa) and not isNumber(otherPA):
                        new_n = other
                        new_n[3] = pa
                        d[key] = new_n
                        
            elif isNumber(major) and not isNumber(otherMajor):
                # this one has a major, the other does not
                d[key] = n
                
                # if the other had a PA and this one does not, take the PA
                if not isNumber(pa) and isNumber(otherPA):
                    new_n = d[key]
                    new_n[3] = otherPA
                    d[key] = new_n
            
            # if piecemeal, grab the median PA and ratio
            if piecemeal:
                ratios = d_ratio[key]
                pas = d_pa[key]
                
                new_n = d[key]

                if len(ratios) >0:
                    med_ratio = median_low(ratios)
                    new_n[2] = med_ratio

                if len(pas) >0:
                    med_pa = median_low(pas)
                    new_n[3] = med_pa

        else:
            # if d does not have this key yet
            d[key] = n
            
    new_options = d.values()
    print 'new_options: ',new_options
    
    found = False
    found_equal_axes = False
    found_zero_axes = False
    found_no_pa = False
    found_no_minor = False
    found_only_major = False
    
    major_order = []
    minor_order = []
    ratio_order = []
    pa_order = []
    
    ratio_only = []
    pa_only = []
    
    ordered_options = []
    
    if verbose:
        print 'd: ',d
    
    for i in order:
        if d.has_key(i):
            ordered_options.append(d[i])
            
            all = d[i]
            survey = all[0]
            major = all[1][0]
            minor = all[1][1]
            ratio = all[2]
            pa = all[3]
            
            if verbose:
                print 'Currently on: ',survey
            
            # record the order of non-zero values for each parameter
            if isNumber(major):
                major_order.append([major,survey])
            
            if isNumber(minor):
                minor_order.append([minor,survey])
            
            if isNumber(ratio):
                ratio_order.append([ratio,survey])
                ratio_only.append(float(ratio))
            
            if isNumber(pa):
                pa_order.append([pa,survey])
                pa_only.append(float(pa))
                
            if not found:
                if isNumber(major) and isNumber(minor) and isNumber(pa):
                    # skip the weird cases where one measurement option has major = minor
                    if float(major) > float(minor) and float(major)>0.0:
                        found = all
                        break
                        
                    elif float(major) == float(minor) and float(major) >0.0:
                        found_equal_axes = all

                    else:
                        found_zero_axes = all
                
                elif isNumber(minor) and not isNumber(pa) and not found_no_pa:
                    found_no_pa = all
            
                elif not isNumber(minor) and isNumber(pa) and not found_no_minor:
                    found_no_minor = all
        
                elif not found_only_major:
                    found_only_major = all

    if verbose:
        print
        print 'major_order: ',major_order
        print 'minor_order: ',minor_order
        print 'ratio_order: ',ratio_order
        print 'pa_order: ',pa_order
        print
        print 'ordered_options: ',ordered_options

    if piecemeal and len(ordered_options) >0:
        major, major_survey_key = major_order[0]
        survey = major_survey_key
    
        if len(ratio_order)>0:
            ratio,ratio_survey_key = ratio_order[0]
        else:
            ratio = None
            ratio_survey_key = None
        
        if len(pa_order)>0:
            pa, pa_survey_key = pa_order[0]
        else:
            pa = None
            pa_survey_key = None

        # if this result has ratio=1.0, find the next highest ranking ratio and pa
        
        # do a KS test to see if the top rank is an outlier
        
        print
        print 'ratio_only: ',ratio_onlyFull
        print 'ratio outliers: ',is_outlier(ratio_onlyFull, thresh=2.5)
        print
        print 'pa_only: ',pa_onlyFull
        print 'pa outliers: ',is_outlier(pa_onlyFull, thresh=2.5)
        
        if isNumber(ratio): 
            ######## check ratio outlier status ######
            outliers = is_outlier(ratio_onlyFull, thresh=2.5)
            
            i = 0
            for r,o in zip(ratio_onlyFull,outliers):
                if r >1.0:
                    outliers[i] = True
                i+=1
                    
            indexRatio = ratio_onlyFull.index(ratio)
            print 'indexRatio = ',indexRatio
        
            # this is the truth value of this particular ratio value
            ratio_tval = outliers[indexRatio]
            print 'ratio_tval = ',ratio_tval
        
            # if the chosen is an outlier  - Ratio
            if ratio_tval:
                # work from the top down to find the highest ranking non-outlier
                for r,s in ratio_order:
                    # r is the ratio, s the survey
                    # here the truth value of r, which must be taken from the Full list
                    r_tval = outliers[ratio_onlyFull.index(r)]
                
                    if not r_tval:
                        ratio = r
                        ratio_survey_key = s

        if isNumber(pa):
            ######## check PA outlier status ######
            outliers = is_outlier(pa_onlyFull, thresh=2.5)
            indexPA = pa_onlyFull.index(pa)
        
            # this is the truth value of this particular PA value
            pa_tval = outliers[indexPA]
                    
            # if the chosen is an outlier - PA
            if pa_tval:
                # work from the top down to find the highest ranking non-outlier
                for p,s in pa_order:
                    # p is the PA, s the survey
                    # here the truth value of p, which must be taken from the Full list
                    p_tval = outliers[pa_onlyFull.index(p)]
                
                    if not p_tval:
                        pa = p
                        pa_survey_key = s
        
#         if len(ratio_only) >0:
#             median_ratio = np.median(ratio_onlyFull)
#             print 'median ratio: ',median_ratio
#             if abs(float(ratio) - median_ratio) >= 0.2:
#                 print '(float(ratio) - median_ratio) = ',(float(ratio) - median_ratio)
#                 for i in ratio_order:
#                     r,surv = i
#                     if abs(float(r) - median_ratio) < 0.2:
#                         ratio = r
#                         ratio_survey_key = surv
#                         break
#                         
#         if len(pa_only) >0:
#             median_pa = np.median(pa_onlyFull)
#             if abs(float(pa) - median_pa) >= 30:
#                 for i in pa_order:
#                     r,surv = i
#                     if abs(float(r) - median_pa) < 30:
#                         pa = r
#                         pa_survey_key = surv
#                         break
# 
#         
#         if ratio >= 1.0 or major == 0.0:
#             good_ratios = []
#             found_better_ratio = False
#             for r,p in zip(ratio_order,pa_order):
#                 if r[0] <= 1.0 and r[0] >= 0.0:
#                     if isNumber(p[0]):
#                         pa = p[0]
#                         pa_survey_key = p[1]
#                         
#                         ratio = r[0]
#                         ratio_survey_key = r[1]
#                         found_better_ratio = True
#                         break
#                     else:
#                         good_ratios.append(r)
#             
#             if not found_better_ratio and len(good_ratios)>0:
#                 ratio,ratio_survey_key = good_ratios[0]
            
        # apply the minor calculation
        try:
            minor = float(major)*float(ratio)
        except Exception,e:
            print
            print 'problem'
            print 'major = {0}, ratio = {1}'.format(major,ratio)
            print
            print 'options: ',options
            print 'corrected_options: ',corrected_options
            print 'd: ',d
            minor = 'x'
                
        to_return = [survey, (major,minor), ratio, pa]
    
    else:
        if found:
            if verbose:
                print 'found! ',found
            
            to_return = found
        
            # keys
            major_survey_key = to_return[0]
            ratio_survey_key = to_return[0]
            pa_survey_key = to_return[0]
        
        
        elif found_equal_axes:
            survey = found_equal_axes[0]
            major = found_equal_axes[1][0]
            minor = found_equal_axes[1][1]
            ratio = found_equal_axes[2]
            pa = found_equal_axes[3]
        
            # go through and find the next best non-unity ratio and apply it
            found_nonzero_ratio = False
            for i in ratio_order:
                i_ratio, i_survey = i
                if float(i_ratio) <= 1.0:
                    found_nonzero_ratio = True
                    new_minor = float(major) * float(i_ratio)
                    to_return = [survey,(major,new_minor),i_ratio,pa]
                
                    # keys
                    major_survey_key = survey
                    ratio_survey_key = i_survey
                    pa_survey_key = survey
                    break
                
            # otherwise, just grab the top ranking one
            if not found_nonzero_ratio:
                i_ratio, i_survey = ratio_order[0]
                new_minor = float(major) * float(i_ratio)
                to_return = [survey,(major,new_minor),i_ratio,pa]
            
                # keys
                major_survey_key = survey
                ratio_survey_key = i_survey
                pa_survey_key = survey
        
        
        elif found_no_pa:
            survey = found_no_pa[0]
            major = found_no_pa[1][0]
            minor = found_no_pa[1][1]
            ratio = found_no_pa[2]
            pa = found_no_pa[3]
    
            # grab the next best pa
            i_pa, i_survey = pa_order[0]
            to_return = [survey,(major,minor),ratio,i_pa]
        
            # keys
            major_survey_key = survey
            ratio_survey_key = survey
            pa_survey_key = i_survey
    
    #         if found_no_minor:
    #             # grab the pa from the next best survey
    #             pa = found_no_minor[3]
    #             pa_key = found_no_minor[0]
    #             
    #             # grab the rest from the one with intact (major, minor), but no pa
    #             to_return = found_no_pa
    #             
    #             # replace the missing pa
    #             to_return[3] = pa
    #             
    #             # return a key saying what came from where
    #             keys = [to_return[0],pa_key]
            
        
        elif found_no_minor:
            survey = found_no_minor[0]
            major = found_no_minor[1][0]
            minor = found_no_minor[1][1]
            ratio = found_no_minor[2]
            pa = found_no_minor[3]
    
            # go through and find the next best non-unity ratio and apply it
            found_nonzero_ratio = False
            for i in ratio_order:
                i_ratio, i_survey = i
                if float(i_ratio) <= 1.0:
                    found_nonzero_ratio = True
                    new_minor = float(major) * float(i_ratio)
                    to_return = [survey,(major,new_minor),i_ratio,pa]
                
                    # keys
                    major_survey_key = survey
                    ratio_survey_key = i_survey
                    pa_survey_key = survey
                    break
        
            # otherwise, just grab the top ranking one
            if not found_nonzero_ratio:
                i_ratio, i_survey = ratio_order[0]
                new_minor = float(major) * float(i_ratio)
                to_return = [survey,(major,new_minor),i_ratio,pa]
            
                # keys
                major_survey_key = survey
                ratio_survey_key = i_survey
                pa_survey_key = survey
    
    #         to_return = found_no_minor
    #         keys = [found_no_minor[0],found_no_minor[0]]
    
        else:
            if len(ordered_options)>0:
                # if everything else fails, just grab all the highest ranking bits
                major, major_survey_key = major_order[0]
                survey = major_survey_key
            
                ratio,ratio_survey_key = ratio_order[0]
            
                minor = float(major)*float(ratio)
            
                pa, pa_survey_key = pa_order[0]
            
                to_return = [survey, (major,minor), ratio, pa]
                        
            else:
#                 sys.stdout.write('\n No fit was found for the following: {0}\n'.format(options))
#                 sys.stdout.write('\n Trimmed options list was: {0}\n'.format(new_options))
#                 sys.stdout.write('\n Exiting....\n')
#                 sys.exit()
                to_return = ['x',('x','x'),'x','x']
                major_survey_key ='x'
                ratio_survey_key ='x'
                pa_survey_key ='x'
                

    keys = [major_survey_key, ratio_survey_key, pa_survey_key]
    
    return to_return, keys, secondary
    
    
def main():
    
#     options = [['r (SDSS Isophotal)', (5.35, 4.91), 0.92, 152.0], ['r (SDSS Petrosian)', (1.18, None), None, None], ['r (SDSS deVaucouleurs)', (0.01, 0.01), 1.0, 163.0], ['r (SDSS Exponential)', (0.0, 0.0), 0.81, 163.0]]
#     options = [['r (SDSS Exponential)', (0.0, 0.0), 0.81, 163.0], ['r (SDSS Petrosian)', (1.18, None), None, None], ['r (SDSS deVaucouleurs)', (0.01, 0.01), 1.0, 163.0]]

#     options = [['POSS1 103a-O', (270.0, 33.0), 0.12, None], ['POSS1 103a-O', (1350.0, 372.0), 0.28, None], ['RC3 D_25, R_25 (blue)', (1312.7, 434.5), 0.33, 162.0], ['RC3 A_e (Johnson B)', (476.6, None), None, None], ['RC3 D_0 (blue)', (1312.7, None), None, None], ['K_s (LGA/2MASS isophotal)', (146.4, 146.4), 1.0, 16.0], ['K_s (LGA/2MASS "total")', (250.2, 250.2), 1.0, 16.0], ['POSS1 103a-O', (1380.0, 480.0), 0.35, 162.0]]

#     options = [['POSS1 103a-O', (1350.0, 372.0), 0.28, None],\
#     ['POSS1 103a-O', (1380.0, 480.0), 0.35, None],\
#     ['POSS1 103a-O', (270.0, 33.0), 0.12, 162.0],\
#     ['RC3 D_25, R_25 (blue)', (1312.7, 434.5), 0.33, 162.0],\
#     ['RC3 A_e (Johnson B)', (476.6, None), None, None],\
#     ['RC3 D_0 (blue)', (1312.7, None), None, None],\
#     ['K_s (LGA/2MASS isophotal)', (146.4, 146.4), 1.0, 16.0],\
#     ['K_s (LGA/2MASS "total")', (250.2, 250.2), 1.0, 16.0]]
    
    
#     options = ['x']
#     options = [['POSS1 103a-O', (54.0, 96.0), 1.78, None], ['RC3 D_25, R_25 (blue)', (101.9, 60.02), 0.59, None], ['RC3 A_e (Johnson B)', (42.5, None), None, None], ['RC3 D_0 (blue)', (106.7, None), None, None]]
#     options = [['ESO-LV "Quick Blue" IIa-O', (54.95, 104.87), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (120.23, 229.45), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (142.89, 272.69), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (94.41, 180.17), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (107.15, 204.48), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (124.45, 237.5), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (71.61, 136.66), 0.52, 104.0], ['ESO-LV "Quick Blue" IIa-O', (60.26, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (28.43, 108.93), 0.26, 104.0], ['ESO-LV IIIa-F', (25.35, None), None, None]]

    options = [['ESO-Uppsala "Quick Blue" IIa-O', (96.0, 96.0), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (97.72, 97.72), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (117.49, 117.49), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (134.9, 134.9), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (74.13, 74.13), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (86.1, 86.1), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (103.51, 103.51), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (51.29, 51.29), 1.0, None], ['ESO-LV "Quick Blue" IIa-O', (51.29, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (12.15, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (5.34, 8.59), 0.62, 173.0], ['ESO-LV IIIa-F', (5.82, None), None, None], ['RC3 D_25, R_25 (blue)', (95.1, 90.82), 0.95, None], ['RC3 D_0 (blue)', (95.1, None), None, None]]
#     options = [['ESO-Uppsala "Quick Blue" IIa-O', (300.0, 60.0), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (197.24, 986.2), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (323.59, 1617.95), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (506.99, 2534.95), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (367.28, 1836.4), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (489.78, 2448.9), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (668.34, 3341.7), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (213.8, 1069.0), 0.2, 31.0], ['ESO-LV "Quick Blue" IIa-O', (49.55, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (34.62, None), None, None], ['ESO-LV "Quick Blue" IIa-O', (46.46, 227.75), 0.2, 94.0], ['ESO-LV IIIa-F', (16.03, None), None, None], ['RC3 D_25, R_25 (blue)', (198.7, 41.53), 0.21, 31.0], ['RC3 D_0 (blue)', (198.7, None), None, None]]
    best,best_keys,secondary = choose_diameter_measurement(options)
    
#     survey,(major,minor),ratio,pa = best
#     
#     if not isNumber(major):
#         if len(secondary) >0:
#             best_maj,best_keys_maj,secondary_maj = choose_diameter_measurement(secondary)
#             
#     survey_maj,(major_maj,minor_maj),ratio_maj,pa_maj = best_maj
#     
#     if not isNumber(minor):
#         if len(secondary) >0:
#             best_min,best_keys_min,secondary_min = choose_diameter_measurement(secondary)
    
    
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
    