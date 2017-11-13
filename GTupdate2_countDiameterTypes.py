#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: GTupdate2_countDiameterTypes.py  v 1.0  03/17/15

Count how many of each diameter types there are in the table

"""

import getpass
import sys
import os
import csv
import string
from math import *

# programName = 'convertTable2'
# version = '2.0'
# date = '04/08/14'
    
################################################################

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def trunc(x,n):
    # truncates a number, x, to n places without rounding
    if isNumber(x):
        x = float(x)
        slen = len('%.*f' %(n,x))
        return str(x)[:slen]
    else:
        return x
    

def main():
    
    # Decide which directory to use
    user = getpass.getuser()
    if user == 'frenchd':
#         galaxyFileName = '/usr/users/frenchd/gt/NewGalaxyTable3.csv'
        galaxyFileName = '/usr/data/moosejaw/frenchd/GT_update2/returnPA11111111111111.csv'
    
    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()    
        
    galaxyFile = open(galaxyFileName,'rU')
    galaxyReader = csv.DictReader(galaxyFile)
    
    d = {}
    for line in galaxyReader:
#         dKey = line['diameterKey']
        dList = eval(line['completeList'])
        for i in dList:
            dKey = i[0]
        
            if d.has_key(dKey):
                i = d[dKey]
                i +=1
                d[dKey] = i
                    
            # otherwise make a new dictionary entry for this
            else:
                d[dKey] = 1


    keys = d.keys()
    values = d.values()
    
    results = zip(values,keys)
    results.sort(reverse=True)
    
    # unzip
    numbers,dTypes = zip(*results)

    for s,n in zip(dTypes,numbers):
        print s, ' : ',n
    
    galaxyFile.close()


if __name__=="__main__":
    main()
