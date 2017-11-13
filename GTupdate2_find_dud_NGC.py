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
    
    dudList = []
    count = 0
    for line in galaxyReader:
        count+=1
        sys.stdout.write("\r Count = {0}".format(count))
        sys.stdout.flush()
    
        name = line['oldName']
        completeList = eval(line['completeList'])
        
        dud = False
        if name.find('NGC') !=-1:
            if len(completeList) <2:
                dudList.append(name)
        
        if name.find("MESSIER") !=-1:
            if len(completeList) <2:
                dudList.append(name)

    print 'Here are the duds: '
    for i in dudList:
        print i        
        
    galaxyFile.close()


if __name__=="__main__":
    main()
