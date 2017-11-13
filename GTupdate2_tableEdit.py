#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_tableEdit.py, v1 09/12/17

Read in the fixed width galaxy table and left justify it

'''

import sys
import os
import csv
# import string
import warnings
import numpy
# import atpy
import getpass
from utilities import *
import math
from astropy.io import ascii
from astropy.table import Table
import pandas as pd

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################

    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':        
        inFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group_ascii2_1.dat'

        outFilename = '/usr/data/moosejaw/frenchd/GT_update2/rejected_final_combined_group_ascii2_1_edit.dat'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
        
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
        
    # open the file
#     inFile = ascii.read(inFilename, format='fixed_width', names = fieldnames)
#     inFile = Table.read(inFilename,format='ascii.fixed_width')
    
    inFile = pd.read_fwf(inFilename)
    
    print inFile.head()
    
#     print
#     print 'colnames: ',inFile.colnames
#     print
#     print "inFile['NEDname'][2] = ",inFile['NEDname'][2]

    nullFloat = -99.99
    nullInt = -99
    nullStr = 'x'


        
    
    d = {'Name':'U29',\
    'NEDname':'U29',\
    'z':'f8',\
    'RAdeg':'f8',\
    'DEdeg':'f8',\
    'RAh':'i4',\
    'RAm':'i4',\
    'RAs':'f8',\
    'DE-':'U1',\
    'DEd':'i4',\
    'DEm':'i4',\
    'DEs':'f8',\
    'GLON':'f8',\
    'GLAT':'f8',\
    'Vhel':'i4',\
    'vcorr':'i4',\
    'distvcorr':'f4',\
    'RID_mean':'f4',\
    'RID_std':'f4',\
    'RID_min':'f4',\
    'RID_max':'f4',\
    'bestDist':'f4',\
    'e_bestDist':'f4',\
    'MajDiam_ang':'f4',\
    'MinDiam_ang':'f4',\
    'e_MajDiam_ang':'f4',\
    'e_MinDiam_ang':'f4',\
    'MajDiam':'f4',\
    'MinDiam':'f4',\
    'e_MajDiam':'f4',\
    'e_MinDiam':'f4',\
    'inc':'i4',\
    'adjustedInc':'i4',\
    'e_inc':'i4',\
    'PA':'i4',\
    'diameterKey':'U14',\
    'ratioKey':'U14',\
    'paKey':'U14',\
    'RC3_type':'U8',\
    'RC3_d25':'f4',\
    'RC3_r25':'f4',\
    'RC3_pa':'i4',\
    'group_num':'i4',\
    'group_mem':'i4',\
    'group_dist':'f4',\
    'MType':'U20',\
    'distIndicator':'U20',\
    'lumClass':'U10',\
    'E(B-V)':'f4',\
    'Bmag':'f4',\
    'Bmag_key':'U20',\
    'Bmag_max':'f4',\
    'Bmag_max_key':'U20',\
    'Bmag_min':'f4',\
    'Bmag_min_key':'U20',\
    'Bmag_sdss':'f4',\
    'gmag_sdss':'f4',\
    'rmag_sdss':'f4',
    'zmag_sdss':'f4',\
    'Lstar_med':'f4',\
    'e_Lstar_med':'f4',\
    'Lstar_max':'f4',\
    'e_Lstar_max':'f4',\
    'Lstar_min':'f4',\
    'e_Lstar_min':'f4',\
    'Lstar_sdss':'f4',\
    'e_Lstar_sdss':'f4',\
    'altNames':'U110'}

    
    # initiate the table
 #    t = table.Table(names=fieldnames,dtype=(d['Name'],\
#     d['NEDname'],\
#     d['z'],\
#     d['RAdeg'],\
#     d['DEdeg'],\
#     d['RAh'],\
#     d['RAm'],\
#     d['RAs'],\
#     d['DE-'],\
#     d['DEd'],\
#     d['DEm'],\
#     d['DEs'],\
#     d['GLON'],\
#     d['GLAT'],\
#     d['HRV'],\
#     d['vcorr'],\
#     d['distvcorr'],\
#     d['RID_mean'],\
#     d['RID_std'],\
#     d['RID_min'],\
#     d['RID_max'],\
#     d['bestDist'],\
#     d['e_bestDist'],\
#     d['MajDiam_ang'],\
#     d['MinDiam_ang'],\
#     d['e_MajDiam_ang'],\
#     d['e_MinDiam_ang'],\
#     d['MajDiam'],\
#     d['MinDiam'],\
#     d['e_MajDiam'],\
#     d['e_MinDiam'],\
#     d['inc'],\
#     d['adjustedInc'],\
#     d['e_inc'],\
#     d['PA'],\
#     d['diameterKey'],\
#     d['ratioKey'],\
#     d['paKey'],\
#     d['RC3_type'],\
#     d['RC3_d25'],\
#     d['RC3_r25'],\
#     d['RC3_pa'],\
#     d['group_num'],\
#     d['group_mem'],\
#     d['group_dist'],\
#     d['MType'],\
#     d['distIndicator'],\
#     d['lumClass'],\
#     d['E(B-V)'],\
#     d['Bmag'],\
#     d['Bmag_key'],\
#     d['Bmag_max'],\
#     d['Bmag_max_key'],\
#     d['Bmag_min'],\
#     d['Bmag_min_key'],\
#     d['Bmag_sdss'],\
#     d['gmag_sdss'],\
#     d['rmag_sdss'],
#     d['zmag_sdss'],\
#     d['Lstar_med'],\
#     d['e_Lstar_med'],\
#     d['Lstar_max'],\
#     d['e_Lstar_max'],\
#     d['Lstar_min'],\
#     d['e_Lstar_min'],\
#     d['Lstar_sdss'],\
#     d['e_Lstar_sdss'],\
#     d['altNames']))

  

##########################################################################################
##########################################################################################
        # now write it to the new file
        
#         if not count % 1:
#             ascii.write(t, outFilename,format='fixed_width',overwrite=True,delimiter=' ',\
#             bookend=False,delimiter_pad=None)
            

    
    # close the files
#     inFile.close()

    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
