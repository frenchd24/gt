#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTSearch2.py, v 2 11/11/2017

A program to grab galaxy data from the galaxy table (v1.0 09/12/2013)

v1.1: updated for NewGalaxyTable5.csv , and moved to /usr/users/frenchd/gt/

v2: updated for FinalGalaxyTable11_filtered.csv

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


###########################################################################

def absoluteMag_noExtinc(m,d):
    # m is apparent magnitude, d is distance in Mpc
    M = float(m) - 5*math.log10((float(d)*10**6)/10)
    return M


def absoluteMag(m,d,e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)
    return M


def printOutInfo(line,t):
    # grab all the information in the table for an object and return it all if t == f, 
    # or just the basics if t == b
    
    Name = line['Name']
    NEDname = line['NEDname']
    z = line['z']
    RAdeg = line['RAdeg']
    DEdeg = line['DEdeg']
    RA = line['RAh'], line['RAm'], line['RAs']
    Dec = line['DE-'],line['DEd'], line['DEm'], line['DEs']
    GLON = line['GLON']
    GLAT = line['GLAT']
    Vhel = line['Vhel']
    vcorr = line['vcorr']
    distvcorr = line['distvcorr']
    RID_mean = line['RID_mean']
    RID_median = line['RID_median']
    RID_std = line['RID_std']
    RID_min = line['RID_min']
    RID_max = line['RID_max']
    bestDist = line['bestDist']
    e_bestDist = line['e_bestDist']
    MajDiam_ang = line['MajDiam_ang']
    MinDiam_ang = line['MinDiam_ang']
    MajDiam = line['MajDiam']
    MinDiam = line['MajDiam']
    e_MajDiam_ang = line['e_MajDiam_ang']
    e_MinDiam_ang = line['e_MinDiam_ang']
    e_MajDiam = line['e_MajDiam']
    e_MinDiam = line['e_MajDiam']
    inc = line['inc']
    adjustedInc = line['adjustedInc']
    e_inc = line['e_inc']
    R_vir = line['R_vir']
    PA = line['PA']
    diameterKey = line['diameterKey']
    ratioKey = line['ratioKey']
    paKey = line['paKey']
    RC3_type = line['RC3_type']
    RC3_r25= line['RC3_r25']
    RC3_pa = line['RC3_pa']
    RC3_d25 = line['RC3_d25']
    group_num = line['group_num']
    group_mem = line['group_mem']
    group_dist = line['group_dist']
    MType = line['MType']
    flag = line['flag']
    distIndicator = line['distIndicator']
    lumClass = line['lumClass']
    EBminusV = line['E(B-V)']
    Bmag = line['Bmag']
    Bmag_key = line['Bmag_key']
    Bmag_max = line['Bmag_max']
    Bmag_max_key = line['Bmag_max_key']
    Bmag_min = line['Bmag_min']
    Bmag_min_key = line['Bmag_min_key']
    Bmag_sdss = line['Bmag_sdss']
    gmag_sdss = line['gmag_sdss']
    rmag_sdss = line['rmag_sdss']
    zmag_sdss = line['zmag_sdss']
    Lstar_med = line['Lstar_med']
    e_Lstar_med = line['e_Lstar_med']
    Lstar_max = line['Lstar_max']
    e_Lstar_max = line['e_Lstar_max']
    Lstar_min = line['Lstar_min']
    e_Lstar_min = line['e_Lstar_min']
    Lstar_sdss = line['Lstar_sdss']
    e_Lstar_sdss = line['e_Lstar_sdss']
    altNames = eval(line['altNames'])

    # if t == 'b':
    print
    print 'Name: ',Name
    print 'RA, Dec: ',RAdeg, DEdeg
    print 'z: ', z
    print 'Vhel',Vhel
    print 'vcorr: ',vcorr
    print 'Best Distance: ',bestDist, ' +/- ',e_bestDist
    print 'MajDiam: ',MajDiam,' +/- ',e_MajDiam
    print 'MajDiam_ang: ',MajDiam_ang, ' +/- ',e_MajDiam_ang
    print 'DiameterKey: ',diameterKey
    print 'R_vir: ',R_vir
    print
    print 'Inc: ',inc
    print 'adjustedInc: ',adjustedInc, ' +/- ',e_inc
    print 'ratioKey: ',ratioKey
    print 'Position angle: ',PA
    print 'paKey: ',paKey
    print
    print 'Morphology: ',MType
    print 'Flag: ',flag
    print
    print 'Best L* estimates: '
    print 'Bmag: ',Bmag
    print 'Lstar_med: ',Lstar_med, ' +/- ', e_Lstar_med
    print
    print 'All estimates: '
    print 'Bmag_min: ',Bmag_min, '({0})'.format(Bmag_min_key)
    print 'Bmag_max: ',Bmag_max, '({0})'.format(Bmag_max_key)
    print 'Bmag_sdss: ',Bmag_sdss
    print 'Lstar_sdss: ',Lstar_sdss, ' +/- ',e_Lstar_sdss
    print 'Min, Max Lstar: ',Lstar_min, ', ',Lstar_max
    print
    
    if t == 'f':
        # print some additional info
        print
        print '------ Additional Info -------'
        print 'NED name: ',NEDname
        print 'RA: ',RA
        print 'Dec: ',Dec
        print 'Galactic Coordinates: ',GLON,GLAT
        print 'RID_mean: ',RID_mean
        print 'RID_median: ',RID_median
        print 'RID_min: ',RID_min
        print 'RID_max: ',RID_max
        print 'RID_std: ',RID_std
        print 'RC3_d25: ',RC3_d25
        print 'RC3_r25: ',RC3_r25
        print 'RC3_pa: ',RC3_pa
        print 'RC3_type: ',RC3_type
        print 'group number: ',group_num
        print 'group dist: ',group_dist
        print 'group members: ',group_mem
        print 'distIndicator: ',distIndicator
        print 'Luminosity Class: ',lumClass
        print 'E(B-V): ',EBminusV
        print 'gmag_sdss: ', gmag_sdss
        print 'rmag_sdss: ', rmag_sdss
        print 'zmag_sdss: ', zmag_sdss
        print
        print 'Alternative names: '
        for a in altNames:
            print '\t ',a
            
        print
        
    return
    
    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()
	if user == 'frenchd':
        filename = '/Users/frenchd/Research/gt/FinalGalaxyTable13_filtered.csv'
        filename_altNames = '/Users/frenchd/Research/gt/FinalGalaxyTable13_altNames.csv'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the file and read out the table
    theFile = open(filename,'rU')
    tableReader = csv.DictReader(theFile)
    
    
    table_dict = {}
	table_dict['Name']
	table_dict['NEDname']
	table_dict['z']
	table_dict['RAdeg']
	table_dict['DEdeg']
	table_dict['RAh']
	table_dict['RAm']
	table_dict['RAs']
	table_dict['DE-']
	table_dict['DEd']
	table_dict['DEm']
	table_dict['DEs']
	table_dict['GLON']
	table_dict['GLAT']
	table_dict['Vhel']
	table_dict['vcorr']
	table_dict['distvcorr']
	table_dict['RID_mean']
	table_dict['RID_median']
	table_dict['RID_std']
	table_dict['RID_min']
	table_dict['RID_max']
	table_dict['bestDist']
	table_dict['e_bestDist']
	table_dict['MajDiam_ang']
	table_dict['MinDiam_ang']
	table_dict['MajDiam']
	table_dict['MajDiam']
	table_dict['e_MajDiam_ang']
	table_dict['e_MinDiam_ang']
	table_dict['e_MajDiam']
	table_dict['e_MajDiam']
	table_dict['inc']
	table_dict['adjustedInc']
	table_dict['e_inc']
	table_dict['R_vir']
	table_dict['PA']
	table_dict['diameterKey']
	table_dict['ratioKey']
	table_dict['paKey']
	table_dict['RC3_type']
	table_dict['RC3_r25']
	table_dict['RC3_pa']
	table_dict['RC3_d25']
	table_dict['group_num']
	table_dict['group_mem']
	table_dict['group_dist']
	table_dict['MType']
	table_dict['flag']
	table_dict['distIndicator']
	table_dict['lumClass']
	table_dict['E(B-V)']
	table_dict['Bmag']
	table_dict['Bmag_key']
	table_dict['Bmag_max'] = []
	table_dict['Bmag_max_key'] = []
	table_dict['Bmag_min'] = []
	table_dict['Bmag_min_key'] = []
	table_dict['Bmag_sdss'] = []
	table_dict['gmag_sdss'] = []
	table_dict['rmag_sdss'] = []
	table_dict['zmag_sdss'] = []
	table_dict['Lstar_med'] = []
	table_dict['e_Lstar_med'] = []
	table_dict['Lstar_max'] = []
	table_dict['e_Lstar_max'] = []
	table_dict['Lstar_min'] = []
	table_dict['e_Lstar_min'] = []
	table_dict['Lstar_sdss'] = []
	table_dict['e_Lstar_sdss'] = []
	table_dict['altNames'] = []
		
		    
    
	for line in tableReader:
		Name = line['Name']
		NEDname = line['NEDname']
		z = line['z']
		RAdeg = line['RAdeg']
		DEdeg = line['DEdeg']
		RAh = line['RAh']
		RAm = line['RAm']
		RAs= line['RAs']
		DE = line['DE-']
		DEd = line['DEd']
		DEm = line['DEm']
		DEs = line['DEs']
		GLON = line['GLON']
		GLAT = line['GLAT']
		Vhel = line['Vhel']
		vcorr = line['vcorr']
		distvcorr = line['distvcorr']
		RID_mean = line['RID_mean']
		RID_median = line['RID_median']
		RID_std = line['RID_std']
		RID_min = line['RID_min']
		RID_max = line['RID_max']
		bestDist = line['bestDist']
		e_bestDist = line['e_bestDist']
		MajDiam_ang = line['MajDiam_ang']
		MinDiam_ang = line['MinDiam_ang']
		MajDiam = line['MajDiam']
		MinDiam = line['MajDiam']
		e_MajDiam_ang = line['e_MajDiam_ang']
		e_MinDiam_ang = line['e_MinDiam_ang']
		e_MajDiam = line['e_MajDiam']
		e_MinDiam = line['e_MajDiam']
		inc = line['inc']
		adjustedInc = line['adjustedInc']
		e_inc = line['e_inc']
		R_vir = line['R_vir']
		PA = line['PA']
		diameterKey = line['diameterKey']
		ratioKey = line['ratioKey']
		paKey = line['paKey']
		RC3_type = line['RC3_type']
		RC3_r25= line['RC3_r25']
		RC3_pa = line['RC3_pa']
		RC3_d25 = line['RC3_d25']
		group_num = line['group_num']
		group_mem = line['group_mem']
		group_dist = line['group_dist']
		MType = line['MType']
		flag = line['flag']
		distIndicator = line['distIndicator']
		lumClass = line['lumClass']
		EBminusV = line['E(B-V)']
		Bmag = line['Bmag']
		Bmag_key = line['Bmag_key']
		Bmag_max = line['Bmag_max']
		Bmag_max_key = line['Bmag_max_key']
		Bmag_min = line['Bmag_min']
		Bmag_min_key = line['Bmag_min_key']
		Bmag_sdss = line['Bmag_sdss']
		gmag_sdss = line['gmag_sdss']
		rmag_sdss = line['rmag_sdss']
		zmag_sdss = line['zmag_sdss']
		Lstar_med = line['Lstar_med']
		e_Lstar_med = line['e_Lstar_med']
		Lstar_max = line['Lstar_max']
		e_Lstar_max = line['e_Lstar_max']
		Lstar_min = line['Lstar_min']
		e_Lstar_min = line['e_Lstar_min']
		Lstar_sdss = line['Lstar_sdss']
		e_Lstar_sdss = line['e_Lstar_sdss']
		altNames = eval(line['altNames'])
	
	
		table_dict[Name] = {'Name':Name,\
		'NEDname':NEDname,\
		'z':z,\
		'RAdeg':RAdeg,\
		'DEdeg':DEdeg,\
		'RAh':RAh,\
		'RAm':RAm,\
		'RAs':RAs,\
		'DE-':DE,\
		'DEd':DEd,\
		'DEm':DEm,\
		'DEs':DEs,\
		'GLON':GLON,\
		'GLAT':GLAT,\
		'Vhel':Vhel,\
		'vcorr':vcorr,\
		'distvcorr':distvcorr,\
		'RID_mean':RID_mean,\
		'RID_median':RID_median,\
		'RID_std':RID_std,\
		'RID_min':RID_min,\
		'RID_max':RID_max,\
		'bestDist':bestDist,\
		'e_bestDist':e_bestDist,\
		'MajDiam_ang':MajDiam_ang,\
		'MinDiam_ang':MinDiam_ang,\
		'MajDiam':MajDiam,\
		'MajDiam':MajDiam,\
		'e_MajDiam_ang':e_MajDiam_ang,\
		'e_MinDiam_ang':e_MinDiam_ang,\
		'e_MajDiam':e_MajDiam,\
		'e_MajDiam':e_MajDiam,\
		'inc':inc,\
		'adjustedInc':adjustedInc,\
		'e_inc':e_inc,\
		'R_vir':R_vir,\
		'PA':PA,\
		'diameterKey':diameterKey,\
		'ratioKey':ratioKey,\
		'paKey':paKey,\
		'RC3_type':RC3_type,\
		'RC3_r25':RC3_r25,\
		'RC3_pa':RC3_pa,\
		'RC3_d25':RC3_d25,\
		'group_num':group_num,\
		'group_mem':group_mem,\
		'group_dist':group_dist,\
		'MType':MType,\
		'flag':flag,\
		'distIndicator':distIndicator,\
		'lumClass':lumClass,\
		'E(B-V)':EBminusV,\
		'Bmag':Bmag,\
		'Bmag_key':Bmag_key,\
		'Bmag_max':Bmag_max,\
		'Bmag_max_key':Bmag_max_key,\
		'Bmag_min':Bmag_min,\
		'Bmag_min_key':Bmag_min_key,\
		'Bmag_sdss':Bmag_sdss,\
		'gmag_sdss':gmag_sdss,\
		'rmag_sdss':rmag_sdss,\
		'zmag_sdss':zmag_sdss,\
		'Lstar_med':zmag_sdss,\
		'e_Lstar_med':e_Lstar_med,\
		'Lstar_max':Lstar_max,\
		'e_Lstar_max':e_Lstar_max,\
		'Lstar_min':Lstar_min,\
		'e_Lstar_min':e_Lstar_min,\
		'Lstar_sdss':Lstar_sdss,\
		'e_Lstar_sdss':e_Lstar_sdss,\
		'altNames':altNames}
		
	
 
    theFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
