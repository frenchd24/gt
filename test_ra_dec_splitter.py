#!/usr/bin/env python

# test
from utilities import *

# J2000RA_Dec = '((12, 30, 58.0584), (12, 16, 2.208))'
# 
# # separate coordinates
# ra,dec = eval(J2000RA_Dec)
#     
# RAh, RAm, RAs = ra[0],ra[1],ra[2]
# 
# # figure out the sign
# DE = str(dec[0])[0]
# 
# if isNumber(DE):
#     DEd = DE
#     DE = '+'
# else:
#     DEd = int(str(dec[0])[1:])
#     
# DEm = dec[1]
# DEs = dec[2]
# 
# print 'RAh, RAm, RAs'
# print RAh, RAm, RAs
# 
# print 'DE, DEd, DEm, DEs'
# print DE, DEd, DEm, DEs


##########################################################################################
########################################################################################

# J2000RA_Dec = '((12, 30, 58.0584), (12, 16, 2.208))'

# J2000RA_Dec = '((8, 59, 20.48), (+0, 51, 39.1))'

J2000RA_Dec = '((9, 5, 20.35), (-0, 2, 13.6))'

# separate coordinates
ra,dec = eval(J2000RA_Dec)

print 'dec: ',dec
    
RAh, RAm, RAs = ra[0],ra[1],ra[2]

DEd = dec[0]
DEm = dec[1]
DEs = dec[2]

# figure out the sign
# if not isNumber(DEd[0]):
#     if DEd[0] == '-':
#         DE = '-'
#     if DEd[0] == '+':
#         DE = '+'
#         
# else:
#     if DEd[0] < 0:
#         DE = '-'
#     else:
#         DE = '+'


if DEd < 0:
    DE = '-'
else:
    DE = '+'


print 'RAh, RAm, RAs'
print RAh, RAm, RAs

print 'DE, DEd, DEm, DEs'
print DE, DEd, DEm, DEs