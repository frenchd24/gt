#!/usr/bin/env python


# from astropy.io import fits
# from astropy import wcs
# import matplotlib.pyplot as plt
# import numpy as np
# 
# ras = np.array([184.1755, 184.1755])
# decs = np.array([69.46258, 69.46258])
# 
# w = wcs.WCS(naxis=2)
# w.wcs.ctype = ["RA---AIT", "DEC--AIT"]
# # x, y = w.wcs_world2pix(ras, decs, 0)
# 
# c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)
# x, y = w.wcs_pix2world(ras, decs, 0)
# 
# 
# #     plt.figure()
# fig = plt.figure()
# ax = fig.add_subplot(111, projection="aitoff")
# ax.scatter(x,y)
# plt.grid(True)
# plt.show()

##########################################################################################

# This script plots objects positions in a Aitoff projection using
# matplotlib. For this, random data is generated, interpreted as SkyCoord
# object and then transformed into RA and DEC radian values. If you want to
# use your own data, just import them as SkyCoord objects and continue.

# Import all required packages.
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.pylab as lab
import numpy as np

# Generate random data, for RA between 0 and 360 degrees, for DEC between
# -90 and +90 degrees.
ra_random = np.random.rand(100)*360.0
dec_random = np.random.rand(100)*180.0-90.0

# Transform the data into a SkyCoord object.
c = SkyCoord(ra=ra_random*u.degree, dec=dec_random*u.degree, frame='icrs')

# Matplotlib needs the coordinates in radians, so we have to convert them.
# Furtermore matplotlib needs the RA coordinate in the range between -pi
# and pi, not 0 and 2pi.
ra_rad = c.ra.radian
dec_rad = c.dec.radian
ra_rad[ra_rad > np.pi] -= 2. * np.pi

# Now plot the data in Aitoff projection with a grid.
fig = plt.figure()
lab.subplot(111,projection="aitoff")
# lab.subplot(111)

lab.title("Aitoff projection of our random data")
lab.grid(True)
plt.plot(ra_rad, dec_rad, 'o')
plt.show()