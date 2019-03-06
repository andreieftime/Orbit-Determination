from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits


#im = fits.getdata("LA1950LIGHT00000031.REDUCED.FIT")
im = fits.getdata("LA1950LIGHT00000004.REDUCED.FIT")

plt.imshow(im, vmin = im.mean(), vmax = 2 * im.mean())
plt.gray()
plt.show()
def get_centroid(xa, ya, di):

    st = (di - 1) / 2
    table = im[ya - st:ya + st + 1, xa - st:xa + st + 1]

    #find x centroid
    saux = 0
    stot = np.sum(table)
    for i in range(di):
        saux = saux + i * np.sum(table[0:di, i])

    xc = saux / stot

    #find y centroid
    sauy = 0
    for i in range(di):
        sauy = sauy  + i * np.sum(table[i, 0:di])

    yc = sauy / stot

    #find error in x
    erx = 0;
    for i in range(di):
        erx = erx + (xc - i)**2 * np.sum(table[0:di, i])

    erx = np.sqrt(erx / (stot * (stot - 1)))

    #find error in y
    ery = 0;
    for i in range(di):
        ery = ery + (yc - i)**2 * np.sum(table[i, 0:di])

    ery = np.sqrt(ery / (stot * (stot - 1)))
    return xa + xc - (di - 1) / 2, yc + ya - (di - 1) / 2, erx, ery


#print get_centroid(464, 382, 3)
#print get_centroid(536, 364, 3)
#print get_centroid(602, 381, 5)
#print get_centroid(630, 300, 5)
#print get_centroid(726, 117, 3)
#print get_centroid(766, 99, 3)
#print get_centroid(363, 570, 7)
#print get_centroid(543, 882, 9)
#print get_centroid(586, 838, 7)
#print get_centroid(766, 600, 3)
#print get_centroid(725, 780, 7)
#print get_centroid(927, 777, 7)
##print get_centroid(623, 387, 3) coordinates for 6/30/17 31

print get_centroid(457, 391, 3)
print get_centroid(530, 372, 3)
print get_centroid(596, 389, 5)
print get_centroid(624, 309, 5)
print get_centroid(719, 126, 3)
print get_centroid(760, 108, 3)
print get_centroid(357, 578, 7)
print get_centroid(535, 891, 9)
print get_centroid(580, 846, 7)
print get_centroid(759, 609, 3)
print get_centroid(717, 788, 7)
print get_centroid(920, 785, 7)
print get_centroid(616, 396, 3) #coordinates for 6/29/17 4
