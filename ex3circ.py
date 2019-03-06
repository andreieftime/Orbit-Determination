from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits

#checks whether (x1, y1) is inside the circle (x0, y0, r) 1=inside 0=partially inside -1 = outside
def check_inside(x0, y0, x1, y1, r):
    #we'll count how many vertexes of the pixel are in the circle
    vx = [.5, -.5, -.5, .5]
    vy = [-.5, .5, -.5, .5]

    nrv = 0
    for i in range(4):
        if (x1 - x0 + vx[i])**2 + (y1 - y0 + vy[i])**2<=r**2:
            nrv += 1
    if nrv==0:
        return -1
    if nrv==4:
        return 1
    return 0

#in order to check the fractional pixels, we'll divide them in smaller squares
def get_fraction(x0, y0, x1, y1, r):
    divi = 50
    nrin = 0 #number of subpixels inside
    for i in range(divi):
        for j in range(divi):
            if (x1 - x0 + j / divi - .5)**2 + (y1 - y0 + i / divi - .5)**2<=r**2:
                nrin += 1
    return nrin / divi**2

im = fits.getdata("sampleimage.fits");

plt.imshow(im, vmin = im.mean(), vmax = 2 * im.mean());
##plt.gray();
##plt.show();

mag = input("Enter magnitude:");
xs = input("Enter x coordinate of the star:");
ys = input("Enter y coordinate of the star:");
xa = input("Enter x coordinate of the asteroid:");
ya = input("enter y coordinate of the asteroid:");

r1 = input("Enter apperture radius");
r3 = input("Enter outter annulus radius");

r2 = (int)((r1 + r3) / 2) #inner radius of the annulus

mat = np.array(im[ys - r3:ys + r3 + 1, xs - r3:xs + r3 + 1])

print mat

avs = 0
nro = 0
for i in range(2 * r3 + 1):
    for j in range(2 * r3 + 1):
       if check_inside(r3, r3, j, i, r3)==1 and check_inside(r3, r3, j, i, r2)==-1:
           nro += 1
           avs += mat[i, j]

avs = avs / nro

print avs

sig = 0
nap = 0
for i in range(2 * r3 + 1):
    for j in range(2 * r3 + 1):
        if check_inside(r3, r3, j, i, r1)==1:
            nap += 1
            sig += mat[i][j] - avs
        if check_inside(r3, r3, j, i, r1)==0:
            nap += get_fraction(r3, r3, j, i, r1)
            sig += get_fraction(r3, r3, j, i, r1) * (mat[i][j] - avs)

con = mag + 2.5 * np.log10(sig);

mata = np.array(im[ya - r1:ya + r1 + 1, xa - r1:xa + r1 + 1])

siga = 0

for i in range(2 * r1 + 1):
    for j in range(2 * r1 + 1):
        if check_inside(r1, r1, j, i, r1)==1:
            siga += mat[i][j] - avs
        if check_inside(r1, r1, j, i, r1)==0:
            siga += get_fraction(r1, r1, j, i, r1) * (mat[i][j] - avs)

print siga
print -2.5 * np.log10(siga) + con;
