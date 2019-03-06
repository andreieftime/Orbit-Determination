from __future__ import division
import numpy as np
from math import sin, cos, radians, tan, asin, acos, atan, degrees

def conv_d_to_h(x):
    x = x / 15
    h = int(x)
    x -= h
    x *= 60
    m = int(x)
    x -= m
    s = x * 60
    return h, m, s

def conv_d_to_m(x):
    sign = 0
    if x<0:
        sign = -1
    x = abs(x)
    
    h = int(x)
    x -= h
    x *= 60
    m = int(x)
    x -= m
    s = x * 60

    if sign==0:
        return h, m, s
    else:
        return -h, m, s

#f = np.loadtxt("ourstars4.txt")
f = np.loadtxt("ourstars31.txt")
nstar = 12

m1 = np.array([[sum(f[:, 2])],
               [sum(f[:, 2] * f[:, 0])],
               [sum(f[:, 2] * f[:, 1])]])

m2 = np.array([[12, sum(f[:, 0]), sum(f[:, 1])],
               [sum(f[:, 0]), sum(f[:, 0] * f[:, 0]), sum(f[:, 0] * f[:, 1])],
               [sum(f[:, 1]), sum(f[:, 0] * f[:, 1]), sum(f[:, 1] * f[:, 1])]])

m2 = np.linalg.inv(m2)

m3 = np.matmul(m2, m1)

b1 = m3[0, 0]
a11 = m3[1, 0]
a12 = m3[2, 0]

m4 = np.array([[sum(f[:, 3])],
               [sum(f[:, 3] * f[:, 0])],
               [sum(f[:, 3] * f[:, 1])]])

m5 = np.array([[12, sum(f[:, 0]), sum(f[:, 1])],
               [sum(f[:, 0]), sum(f[:, 0] * f[:, 0]), sum(f[:, 0] * f[:, 1])],
               [sum(f[:, 1]), sum(f[:, 0] * f[:, 1]), sum(f[:, 1] * f[:, 1])]])

m5 = np.linalg.inv(m5)

m6 = np.matmul(m5, m4)

b2 = m6[0, 0]
a21 = m6[1, 0]
a22 = m6[2, 0]

xa = 0
xb = 0
for i in range(nstar):
    xa += (f[i, 2] - b1 - a11 * f[i, 0] - a12 * f[i, 1])**2
    xb += (f[i, 3] - b2 - a21 * f[i, 0] - a22 * f[i, 1])**2

unsra = np.sqrt(xa / (nstar - 3))
unsde = np.sqrt(xb / (nstar - 3))

#x = 623.09810671256457
#y = 386.99397590361446 coordinates for 6/30/2017 31

x = 616.09255079006778
y = 395.80361173814896 # coordinates for 6/29/2017 4

ra = b1 + a11 * x + a12 * y
de = b2 + a21 * x + a22 * y

print "b1:", b1
print "b2:", b2
print "a11:", a11
print "a12:", a12
print "a21:", a21
print "a22:", a22

print "uncertainty in dec:", unsde * 3600
print "uncertainty in ra:", unsra * 3600

print conv_d_to_h(ra), conv_d_to_m(de)

L = 162958.33

#now let's do the optional
def non_flat():
    psi = b1 + a11 * x + a12 * y + x / L
    eta = b2 + a21 * x + a22 * y + y / L

    D = radians(sum(f[:, 3]) / nstar)
    A = radians(sum(f[:, 2]) / nstar)

    delta = cos(D) - eta * sin(D)
    r = np.sqrt(D**2 + psi**2)

    ras = degrees(atan(psi / delta) + A)
    dec = atan((sin(D) + eta * cos(D)) / r)

    return ras, dec
    
