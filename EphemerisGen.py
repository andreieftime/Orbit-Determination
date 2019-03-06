from __future__ import division
import numpy as np
from math import radians, pi, sin, cos, asin, acos, degrees

#a = 1.704987880150487
#e = 0.4518788622141914
#M = radians(287.338829192644)
#O = radians(144.567249166763)
#i = radians(4.620927383463039)
#w = radians(197.6954745176973)
e = 0.09618818215877029
a = 1.844195082375194
i = radians(23.66084024789276)
O = radians(132.0979894264106)
w = radians(129.1630964507269)
M = radians(113.1648904822557)
ep = radians(23.43702)
k = 0.0172020989
day = 118
n = k * (a**(-1.5))

#etos = np.array([[-.207406455],
#                 [.9953010074],
#                 [-3.31375772 * 10**(-5)]])

etos = np.array([[.1579922623],
                 [-.9705438166],
                 [3.780898978 * 10**(-5)]])

#normalize an angle in order to be between 0 and 2 * pi
def norm(x):
    while x>2 * pi:
        x -= 2 * pi
    return x

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
    semn = 0
    if x<0:
        semn = 1
        x = -x
    
    h = int(x)
    x -= h
    x *= 60
    m = int(x)
    x -= m
    s = x * 60

    if semn==1:
        return -h, m, s
    
    return h, m, s

#trans matrixes
transO = np.array([[cos(O), -sin(O), 0],
                   [sin(O), cos(O), 0],
                   [0, 0, 1]])
transi = np.array([[1, 0, 0],
                   [0, cos(i), -sin(i)],
                   [0, sin(i), cos(i)]])
transw = np.array([[cos(w), -sin(w), 0],
                   [sin(w), cos(w), 0],
                   [0, 0, 1]])
transe = np.array([[1, 0, 0],
                   [0, cos(ep), -sin(ep)],
                   [0, sin(ep), cos(ep)]])

def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

#we solve the equation E - e * sin(E) = M for a certain M
def solvekep(x):
    Eguess = x
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = x
    while abs(Mguess - Mtrue) > .0001:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Eguess - e*sin(Eguess) - Mtrue) / (1 - e*cos(Eguess))
    return Eguess


#we update the coordinates depending on time
def update(t):
    Mtrue = n*t + M
    Etrue = solvekep(Mtrue)
    print Etrue
    
    coord1 = np.array([[a * cos(Etrue) - a * e],
                      [a * np.sqrt(1 - e**2) * sin(Etrue)],
                      [0]])

    print coord1
    #convert the coordinates to eclipitc
    m0 = np.matmul(transw, coord1)
    m1 = np.matmul(transi, m0)
    m2 = np.matmul(transO, m1)
    m2 = m2 + etos
    m3 = np.matmul(transe, m2)
    print m3

    return m3

eq_coord = update(day)

r = np.sqrt(eq_coord[0,0]**2 + eq_coord[1, 0]**2 + eq_coord[2, 0]**2)
dec = asin(eq_coord[2, 0] / r)
print conv_d_to_m(degrees(dec))
ra = findQuadrant(eq_coord[1, 0] / r / cos(dec), eq_coord[0, 0] / r / cos(dec))
print conv_d_to_h(degrees(ra))

