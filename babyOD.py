from __future__ import division
import numpy as np
from math import sin, cos, asin, acos, degrees, radians, pi

def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

def get_coord(r, rd):
    k=0.01720209894
    lv = np.linalg.norm(rd)
    lr = np.linalg.norm(r)

    a = lr/(2 - lv**2 * lr)

    h = k * np.cross(r, rd)

    e = np.sqrt(1 - np.linalg.norm(np.cross(r, rd))**2 / a)

    i = acos(h[2] / np.linalg.norm(h))

    O = findQuadrant(h[0] / np.linalg.norm(h) / sin(i), -h[1] / np.linalg.norm(h) / sin(i))

    f = asin(np.sqrt(a * (1 - e**2)) / e / lr * np.dot(r, rd))

    U = findQuadrant(r[2] / lr / sin(i), (r[0] * cos(O) + r[1] * sin(O)) / lr)

    if f<0:
        E = 2 * pi - acos((1 - lr / a) / e)
    else:
        E = acos((1 - lr / a) / e)
    M = E - e * sin(E)
    
    return a, e, degrees(i), degrees(O), degrees(U - f), degrees(M)

r1 = np.array([0.244,2.17,-.455])
r1d = np.array([-.731,-.0041,.0502])

print get_coord(r1, r1d)
