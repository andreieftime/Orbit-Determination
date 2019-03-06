from __future__ import division
from visual import *
import numpy as np

a = 2.773017979589484
e = 0.1750074901308245
M = radians(336.0050001501443)
O = radians(108.032597191534)
i = radians(16.34548466739393)
w = radians(74.95130563682554)
k = 0.0172020989

#we define the scale we want for our orbits with respect to the Sun
scale = 150

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

#we solve the equation E - e * sin(E) = M for a certain M
def solvekep(M):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while abs(Mguess - Mtrue) > 1e-004:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Eguess - e*sin(Eguess) - Mtrue) / (1 - e*cos(Eguess))
    return Eguess

n = k * a**(-1.5)

#we update the coordinates depending on time
def update(t):
    Mtrue = n*t + M
    Etrue = solvekep(Mtrue)

    coord1 = np.array([[a * cos(Etrue) - a * e],
                      [a * np.sqrt(1 - e**2) * sin(Etrue)],
                      [0]])
    
    #convert the coordinates to eclipitc
    m1 = np.matmul(transw, coord1)
    m2 = np.matmul(transi, m1)
    m3 = np.matmul(transO, m2)

    return m3

time = 0
dt = 5

asteroid = sphere(pos=(0, 0, 0), radius=(15), color=color.white)
asteroid.trail = curve(color=color.white)
sun = sphere(pos=(0,0,0), radius=(50), color=color.yellow)

while True:
    rate(200)

    time += dt
    
    asteroid.pos = update(time)*scale
    asteroid.trail.append(pos=asteroid.pos)  
