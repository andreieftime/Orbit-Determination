from __future__ import division
import numpy as np
from math import sin, cos, asin, acos, pi, degrees, radians

k=0.01720209894
ep = radians(23.43702)
c = 173.144633
tolerance = 1e-9
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

def norm(x):
    while x<0:
        x += 2 * pi
    while x>=2 * pi:
        x -= 2 * pi
    return x

def get_coord(r, rd):
    lv = np.linalg.norm(rd)
    lr = np.linalg.norm(r)

    a = lr/(2 - lv**2 * lr)

    h = k * np.cross(r, rd)

    e = np.sqrt(1 - np.linalg.norm(np.cross(r, rd))**2 / a)
    
    i = acos(h[2] / np.linalg.norm(h))

    O = findQuadrant(h[0] / np.linalg.norm(h) / sin(i), -h[1] / np.linalg.norm(h) / sin(i))

    f = findQuadrant(np.sqrt(a * (1 - e**2)) / e / lr * np.dot(r, rd), (a * (1 - e**2) / lr - 1) / e)

    U = findQuadrant(r[2] / lr / sin(i), (r[0] * cos(O) + r[1] * sin(O)) / lr)

    if f>pi:
        E = 2 * pi - acos((1 - lr / a) / e)
    else:
        E = acos((1 - lr / a) / e)
    M = E - e * sin(E)
    
    return a, e, degrees(i), degrees(O), degrees(norm(U - f)), degrees(M)

def give_ra(h, m, s):
    return 15 * (h + m / 60 + s / 3600)

def give_dec(d, m, s):
    if d<0:
        return d - m / 60 - s / 3600
    else:
        return d + m / 60 + s / 3600
def eq_to_ecl(x):
    x = np.matmul(np.linalg.inv(transe), x.transpose())
    return x
def give_JD(y, m, d, ut):
    return 367 * y - int(7 * (y + int( (m + 9) / 12 ) ) / 4) + int(275 * m / 9) + d + 1721013.5 + ut / 24


#f = np.loadtxt("ODfileinputecliptic.txt")
f = np.loadtxt("OD_text.txt")
#f = np.loadtxt("eftimeInput.txt")

def MoG(nr1, nr2, nr3):
    t1 = give_JD(f[nr1, 0], f[nr1, 1], f[nr1, 2], f[nr1, 3] + f[nr1, 4] / 60 + f[nr1, 5] / 3600)
    t2 = give_JD(f[nr2, 0], f[nr2, 1], f[nr2, 2], f[nr2, 3] + f[nr2, 4] / 60 + f[nr2, 5] / 3600)
    t3 = give_JD(f[nr3, 0], f[nr3, 1], f[nr3, 2], f[nr3, 3] + f[nr3, 4] / 60 + f[nr3, 5] / 3600)

    #get tau
    tau3 = k * (t3 - t2)
    tau1 = k * (t1 - t2)
    tau = k * (t3 - t1)

    #get the coordintates for each day
    ra1 = radians(give_ra(f[nr1, 6], f[nr1, 7], f[nr1, 8]))
    dec1 = radians(give_dec(f[nr1, 9], f[nr1, 10], f[nr1, 11]))
    roh1 = np.array([cos(ra1) * cos(dec1),
                     sin(ra1) * cos(dec1),
                     sin(dec1)])
    roh1 = eq_to_ecl(roh1)
    R1 = np.array([f[nr1, 12], f[nr1, 13], f[nr1, 14]])

    ra2 = radians(give_ra(f[nr2, 6], f[nr2, 7], f[nr2, 8]))
    dec2 = radians(give_dec(f[nr2, 9], f[nr2, 10], f[nr2, 11]))
    roh2 = np.array([cos(ra2) * cos(dec2),
                     sin(ra2) * cos(dec2),
                     sin(dec2)])
    roh2 = eq_to_ecl(roh2)
    R2 = np.array([f[nr2, 12], f[nr2, 13], f[nr2, 14]])

    ra3 = radians(give_ra(f[nr3, 6], f[nr3, 7], f[nr3, 8]))
    dec3 = radians(give_dec(f[nr3, 9], f[nr3, 10], f[nr3, 11]))
    roh3 = np.array([cos(ra3) * cos(dec3),
                     sin(ra3) * cos(dec3),
                     sin(dec3)])
    roh3 = eq_to_ecl(roh3)
    R3 = np.array([f[nr3, 12], f[nr3,13], f[nr3,14]])

    #get the D's
    D0 = np.dot(roh1, np.cross(roh2, roh3))

    D11 = np.dot(np.cross(R1, roh2), roh3)
    D12 = np.dot(np.cross(R2, roh2), roh3)
    D13 = np.dot(np.cross(R3, roh2), roh3)

    D21 = np.dot(np.cross(roh1, R1), roh3)
    D22 = np.dot(np.cross(roh1, R2), roh3)
    D23 = np.dot(np.cross(roh1, R3), roh3)

    D31 = np.dot(roh1, np.cross(roh2, R1))
    D32 = np.dot(roh1, np.cross(roh2, R2))
    D33 = np.dot(roh1, np.cross(roh2, R3))
    
    #solving the scalar equation of Lagrange.....pam....pam.....
    A1 = tau3 / tau
    B1 = A1 * (tau**2 - tau3**2) / 6

    A3 = -tau1 / tau
    B3 = A3 * (tau**2 - tau1**2) / 6

    A = (A1 * D21 - D22 + A3 * D23) / (-D0)
    B = (B1 * D21 + B3 * D23) / (-D0)

    E = -2 * np.dot(roh2, R2)
    F = (np.linalg.norm(R2))**2

    a = -(A**2 + A * E + F)
    b = -(2 * A * B + B * E)
    c = -B**2

    coefs = [1, 0, a, 0, 0, b, 0, 0, c]
    
    roo = np.roots(coefs)

    def get_an(r, rd):
        lv = np.linalg.norm(rd)
        lr = np.linalg.norm(r)

        a = lr/(2 - lv**2 * lr)
        n = a**(-1.5)

        return a, n

    def get_pv(f1, g1, f3, g3):#proceed position and velocity from the f&g functions
        c1 = g3 / (f1 * g3 - g1 * f3)
        c2 = -1
        c3 =  - g1 / (f1 * g3 - g1 * f3)

        ro1 = (c1 * D11 + c2 * D12 + c3 * D13) / (c1 * D0)
        ro2 = (c1 * D21 + c2 * D22 + c3 * D23) / (c2 * D0)
        ro3 = (c1 * D31 + c2 * D32 + c3 * D33) / (c3 * D0)

        r1 = ro1 * roh1 - R1
        r2 = ro2 * roh2 - R2
        r3 = ro3 * roh3 - R3

        d1 = -f3 / (f1 * g3 - f3 * g1)
        d3 = f1 / (f1 * g3 - f3 * g1)

        r2d = d1 * r1 + d3 * r3

        return r2, r2d

    def do_r2_trunc(mod):#finding the first r1, r2, r3
        f1 = 1 - tau1**2 / 2 / mod**3
        g1 = tau1 - tau1**3 / 6 /mod**3

        f3 = 1 - tau3**2 / 2 / mod**3
        g3 = tau3 - tau3**3 / 6 /mod**3

        r2, r2d = get_pv(f1, g1, f3, g3)
        ver = 0
        nr2 = 0
        nr2d = 0
        nrit = 0
        while ver==0 or np.linalg.norm(nr2 - r2)>tolerance:
            if ver!=0:
                r2 = nr2
                r2d = nr2d
            ver = 1

            ac = np.linalg.norm(r2)
            
            u = 1 / ac**3
            z = np.dot(r2, r2d) / ac**2
            q = (np.linalg.norm(r2d))**2 / ac**2 - u

            c1 = g3 / (f1 * g3 - g1 * f3)
            c2 = -1
            c3 =  - g1 / (f1 * g3 - g1 * f3)
        
            ro1 = (c1 * D11 + c2 * D12 + c3 * D13) / (c1 * D0)
            ro2 = (c1 * D21 + c2 * D22 + c3 * D23) / (c2 * D0)
            ro3 = (c1 * D31 + c2 * D32 + c3 * D33) / (c3 * D0)

            ta1 = tau1 + k * (ro2 - ro1) / c
            ta3 = tau3 + k * (ro2 - ro3) / c
        
            f1 = 1 - u * ta1**2 / 2 + u * z * ta1**3 / 2 + (3 * u * q - 15 * u * z**2 + u**2) * ta1**4 / 24
            f3 = 1 - u * ta3**2 / 2 + u * z * ta3**3 / 2 + (3 * u * q - 15 * u * z**2 + u**2) * ta3**4 / 24

            g1 = ta1 - u * ta1**3 / 6 + u * z * ta1**4 / 4
            g3 = ta3 - u * ta3**3 / 6 + u * z * ta3**4 / 4
        
            nr2, nr2d = get_pv(f1, g1, f3, g3)

            nrit += 1
            if nrit>1000:
                break

        if nrit>1000:
            return 0
        else:
            return get_coord(r2, r2d)

    def get_E(r, rd, tau, a, n):
        x = n * tau
        lr = np.linalg.norm(r)
        nx = 0
        while abs(nx - x)>tolerance:
            if nx!=0:
                x = nx
            f = x - (1 - lr / a) * sin(x) + np.dot(r, rd) * (1 - cos(x)) / n / a**2 - n * tau
            fp = 1 - (1 - lr / a) * cos(x) + np.dot(r, rd) * sin(x) / n / a**2
            nx = x - f / fp
        return x
        
    def do_r2_closed(mod):
        f1 = 1 - tau1**2 / 2 / mod**3
        g1 = tau1 - tau1**3 / 6 /mod**3

        f3 = 1 - tau3**2 / 2 / mod**3
        g3 = tau3 - tau3**3 / 6 /mod**3

        r2, r2d = get_pv(f1, g1, f3, g3)
        nr2 = 0
        nr2d = 0
        ver = 0
        nrit = 0
        while ver==0 or np.linalg.norm(r2 - nr2)>tolerance:
            if ver!=0:
                r2 = nr2
                r2d = nr2d
            ver = 1
            ac = np.linalg.norm(r2)
            
            c1 = g3 / (f1 * g3 - g1 * f3)
            c2 = -1
            c3 =  - g1 / (f1 * g3 - g1 * f3)
        
            ro1 = (c1 * D11 + c2 * D12 + c3 * D13) / (c1 * D0)
            ro2 = (c1 * D21 + c2 * D22 + c3 * D23) / (c2 * D0)
            ro3 = (c1 * D31 + c2 * D32 + c3 * D33) / (c3 * D0)

            nt1 = k * (t1 - ro1 / c)
            nt2 = k * (t2 - ro2 / c)
            nt3 = k * (t3 - ro3 / c)
            a, n = get_an(r2, r2d)

            E1 = get_E(r2, r2d, nt1 - nt2, a, n)
            E3 = get_E(r2, r2d, nt3 - nt2, a, n)

            f1 = 1 - a * (1 - cos(E1)) / ac
            f3 = 1 - a * (1 - cos(E3)) / ac

            g1 = nt1 - nt2 + (sin(E1) - E1) / n
            g3 = nt3 - nt2 + (sin(E3) - E3) / n

            nr2, nr2d = get_pv(f1, g1, f3, g3)
            
            nrit += 1
            if nrit>1000:
                break

        if nrit>1000:
            return 0
        else:
            return get_coord(r2, r2d)
    def display_trunc(x):
        a, e, i, O, mu, M = do_r2_trunc(x)
        if a>1.65 and a<1.75:
            print "Truncated version for combination ", nr2
            print "Semi-major axis, a (au): ", a
            print "Eccentricity, e: ", e
            print "Inclination w.r.t XY-plane, i (degrees): ", i
            print "Longitude of Ascending Node, OMEGA, (degrees): ", O
            print "Argument of Perifocus, w (degrees): ", mu
            print "Mean anomaly, M (degrees): ", M
            print " "

    def display_closed(x):
        a, e, i, O, mu, M = do_r2_trunc(x)
        if a>1.65 and a<1.75:
            print "Closed version for combination ", nr2
            print "Semi-major axis, a (au): ", a
            print "Eccentricity, e: ", e
            print "Inclination w.r.t XY-plane, i (degrees): ", i
            print "Longitude of Ascending Node, OMEGA, (degrees): ", O
            print "Argument of Perifocus, w (degrees): ", mu
            print "Mean anomaly, M (degrees): ", M
            print " "
    for i in roo:
        if np.isreal(i)==True and i>0:#let's hope it's only one
            i = np.real(i)
            if do_r2_trunc(i)!=0:
                display_trunc(i)
            if do_r2_closed(i)!=0:
                display_closed(i)

MoG(0, 1, 3)
MoG(0, 2, 3)
