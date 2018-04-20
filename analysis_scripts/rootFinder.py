import pylab as pb
import numpy as np
import math
import cmath 

def approximateDisk(beam=50, maxdisk=2.0, k1=0, k2=0.367, k3=0.004, 
                    ylim=[-0.3,0.3]):
    """
    Makes a plot of the error in the value returned by approximateFWHM 
    compared to computeExpectedFWHM, as a function of relative
    size of the uniform disk compared to the Gaussian FWHM beam.
    beam: in arcsec
    maxdisk: in units of beam
    k0, k1: coefficients to pass to approximateFWHM
    """
    x = []
    y = []
    beam = np.float(beam)
    maxdisk *= beam
    for disk in np.arange(0,maxdisk,maxdisk/100.):
        best = computeExpectedFWHM(beam, disk)
        approx = approximateFWHM(beam, disk, k1, k2, k3)
        x.append(disk/beam)
        y.append(100*(best-approx)/best)
    pb.clf()
    pb.plot(x,y,'k-')
    pb.xlabel('Disk diameter / Beam_FWHM')
    pb.ylabel('Percent error in approximation')
    pb.ylim(ylim)
    pb.draw()

def quadratic(a, b, c=None): 
    if c: # (ax^2 + bx + c = 0) 
        a, b = b / float(a), c / float(a) 
    t = a / 2.0 
    r = t**2 - b 
    if r >= 0: # real roots 
        y1 = math.sqrt(r) 
    else: # complex roots 
        y1 = cmath.sqrt(r) 
    y2 = -y1 
    return y1 - t, y2 - t

def cbrt(x): 
    if x >= 0: 
        return pow(x, 1.0/3.0) 
    else: 
        return -pow(abs(x), 1.0/3.0)

def polar(x, y, deg=False):
    """ 
    Convert from rectangular (x,y) to polar (r,w) 
    r = sqrt(x^2 + y^2) w = arctan(y/x) = [-\pi,\pi] = [-180,180] 
    """ 
    if deg: 
        return math.hypot(x, y), 180.0 * math.atan2(y, x) / pi 
    else: 
        return math.hypot(x, y), math.atan2(y, x)


def cubic(a, b, c, d=None): 
    if d: # (ax^3 + bx^2 + cx + d = 0) 
        a, b, c = b / float(a), c / float(a), d / float(a) 
    t = a / 3.0 
    p, q = b - 3 * t**2, c - b * t + 2 * t**3 
    u, v = quadratic(q, -(p/3.0)**3) 
    if type(u) == type(0j): # complex cubic root 
        r, w = polar(u.real, u.imag) 
        y1 = 2 * cbrt(r) * np.cos(w / 3.0) 
    else: # real root 
        y1 = cbrt(u) + cbrt(v) 
    y2, y3 = quadratic(y1, p + y1**2) 
    return y1 - t, y2 - t, y3 - t

def approximateDiskSize(beam_fwhm, fitted_size, k1=0.78, k2=0.336, 
                        k3=0.00427):
#    0 = k3*x**3 + k2*x**2 + k1*x - size
    size = fitted_size**2 - beam_fwhm**2
    a,b,c = cubic(k3,k2,k1,-size)
    if (type(a) != complex):
        if (a > 0):
            return a
    elif (type(b) != complex):
        if (b > 0):
            return b
    elif (type(c) != complex):
        if (c > 0):
            return c
    else:
        return None

def approximateFWHM(beam_fwhm, disk_diameter_arcsec, k1=0.78, k2=0.336, 
                    k3=0.00427):
    """
    Approximation to the convolution of a uniform disk with Gaussian
    beam (inspired by Baars 2007, but improved with more terms).
    It is accurate to +-0.15% over the disk/beam ratio range of 0-2.
    Baars uses k1=0, k2=log(2)/2, k3=0., which is acccurate to 7%
    over the range of 0-1.
    -Todd Hunter
    """
    return((beam_fwhm**2 + k1*disk_diameter_arcsec
                         + k2*disk_diameter_arcsec**2
                         + k3*disk_diameter_arcsec**3)**0.5)
    
