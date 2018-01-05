#!/usr/bin/env python3.5
# my_orbdyn_proc.py
# 
# Author: Kevin Graves
# Last Modified: Jan 5, 2018
# Description: Module for my commonly used functions for orbital dynamics
# 
#
###############################################################################

import numpy as np
import sys


def el2xv(a, ecc, inc, capom, omega, capm, mu):
    """
    Function to change orbital elements semimajor axis (AU), eccentricity, 
    inclination (rad), longitude of ascending node - capom (rad),
    argument of pericenter - omega (rad), and the mean anomaly - capm (rad)
    to position (AU) and velocity (AU/day) using mu (AU^3/day^2 system)
    Only tested for the elliptical case
    """
    # acceptable error
    eps = 1e-12
    # Find Eccentric Anomaly
    cape = capm
    n = 0
    while True: 
        capep1 = cape - (cape - ecc*np.sin(cape) - capm)/(1. - ecc*np.cos(cape))
        if n > 10000:
            print(n,cape,capep1)
            print(abs(capep1 - cape) , eps,abs(capep1 - cape) < eps )
            print("No Convergence in finding eccentric anamoly")
            sys.exit(1)
        if ( abs(capep1 - cape) < eps):
            cape = capep1
            break
        cape = capep1
        n += 1


    r = np.zeros(3)
    v = np.zeros(3)
    nu = 2.*np.arctan(np.sqrt((1. + ecc)/(1. - ecc))*np.tan(cape/2.))
    rn = a*(1. - ecc*np.cos(cape))
    h = np.sqrt(mu*a*(1. - ecc*ecc))
    r[0] = rn*(np.cos(capom)*np.cos(omega+nu) - np.sin(capom)*np.sin(omega+nu)*np.cos(inc))
    r[1] = rn*(np.sin(capom)*np.cos(omega+nu) + np.cos(capom)*np.sin(omega+nu)*np.cos(inc))
    r[2] = rn*(np.sin(inc)*np.sin(omega+nu))
    v[0] = r[0]*h*ecc/(rn*a*(1.-ecc*ecc))*np.sin(nu) - h/rn*(np.cos(capom)*np.sin(omega+nu) + np.sin(capom)*np.cos(omega+nu)*np.cos(inc))
    v[1] = r[1]*h*ecc/(rn*a*(1.-ecc*ecc))*np.sin(nu) - h/rn*(np.sin(capom)*np.sin(omega+nu) - np.cos(capom)*np.cos(omega+nu)*np.cos(inc))
    v[2] = r[2]*h*ecc/(rn*a*(1.-ecc*ecc))*np.sin(nu) + h/rn*np.sin(inc)*np.cos(omega+nu)

    return r, v


def xv2el(x,v,mu):
    """
    Function to change the state vectors (position and velocity)
    to the orbital elements [semimajor axis, eccentricity, 
    inclination, longitude of ascending node (capom),
    argument of pericenter (omega), and the mean anomaly (capm)]
    Only tested for the elliptical case
    """

    # Magnitudes
    xmag = np.sqrt(np.dot(x,x))
    vmag = np.sqrt(np.dot(v,v))

    # Calculate the angular momentum h
    h = np.cross(x,v)
    hmag = np.sqrt(np.dot(h,h))

    # Find eccentricity
    ecc_v = np.cross(v,h)/mu - x/np.sqrt(np.dot(x,x))
    ecc = np.sqrt(np.dot(ecc_v,ecc_v))

    # Find vector pointing towards ascending node
    n = np.cross(np.array([0,0,1]),h)
    nmag = np.sqrt(np.dot(n,n))

    # Energy of system
    energy = vmag**2./2.-mu/xmag

    if abs(ecc-1.0) > 10.**(-12.):
       a = -mu/(2.*energy)
       p = a*(1.-ecc**2.)
    else:
       p = hmag^2/mu
       a = np.nan

    inc = np.arccos(h[2]/hmag)

    capom = np.arccos(n[0]/nmag)

    if n[1] < 0.:
       capom = 2.*np.pi - capom

    # Argument of periapsis
    omega = np.arccos(np.dot(n,ecc_v)/(nmag*ecc))

    if ecc_v[2] < 0:
       omega = 2.*np.pi - omega

    # True Anomaly
    nu = np.arccos(np.dot(ecc_v,x)/(ecc*xmag))

    if np.dot(x,v) < 0.:
       nu = 2.*np.pi - nu


    # Eccentric Anomaly
    cape = 2.*np.arctan(np.tan(nu/2.)/np.sqrt((1.+ecc)/(1.-ecc)))

    # Mean Anamoly
    capm = cape - ecc*np.sin(cape)

    # Semi-major Axis
    a = 1./(2./xmag - vmag**2./mu)


    return a, ecc, inc, capom, omega, capm

# Orbital Similarity Functions
def orb_dist_Dd(el1,el2, mu):
    """
    Function to calculate the Drummond (1981) Orbital Similarity
    distance (Dd) given 6 orbital elements of two bodies:
    Semimajor axis (AU), Eccentricity, Inclination (rad), 
    Longitude of Ascending Node (rad), & Argument of periapsis (rad),
    & Mean Anomaly (rad)
    (transferred in tuples) and the mass of the sun. 
    This equation was taken from Rozek et al. (2010)
    """
    a1, ecc1, inc1, capom1, omega1, capm1 = el1
    a2, ecc2, inc2, capom2, omega2, capm2 = el2
    term1 = ( (ecc2 - ecc1)/(ecc2 + ecc1) )**2.

    q1 = a1*(1. - ecc1)
    q2 = a2*(1. - ecc2)
    term2 = ( (q2 - q1)/(q2 + q1) )**2.
        
    r1, v1 = el2xv(a1, ecc1, inc1, capom1, omega1, capm1, mu)
    r2, v2 = el2xv(a2, ecc2, inc2, capom2, omega2, capm2, mu)
  
    h1 = np.cross(r1,v1)  
    h2 = np.cross(r2,v2)  
    h1_mag = np.sqrt(np.dot(h1,h1))
    h2_mag = np.sqrt(np.dot(h2,h2))

    try:
        I21   = np.arccos(np.dot(h1,h2)/( h1_mag*h2_mag ) )
    except FloatingPointError:
        print("FloatingPointError while trying to calculate I21")
        print(el1)
        print(el2)
        sys.exit(1)
    term3 = (I21/np.pi)*(I21/np.pi)

    vxh1 = np.cross(v1,h1)
    vxh2 = np.cross(v2,h2)

    r1_mag = np.sqrt(np.dot(r1,r1))
    r2_mag = np.sqrt(np.dot(r2,r2))
    e1vec = vxh1/mu - r1/r1_mag
    e2vec = vxh2/mu - r2/r2_mag
    th21  = np.arccos(np.dot(e1vec,e2vec)/ (ecc1*ecc2) )
    term4 = ((ecc2 + ecc1)/2.)**2. * (th21/np.pi)*(th21/np.pi)

    Dd2 = term1 + term2 + term3 + term4
    return np.sqrt(Dd2)

def v6Res(a):
    """
    Function to calculate the inclination of the nu6 resonance for a
    range of semimajor axis values (a) and eccentricities (e)
    NOT FINISHED
    """
    # Solar System Constants
    J2 = 2e-7
    J4 = -4e9
    RSun = 6.957e5
    Msun = 1.9885e30
    MPlan = np.array([0.3301e24,4.8676e24,(5.9726e24+0.07342e24),0.64174e24,
                      1898.3e24,568.36e24,86.816e24,102.42e24])
    aPlan = np.array([0.387,0.723,1.00000011,1.524,
                      5.204,9.582,19.201,30.047])
    G = 6.672e-11

    n = np.sqrt(G*Msun/a**3.)
    def alpha(atp,apl):
        if atp < apl:
            return atp/apl
        else:
            return apl/atp
    def alphabar(atp,apl):
        if atp < apl:
            return 1.
        else:
            return apl/atp

    def b32to1(alph):
        "Eq. 6.67 in Murray & Dermott"
        term1 = (1.5*(1.5+1.))/(1.)*alph
        term2 = 1. + (1.5*(1.5+1.))/(2.)*alph**2.
        term3 = (1.5*(1.5+1.)*(1.5+1.)*(1.5+1.+1.))/(2.*(1.+1.)*(1.+2.))*alph**4.
        return 2.*(term1 * (1. + term2 + term3))



















