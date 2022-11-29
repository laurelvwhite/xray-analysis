#!/usr/bin/python

import argparse
import astropy.constants as const
import astropy.units as u
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--z',
                    help='redshift')
parser.add_argument('--M500',
                    help='M500 in units of 10**14 Msun')
args = parser.parse_args()

z = args.z
M500 = args.M500

def M500toR500(z,m500):
    z = float(z)
    m500 = float(m500)
    ## Calculate critical density at cluster redshift
    H_0 = 70*3.24078*10**(-20)                                                   ## [s^-1]
    Hubble = H_0*np.sqrt(0.3*(1+z)**3+0.7)                                       ## [s^-1]
    pcrit = 3*Hubble**2/(8*np.pi*const.G.value)                                  ## [kg*m^-3]
    ## Convert m500 to r500
    r500 = (m500*10**14*const.M_sun.value/(500*pcrit)*(3/(4*np.pi)))**(1.0/3.0)  ## [m]
    r500 *= 3.241*10**-20                                                        ## [kpc]
    return r500

def M500toR5002(z,m500):
    z = float(z)
    m500 = float(m500)
    H_0 = 70 * u.km / u.s / u.Mpc
    Hubble = H_0*np.sqrt(0.3*(1+z)**3+0.7)
    pcrit = 3*Hubble**2/(8*np.pi*const.G)
    r500 = (m500*10**14*const.M_sun/(500*pcrit)*(3/(4*np.pi)))**(1.0/3.0)
    return r500.to(u.kpc)

print('R500 in kpc is {}'.format(M500toR500(z,M500)))
