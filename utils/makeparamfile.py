import argparse
import astropy.constants as const
import numpy as np

from ciao_contrib.runtool import *

parser = argparse.ArgumentParser()
parser.add_argument('--file', required=True,
                    help='File containing cluster info')
args = parser.parse_args()
file = args.file

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

def listobsids(ra, dec):
    obsidoutput = find_chandra_obsid(ra,dec,instrument='acisi',detail='obsid')
    obsidlist = str(obsidoutput).split()
    obsidlist = obsidlist[2:len(obsidlist)]
    return obsidlist

with open(file, 'r') as f:
    lines = f.readlines()

for line in lines:
    name = line.split()[0]
    ra = line.split()[1]
    dec = line.split()[2]
    z = line.split()[3]
    M500 = line.split()[4]
    R500 = M500toR500(z,M500)
    obsidlist = listobsids(ra,dec)
    obsids = ''
    for obsid in obsidlist:
        obsids += '{},'.format(obsid)
    obsids = obsids.strip(',')
    print('{} has R500 {} and obsids {}'.format(name,R500,obsids))
    with open('../params/param_{}.txt'.format(name), 'w') as f:
        f.write('source_name = "{}"\n'.format(name))
        f.write('obsids = "{}"\n'.format(obsids))
        f.write('z = {}\n'.format(z))
        f.write('R500 = {}\n'.format(R500))
        f.write('use_peak = True\n')
        f.write('fixed_coord = None\n')
        f.write('fast_annuli = True\n')
        f.write('Ysz = [6e-05, 1e-05]\n')
        f.write('single_ann_spec = False\n')
        f.write('fixed_spec_annuli = False\n')
        f.write('fit_kT_profile_directly = True\n')
        f.write('file_ACCEPT = None\n')
        f.write('compute_Lcool = False\n')
        f.write('tcool_th = 7.7\n')
        f.write('input_XSZ_file = None\n')
        f.write('do_err = True\n')
        f.write('tab_obsid = obsids.split(",")\n')
        f.write('if len(tab_obsid) > 1:\n')
        f.write('    multiobs = True\n')
        f.write('else:\n')
        f.write('    multiobs = False\n')
