import os
import numpy as np
import subprocess as sp

from ciao_contrib.runtool import *
from cosmocalc import cosmocalc
from termcolor import colored

def calc_concentration(clus, z, R500):
    print(colored("Calculating concentration...", "blue", None, ["bold"]))
    print("------------------------------------------------------------")

    base_dir = os.environ['CXO_RES_DIR']
    res_dir = base_dir + clus + '/results/'
    conc_dir = res_dir + 'concentration/'

    if not os.path.isdir(conc_dir):
        os.mkdir(conc_dir)

#    with open(base_dir.split('results/')[0] + 'params/param_{}.txt'.format(clus), 'r') as f:
#        lines = f.readlines()

#    for line in lines:
#        if line.split(' ')[0] == 'z':
#            z = float(line.split(' ')[2])
#        elif line.split(' ')[0] == 'R500':
#            R500_kpc = float(line.split(' ')[2])

    cosmo = cosmocalc(z, H0=70, WM=.3)
    DA = cosmo['DA_Mpc']
    kpc_to_arcsec = 1.0 / 1000 / DA * 360 / (2 * np.pi) * 3600
    arcsec_to_pixel = 1.0 / 0.492
    rad_40_kpc_pix = 40 * kpc_to_arcsec * arcsec_to_pixel
    rad_400_kpc_pix = 400 * kpc_to_arcsec * arcsec_to_pixel

    peak_loc_file = res_dir + 'peak_cent_pos.txt'

    with open(peak_loc_file, 'r') as f:
        lines = f.readlines()
        line1 = lines[1]
        line2 = lines[2]
        peak_x = float(line1.split(' ')[0])
        peak_y = float(line2.split(' ')[0])

    infile = res_dir + 'efile_repro_raw_clean.fits'

    dmstat.punlearn()
    dmstat.infile = '{}[(x,y)=circle({},{},{})]'.format(res_dir + 'wide_broad_flux.img', peak_x, peak_y, rad_40_kpc_pix)
    dmstat.centroid = 'no'
    dmstat()
    inner_flux = float(dmstat.out_sum)

    dmstat.punlearn()
    dmstat.infile = '{}[(x,y)=circle({},{},{})]'.format(res_dir + 'wide_broad_flux.img', peak_x, peak_y, rad_400_kpc_pix)
    dmstat.centroid = 'no'
    dmstat()
    outer_flux = float(dmstat.out_sum)

    conc = inner_flux / outer_flux

    with open(conc_dir + 'conc.txt', 'w') as f:
        f.write('{}'.format(conc))


