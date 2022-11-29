import os
import numpy as np

from ciao_contrib.runtool import *
from cosmocalc import cosmocalc

clus = 'PSZ2G000.13+78.04'
z = 0.171000
R500_kpc = 1150.3347423779778
res_dir = os.environ['CXO_RES_DIR'] + clus + '/results/'

cosmo = cosmocalc(z, H0=70, WM=.3)
DA = cosmo['DA_Mpc']
R500_arcsec = R500_kpc / 1000 / DA * 360/(2*np.pi) * 3600
R500 = R500_arcsec / 0.492

peak_loc_file = res_dir + 'peak_cent_pos.txt'

with open(peak_loc_file, 'r') as f:
    lines = f.readlines()
    line1 = lines[1]
    line2 = lines[2]
    peak_x = float(line1.split(' ')[0])
    peak_y = float(line2.split(' ')[0])

csmooth.punlearn()
csmooth.infile = res_dir + 'wide_broad_thresh.img'
csmooth.outfile = res_dir + 'figures/smoothed.img'
csmooth.sigmin = 2.5
csmooth.sigmax = 3.5
csmooth.sclmin = 1
csmooth.sclmax = 200
csmooth.clobber = 'yes'
csmooth()

dmcopy.punlearn()
dmcopy.infile = res_dir + 'figures/smoothed.img[x={}:{},y={}:{}]'.format(peak_x - R500/2, peak_x + R500/2, peak_y - R500/2, peak_y + R500/2)
dmcopy.outfile = res_dir + 'figures/cropped_smoothed.img'
dmcopy.clobber = 'yes'
dmcopy()

