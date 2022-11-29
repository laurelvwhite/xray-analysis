import os
import subprocess as sp

from ciao_contrib.runtool import *
from cosmocalc import cosmocalc

clus = 'PSZ2G000.13+78.04'
z = 0.171000
R500_kpc = 1150.3347423779778

res_dir = os.environ['CXO_RES_DIR'] + clus + '/results/'
conc_dir = res_dir + 'concentration/'

if not os.path.isdir(conc_dir):
    os.mkdir(conc_dir)

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

infile = res_dir + 'efile_repro_raw_clean.fits'
sp.call(['bash', 'shell/make_images.sh', infile, conc_dir, '0.5:5:1'])

## CHANGE TO 40 AND 400 KPC
dmstat.punlearn()
dmstat.infile = '{}[(x,y)=circle({},{},{})]'.format(conc_dir + , peak_x, peak_y, R500)
dmstat.centroid = 'no'
dmstat()

flux = float(dmstat.out_sum)

