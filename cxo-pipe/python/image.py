import os
import numpy as np
import subprocess as sp

from ciao_contrib.runtool import *
from cosmocalc import cosmocalc
from termcolor import colored

def make_labeled_image(clus, z, R500):
    print(colored("Making cropped and labeled image...", "blue", None, ["bold"]))
    print("------------------------------------------------------------")

    base_dir = os.environ['CXO_RES_DIR']
    res_dir = base_dir + clus + '/results/'

    ## Convert R500 to pixels
    cosmo = cosmocalc(z, H0=70, WM=.3)
    DA = cosmo['DA_Mpc']
    R500_arcsec = R500 / 1000 / DA * 360/(2*np.pi) * 3600
    R500 = R500_arcsec / 0.492

    ## Get peak coordinates
    peak_loc_file = res_dir + 'peak_cent_pos.txt'
    with open(peak_loc_file, 'r') as f:
        lines = f.readlines()
        line1 = lines[1]
        line2 = lines[2]
        peak_x = float(line1.split(' ')[0])
        peak_y = float(line2.split(' ')[0])

    ## Smooth image
    csmooth.punlearn()
    csmooth.infile = res_dir + 'wide_broad_thresh.img'
    csmooth.outfile = res_dir + 'figures/counts_smoothed.img'
    csmooth.sigmin = 2.5
    csmooth.sigmax = 3.5
    csmooth.sclmin = 1
    csmooth.sclmax = 200
    csmooth.outsclfile = res_dir + 'figures/smoothed_scl.fits'
    csmooth.clobber = 'yes'
    csmooth()

    csmooth.punlearn()
    csmooth.infile = res_dir + 'wide_broad_flux.img'
    csmooth.outfile = res_dir + 'figures/smoothed.img'
    csmooth.sclmode = 'user'
    csmooth.sclmap = res_dir + 'figures/smoothed_scl.fits'
    csmooth.clobber = 'yes'
    csmooth()

    limit = 500

    ## Crop image
    dmcopy.punlearn()
    dmcopy.infile = res_dir + 'figures/smoothed.img[x={}:{},y={}:{}][bin (x,y)=::2]'.format(peak_x - limit, peak_x + limit, peak_y - limit, peak_y + limit)
    dmcopy.outfile = res_dir + 'figures/cropped_smoothed.img'
    dmcopy.clobber = 'yes'
    dmcopy()

    ## Read in centroid shift and concentration
    with open(res_dir + 'centroid/w.txt', 'r') as f:
        w = float(f.readlines()[0].strip('\n'))
    with open(res_dir + 'concentration/conc.txt', 'r') as f:
        conc = float(f.readlines()[0].strip('\n'))

    ## Make label region
    with open(res_dir + 'label.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\n')
        f.write('global color=white dashlist=8 3 width=1 font="helvetica 16 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('image\n')
        f.write('# text(140,470) text={%s}\n'%clus)
        f.write('# text(140,440) text={w=%f}\n'%w)
        f.write('# text(140,410) text={conc=%f}\n'%conc)

    sp.call(['bash', 'shell/labeled_image.sh', res_dir + 'figures/cropped_smoothed.img', res_dir + 'label.reg', res_dir + 'figures/labeled.png'])
