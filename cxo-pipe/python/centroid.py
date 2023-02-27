import python.exp_correct as exp_correct
import os

import numpy as np

from ciao_contrib.runtool import *
from cosmocalc import cosmocalc
from termcolor import colored

def R500_pixels(z,R500):
    z = float(z)
    R500 = float(R500)/1000                         ## Mpc
    cosmo = cosmocalc(z, H0=70, WM=.3)
    DA = cosmo['DA_Mpc']
    R500_arcsec = R500 / DA * 360/(2*np.pi) * 3600
    R500_pixels = R500_arcsec / 0.492
    return R500_pixels

def c_smooth(img, img_dir):
    smoothed_img = 'smoothed_{}'.format(img)

    csmooth.punlearn()
    csmooth.infile = img_dir + img
    csmooth.outfile = img_dir + smoothed_img
    csmooth.sigmin = 3.0
    csmooth.sigmax = 4.0
    csmooth.sclmin = 1
    csmooth.sclmax = 200
    csmooth.clobber = 'yes'
    csmooth()

    return smoothed_img

def centroid_shift(res_dir, obsids, z, R500_kpc, frac=1.0):
    img_dir = res_dir + 'centroid/'
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)

    raw_img = 'wide_broad_thresh.img'
    peak_loc_file = res_dir + 'peak_cent_pos.txt'

    with open(peak_loc_file, 'r') as f:
        lines = f.readlines()
        line1 = lines[1]
        line2 = lines[2]
        peak_x = float(line1.split(' ')[0])
        peak_y = float(line2.split(' ')[0])

    R500 = R500_pixels(z, R500_kpc)

    apertures = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]) * R500

    downsampled_img = exp_correct.expcorrect(res_dir, obsids, frac)
    smoothed_img = c_smooth(downsampled_img, img_dir)

    centroid_xs = []
    centroid_ys = []

    i = 1
    for aperture in apertures:
        region = 'circle({},{},{})'.format(peak_x, peak_y, aperture)

        dmstat.punlearn()
        dmstat.infile = '{}[sky={}]'.format(img_dir + smoothed_img, region)
        dmstat()

        centroid = dmstat.out_cntrd_phys

        centroid_x = float(centroid.split(',')[0])
        centroid_y = float(centroid.split(',')[1])

        centroid_xs.append(centroid_x)
        centroid_ys.append(centroid_y)

        i += 1

    centroid_xs = np.array(centroid_xs)
    centroid_ys = np.array(centroid_ys)

    shifts_xs = centroid_xs - peak_x
    shifts_ys = centroid_ys - peak_y

    distances = np.sqrt(shifts_xs**2 + shifts_ys**2)

    centroid_shift = np.std(distances) / R500

    return centroid_shift

def calc_centroid_shift(clus, z, R500):
    print(colored("Calculating centroid shift...", "blue", None, ["bold"]))
    print("------------------------------------------------------------")

    base_dir = os.environ['CXO_RES_DIR']
    res_dir = base_dir + clus + '/results/'

    with open(base_dir.split('results/')[0] + 'params/param_{}.txt'.format(clus), 'r') as f:
        lines = f.readlines()

    for line in lines:
#        if line.split(' ')[0] == 'z':
#            z = float(line.split(' ')[2])
#        elif line.split(' ')[0] == 'R500':
#            R500_kpc = float(line.split(' ')[2])
#        elif line.split(' ')[0] == 'obsids':
        if line.split(' ')[0] == 'obsids':
            obsid_str = line.split(' ')[2]

    obsids = obsid_str.split('"')[1].split(',')

    frac = 1.0
    w = centroid_shift(res_dir, obsids, z, R500, frac)

    with open(res_dir + 'centroid/w.txt', 'w') as f:
        f.write('{}'.format(w))
