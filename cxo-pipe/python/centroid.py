import exp_correct
import os

import numpy as np

from ciao_contrib.runtool import *
from cosmocalc import cosmocalc

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
    csmooth.sigmin = 2.5
    csmooth.sigmax = 3.5
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

    print('R500 in pixels is ', R500)

    apertures = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]) * R500

    print('Downsampling and exposure correcting...')
    downsampled_img = exp_correct.expcorrect(res_dir, obsids, frac)
    print('Smoothing...')
    smoothed_img = c_smooth(downsampled_img, img_dir)

    centroid_xs = []
    centroid_ys = []

    i = 1
    for aperture in apertures:
        print('Doing aperture {}...'.format(i))
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

## 0
#res_dir = '../results/0159+0030/results/'
#obsids = [5777]
#z = 0.3856
#R500_kpc = 837.9526307043083

## 0.10
res_dir = '../results/1222+2709/results/'
obsids = [5766]
z = 0.4720
R500_kpc = 762.2990624071908

## 0.26
#res_dir = '../results/0333−2456/results/'
#obsids = [5764]
#z = 0.4751
#R500_kpc = 737.5635994144291

## 0.52
#res_dir = '../results/0542−4100/results/'
#obsids = [914]
#z = 0.6420
#R500_kpc = 889.8921148221584

## 1.37
#res_dir = '../results/0230+1836/results/'
#obsids = [5754]
#z = 0.8106
#R500_kpc = 788.2650382921674

## 2.90
#res_dir = '../results/0152−1358/results/'
#obsids = [913,21703,22856]
#z = 0.8325
#R500_kpc = 813.9091318386161

#res_dir = '../results/PSZ2G000.13+78.04/results/'
#obsids = [17159]
#z = 0.171000
#R500_kpc = 1150.3347423779778

#res_dir = '../results/PSZ2G021.10+33.24/results/'
#obsids = [6104,7940]
#z = 0.151400
#R500_kpc = 1331.756149120992

ws = []

frac = 1.0
w = centroid_shift(res_dir, obsids, z, R500_kpc, frac)
ws.append(w)
print('centroid shift downsampled to {} is {}'.format(frac,w))

frac = 0.75
w = centroid_shift(res_dir, obsids, z, R500_kpc, frac)
ws.append(w)
print('centroid shift downsampled to {} is {}'.format(frac,w))

frac = 0.5
w = centroid_shift(res_dir, obsids, z, R500_kpc, frac)
ws.append(w)
print('centroid shift downsampled to {} is {}'.format(frac,w))

frac = 0.25
w = centroid_shift(res_dir, obsids, z, R500_kpc, frac)
ws.append(w)
print('centroid shift downsampled to {} is {}'.format(frac,w))

frac = 0.1
w = centroid_shift(res_dir, obsids, z, R500_kpc, frac)
ws.append(w)
print('centroid shift downsampled to {} is {}'.format(frac,w))

frac = 0.05
w = centroid_shift(res_dir, obsids, z, R500_kpc, frac)
ws.append(w)
print('centroid shift downsampled to {} is {}'.format(frac,w))

frac = 0.02
w = centroid_shift(res_dir, obsids, z, R500_kpc, frac)
ws.append(w)
print('centroid shift downsampled to {} is {}'.format(frac,w))

print(ws)
