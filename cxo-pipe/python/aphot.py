import os

import numpy as np

from astropy.io import fits
from ciao_contrib.runtool import *

def expcorrect(indir,obsids):
    outdir = indir + 'aphot/'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if len(obsids)==1:
        pt_src = indir + 'wide_wdect_expmap_src_check.reg'
        corr_cts_img = outdir + 'corrected_counts.img'
        corr_bkg_cts_img = outdir + 'corrected_bkg_counts.img'
        cts_img = indir + 'wide_broad_thresh.img'
        exp_map = indir + 'wide_broad_thresh.expmap'
        cts_img_no_pts = outdir + 'counts_no_src.img'
        bkg_img = indir + 'blank_sky_{}.evt'.format(obsids[0])

        dmcopy.punlearn()
        dmcopy.infile = "{}[exclude sky=region({})]".format(cts_img,pt_src)
        dmcopy.outfile = cts_img_no_pts
        dmcopy.clobber = 'yes'
        dmcopy()

        with fits.open(cts_img_no_pts) as f:
            arr = np.array(f[0].data)
            hdr = f[0].header

        blanksky_image.punlearn()
        blanksky_image.bkgfile = bkg_img
        blanksky_image.outroot = outdir + 'bkg_img'
        blanksky_image.imgfile = cts_img_no_pts
        blanksky_image.clobber = True
        blanksky_image()

        with fits.open(outdir + 'bkg_img_particle_bgnd.img') as f:
            bkg_arr = np.array(f[0].data)
            bkg_hdr = f[0].header

        with fits.open(exp_map) as f:
            exp_arr = np.array(f[0].data)

        corr_cts = arr / (exp_arr/np.mean(exp_arr))
        corr_bkg_cts = bkg_arr / (exp_arr/np.mean(exp_arr))

        fits.writeto(corr_cts_img, corr_cts, hdr, overwrite=True)
        fits.writeto(corr_bkg_cts_img, corr_bkg_cts, bkg_hdr, overwrite=True)
    else:
        pt_src = indir + 'wide_wdect_expmap_src_check.reg'
        corr_cts_img = outdir + 'corrected_counts.img'
        corr_bkg_cts_img = outdir + 'corrected_bkg_counts.img'
        for obsid in obsids:
            cts_img = indir + 'wide_{}_broad_thresh.img'.format(obsid)
            exp_map = indir + 'wide_{}_broad_thresh.expmap'.format(obsid)
            cts_img_no_pts = outdir + '{}_counts_no_src.img'.format(obsid)
            corr_cts_img_obsid = outdir + '{}_corrected_counts.img'.format(obsid)
            bkg_img = indir + 'blank_sky_{}.evt'.format(obsid)
            corr_bkg_img_obsid = outdir + '{}_corrected_bkg.img'.format(obsid)

            dmcopy.punlearn()
            dmcopy.infile = "{}[exclude sky=region({})]".format(cts_img,pt_src)
            dmcopy.outfile = cts_img_no_pts
            dmcopy.clobber = 'yes'
            dmcopy()

            with fits.open(cts_img_no_pts) as f:
                arr = np.array(f[0].data)
                hdr = f[0].header

            blanksky_image.punlearn()
            blanksky_image.bkgfile = bkg_img
            blanksky_image.outroot = outdir + '{}_bkg_img'.format(obsid)
            blanksky_image.imgfile = cts_img_no_pts
            blanksky_image.clobber = True
            blanksky_image()

            with fits.open(outdir + '{}_bkg_img_particle_bgnd.img'.format(obsid)) as f:
                bkg_arr = np.array(f[0].data)
                bkg_hdr = f[0].header

            with fits.open(exp_map) as f:
                exp_arr = np.array(f[0].data)

            corr_cts_obsid = arr / (exp_arr/np.mean(exp_arr))
            corr_bkg_obsid = bkg_arr / (exp_arr/np.mean(exp_arr))
            fits.writeto(corr_cts_img_obsid, corr_cts_obsid, hdr, overwrite=True)
            fits.writeto(corr_bkg_img_obsid, corr_bkg_obsid, bkg_hdr, overwrite=True)

            if obsid == obsids[0]:
                corr_cts = arr*0.0
                corr_bkg_cts = bkg_arr*0.0

            corr_cts += corr_cts_obsid
            corr_bkg_cts += corr_bkg_obsid

        fits.writeto(corr_cts_img, corr_cts, hdr, overwrite=True)
        fits.writeto(corr_bkg_cts_img, corr_bkg_cts, bkg_hdr, overwrite=True)

def bkgrate(indir):
    outdir = indir + 'aphot/'
    bkg_img = outdir + 'corrected_bkg_counts.img'
    bkg_reg_file = outdir + 'bkg.reg'

    with open(bkg_reg_file, 'r') as f:
        line = f.readlines()[0]

    center_x = line.split('(')[1].split(',')[0]
    center_y = line.split('(')[1].split(',')[1]
    radius = float(line.split(',')[2].split(')')[0])

    area = np.pi*radius**2

    dmstat.punlearn()
    dmstat.infile = '{}[sky=region({})]'.format(bkg_img,bkg_reg_file)
    dmstat.centroid = 'no'
    dmstat()

    total = float(dmstat.out_sum)
    bkgrate = total/area

    return bkgrate

#indir = '../../results/PSZ2G021.10+33.24/results/'
#obsids = [6104,7940]

indir = '../../results/PSZ2G000.13+78.04/results/'
obsids = [17159]

expcorrect(indir,obsids)
print(bkgrate(indir))

