import numpy as np

from astropy.io import fits
from ciao_contrib.runtool import *

def expcorrect(indir,obsids):
    if len(obsids)==1:
        cts_img = indir + 'wide_broad_thresh.img'
        exp_map = indir + 'wide_broad_thresh.expmap'
        pt_src = indir + 'wide_wdect_expmap_src_check.reg'
        cts_img_no_pts = indir + 'counts_no_src.img'
        corr_cts_img = indir + 'corrected_counts.img'

        dmcopy.punlearn()
        dmcopy.infile = "{}[exclude sky=region({})]".format(cts_img,pt_src)
        dmcopy.outfile = cts_img_no_pts
        dmcopy()

        with fits.open(cts_img_no_pts) as f:
            arr = np.array(f[0].data)
            hdr = f[0].header

        with fits.open(exp_map) as f:
            bkg_arr = np.array(f[0].data)

        corr_cts_obsid = arr / (bkg_arr/np.mean(bkg_arr))

    else:
        corr_cts = []
        for obsid in obsids:
            cts_img = indir + 'wide_{}_broad_thresh.img'.format(obsid)
            exp_map = indir + 'wide_{}_broad_thresh.expmap'.format(obsid)
            pt_src = indir + 'wide_wdect_expmap_src_check.reg'
            cts_img_no_pts = indir + '{}_counts_no_src.img'.format(obsid)
            corr_cts_img_obsid = indir + '{}_corrected_counts.img'.format(obsid)

            dmcopy.punlearn()
            dmcopy.infile = "{}[exclude sky=region({})]".format(cts_img,pt_src)
            dmcopy.outfile = cts_img_no_pts
            dmcopy()

            with fits.open(cts_img_no_pts) as f:
                arr = np.array(f[0].data)
                hdr = f[0].header

            with fits.open(exp_map) as f:
                bkg_arr = np.array(f[0].data)

            corr_cts_obsid = arr / (bkg_arr/np.mean(bkg_arr))

            fits.writeto(corr_cts_img_obsid, corr_cts_obsid, hdr)


indir = '../../results/PSZ2G021.10+33.24/results/'
obsids = [6104,7940]

expcorrect(indir,obsids)

