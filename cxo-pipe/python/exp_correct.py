#!/usr/bin/env python3

import os
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
import subprocess as sp

from astropy.io import fits
from ciao_contrib.runtool import *

def expcorrect(indir,obsids,frac):
    outdir = indir + 'centroid/'
    if len(obsids)==1:
        pt_src = indir + 'wide_wdect_expmap_src_check.reg'
        corr_cts_img = outdir + 'corrected_counts_{}.img'.format(frac)
        cts_img = indir + 'wide_broad_thresh.img'
        exp_map = indir + 'wide_broad_thresh.expmap'
        cts_img_no_pts = outdir + 'counts_no_src.img'

        ## Remove point sources
        dmcopy.punlearn()
        dmcopy.infile = "{}[exclude sky=region({})]".format(cts_img,pt_src)
        dmcopy.outfile = cts_img_no_pts
        dmcopy.clobber = 'yes'
        dmcopy()

        ## Load image 
        with fits.open(cts_img_no_pts) as f:
            arr = np.array(f[0].data)
            hdr = f[0].header

        ## Load exposure map
        with fits.open(exp_map) as f:
            exp_arr = np.array(f[0].data)
            exp_hdr = f[0].header

        if frac == 1.0:
            new_arr = arr
        else:
            ## Randomly select half of pixels
            xs = []
            ys = []
            indices = []
            i = 0
            for x in range(arr.shape[0]):
                for y in range(arr.shape[1]):
                    for z in range(arr[x][y]):
                        xs.append(x)
                        ys.append(y)
                        indices.append(i)
                        i += 1
            nchoices = int(len(indices)*frac)
            choices = np.random.choice(indices, size=nchoices, replace=False)
            new_arr = np.zeros_like(arr)
            for choice in choices:
                new_arr[xs[choice]][ys[choice]] += 1

        ## Exposure correct image
        norm_exp_arr = exp_arr/np.max(exp_arr)
        corr_cts = np.divide(new_arr, norm_exp_arr, out=np.zeros(new_arr.shape, dtype=float), where=norm_exp_arr!=0)
        fits.writeto(corr_cts_img, corr_cts, hdr, overwrite=True)
    else:
        pt_src = indir + 'wide_wdect_expmap_src_check.reg'
        corr_cts_img = outdir + 'corrected_counts_{}.img'.format(frac)
        for obsid in obsids:
            cts_img = indir + 'wide_{}_broad_thresh.img'.format(obsid)
            exp_map = indir + 'wide_{}_broad_thresh.expmap'.format(obsid)
            cts_img_no_pts = outdir + '{}_counts_no_src.img'.format(obsid)
            corr_cts_img_obsid = outdir + '{}_corrected_counts_{}.img'.format(obsid,frac)

            ## Remove point sources
            dmcopy.punlearn()
            dmcopy.infile = "{}[exclude sky=region({})]".format(cts_img,pt_src)
            dmcopy.outfile = cts_img_no_pts
            dmcopy.clobber = 'yes'
            dmcopy()

            ## Load image
            with fits.open(cts_img_no_pts) as f:
                arr = np.array(f[0].data)
                hdr = f[0].header

            ## Load exposure map
            with fits.open(exp_map) as f:
                exp_arr = np.array(f[0].data)
                exp_hdr = f[0].header

            if frac == 1.0:
                new_arr = arr
            else:
                ## Randomly select half of pixels
                xs = []
                ys = []
                indices = []
                i = 0
                for x in range(arr.shape[0]):
                    for y in range(arr.shape[1]):
                        for z in range(arr[x][y]):
                            xs.append(x)
                            ys.append(y)
                            indices.append(i)
                            i += 1
                nchoices = int(len(indices)*frac)
                choices = np.random.choice(indices, size=nchoices, replace=False)
                new_arr = np.zeros_like(arr)
                for choice in choices:
                    new_arr[xs[choice]][ys[choice]] += 1
    
            ## Exposure correct image
            norm_exp_arr = exp_arr/np.max(exp_arr)
            corr_cts_obsid = np.divide(new_arr, norm_exp_arr, out=np.zeros(new_arr.shape, dtype=float), where=norm_exp_arr!=0)
            fits.writeto(corr_cts_img_obsid, corr_cts_obsid, hdr, overwrite=True)

            ## Add images for all obsids
            if obsid == obsids[0]:
                corr_cts = arr*0.0

            corr_cts += corr_cts_obsid

        ## Save summed image and blanksky image
        fits.writeto(corr_cts_img, corr_cts, hdr, overwrite=True)

    return 'corrected_counts_{}.img'.format(frac)
