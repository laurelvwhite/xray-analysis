import os
import time

import matplotlib.pyplot as plt
import numpy as np
import subprocess as sp

from astropy.io import fits
from ciao_contrib.runtool import *
from cosmocalc import cosmocalc

def getR500inpixels(z,R500):
    z = float(z)
    R500 = float(R500)/1000               ## Mpc
    cosmo = cosmocalc(z, H0=70, WM=.3)
    DA = cosmo['DA_Mpc']
    R500_arcsec = R500 / DA * 360/(2*np.pi) * 3600
    R500_pixels = R500_arcsec / 0.492
    return R500_pixels

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

        ## Convert blanksky background into an image
        blanksky_image.punlearn()
        blanksky_image.bkgfile = bkg_img
        blanksky_image.outroot = outdir + 'bkg_img'
        blanksky_image.imgfile = cts_img_no_pts
        blanksky_image.clobber = True
        blanksky_image()

        ## Load blanksky image
        with fits.open(outdir + 'bkg_img_particle_bgnd.img') as f:
            bkg_arr = np.array(f[0].data)
            bkg_hdr = f[0].header

        ## Load exposure map
        with fits.open(exp_map) as f:
            exp_arr = np.array(f[0].data)

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
        nchoices = int(len(indices)/2)
        choices = np.random.choice(indices, size=nchoices, replace=False)
        new_arr = np.zeros_like(arr)
        for choice in choices:
            new_arr[xs[choice]][ys[choice]] += 1

        ## Repeat for background
        bkg_arr *= 1000
        xs = []
        ys = []
        indices = []
        i = 0
        for x in range(bkg_arr.shape[0]):
            for y in range(bkg_arr.shape[1]):
                for z in range(int(bkg_arr[x][y])):
                    xs.append(x)
                    ys.append(y)
                    indices.append(i)
                    i += 1
        nchoices = int(len(indices)/2)
        choices = np.random.choice(indices, size=nchoices, replace=False)
        new_bkg_arr = np.zeros_like(bkg_arr)
        for choice in choices:
            new_bkg_arr[xs[choice]][ys[choice]] += 1
        new_bkg_arr /= 1000

        ## Exposure correct image and blanksky image
        corr_cts = new_arr / (exp_arr/np.mean(exp_arr))
        corr_bkg_cts = new_bkg_arr / (exp_arr/np.mean(exp_arr))
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

            ## Convert blanksky background into an image
            blanksky_image.punlearn()
            blanksky_image.bkgfile = bkg_img
            blanksky_image.outroot = outdir + '{}_bkg_img'.format(obsid)
            blanksky_image.imgfile = cts_img_no_pts
            blanksky_image.clobber = True
            blanksky_image()

            ## Load blanksky image
            with fits.open(outdir + '{}_bkg_img_particle_bgnd.img'.format(obsid)) as f:
                bkg_arr = np.array(f[0].data)
                bkg_hdr = f[0].header

            ## Load exposure map
            with fits.open(exp_map) as f:
                exp_arr = np.array(f[0].data)

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
            nchoices = int(len(indices)/2)
            choices = np.random.choice(indices, size=nchoices, replace=False)
            new_arr = np.zeros_like(arr)
            for choice in choices:
                new_arr[xs[choice]][ys[choice]] += 1
    
            ## Repeat for background
            bkg_arr *= 1000
            xs = []
            ys = []
            indices = []
            i = 0
            for x in range(bkg_arr.shape[0]):
                for y in range(bkg_arr.shape[1]):
                    for z in range(int(bkg_arr[x][y])):
                        xs.append(x)
                        ys.append(y)
                        indices.append(i)
                        i += 1
            nchoices = int(len(indices)/2)
            choices = np.random.choice(indices, size=nchoices, replace=False)
            new_bkg_arr = np.zeros_like(bkg_arr)
            for choice in choices:
                new_bkg_arr[xs[choice]][ys[choice]] += 1
            new_bkg_arr /= 1000

            ## Exposure correct image and blanksky image
            corr_cts_obsid = arr / (exp_arr/np.mean(exp_arr))
            corr_bkg_obsid = bkg_arr / (exp_arr/np.mean(exp_arr))
            fits.writeto(corr_cts_img_obsid, corr_cts_obsid, hdr, overwrite=True)
            fits.writeto(corr_bkg_img_obsid, corr_bkg_obsid, bkg_hdr, overwrite=True)

            ## Add images for all obsids
            if obsid == obsids[0]:
                corr_cts = arr*0.0
                corr_bkg_cts = bkg_arr*0.0

            corr_cts += corr_cts_obsid
            corr_bkg_cts += corr_bkg_obsid

        ## Save summed image and blanksky image
        fits.writeto(corr_cts_img, corr_cts, hdr, overwrite=True)
        fits.writeto(corr_bkg_cts_img, corr_bkg_cts, bkg_hdr, overwrite=True)

def getlimit(indir,obsids,center_x,center_y):
    if len(obsids) == 1:
        cts_img = indir + 'wide_broad_thresh.img'

        dmcoords.punlearn()
        dmcoords.infile = cts_img
        dmcoords.x = center_x
        dmcoords.y = center_y
        dmcoords.opt = 'sky'
        dmcoords()

        chip = int(dmcoords.chip_id)
        if chip == 0 or chip == 3:
            minx = dmcoords.chipx
            miny = dmcoords.chipy
        elif chip == 1 or chip == 2:
            minx = 1024 - dmcoords.chipx
            miny = dmcoords.chipy

        limit = min(minx,miny)

    else:
        limit = 100000
        for obsid in obsids:
            cts_img = indir + 'wide_{}_broad_thresh.img'.format(obsid)

            dmcoords.punlearn()
            dmcoords.infile = cts_img
            dmcoords.x = center_x
            dmcoords.y = center_y
            dmcoords.opt = 'sky'
            dmcoords()

            chip = int(dmcoords.chip_id)
            if chip == 0 or chip == 3:
                minx = dmcoords.chipx
                miny = dmcoords.chipy
            elif chip == 1 or chip == 2:
                minx = dmcoords.chipx
                miny = 1024 - dmcoords.chipy
        
            limit = min(minx,miny,limit)

    return limit

def bkgrate(indir):
    outdir = indir + 'aphot/'
    bkg_img = outdir + 'corrected_bkg_counts.img'
    bkg_reg_file = outdir + 'bkg.reg'

    ## Get counts for blanksky image
    dmstat.punlearn()
    dmstat.infile = '{}[sky=region({})]'.format(bkg_img,bkg_reg_file)
    dmstat.centroid = 'no'
    dmstat()

    ## Calculate background rate for all good pixels in cts/pixels**2
    total = float(dmstat.out_sum)
    area = float(dmstat.out_good)
    bkgrate = total/area
    return bkgrate

def calc(indir,center_x,center_y,R500,bkgrate,limit):
    outdir = indir + 'aphot/'

    bins = np.array([0.05, 0.12, 0.2, 0.30, 1])*R500
    angles = [0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360]

    if limit <= R500:
        bins[-1] = limit

    ## Store regions for later reference
    with open('{}annuli.reg'.format(outdir), 'w') as f:
        f.write('annulus({},{},{},{},{},{},{})'.format(center_x,center_y,bins[0],bins[1],bins[2],bins[3],bins[4]))

    ## Get counts in all angular slices for all bins
    for i in range(len(bins)-1):
        print('Annulus {}...'.format(i+1))
        inner = bins[i]
        outer = bins[i+1]
        counts = []
        total = 0
        for j in range(len(angles)-1):
            print('    Angle {}...'.format(j))
            region = 'pie({},{},{},{},{},{})'.format(center_x,center_y,inner,outer,angles[j],angles[j+1])

            dmstat.punlearn()
            dmstat.infile = '{}corrected_counts.img[(x,y)={}]'.format(outdir,region)
            dmstat.centroid = 'no'
            dmstat()

            area = float(dmstat.out_good)

            total += float(dmstat.out_sum) - bkgrate*area

            counts.append(total)

        fig, ax = plt.subplots()
        ax.plot(np.array(counts)/total)
        plt.savefig('{}plot_{}.png'.format(outdir,i))

#indir = '../results/PSZ2G021.10+33.24/results/'
#obsids = [6104,7940]
#z = 0.151400
#R500 = 1331.756149120992
#center_x = 3910
#center_y = 4072
#N = 1

#indir = '../results/PSZ2G000.13+78.04/results/'
#obsids = [17159]
#z = 0.171000
#R500 = 1150.3347423779778
#center_x = 3926
#center_y = 3766
#N = 1

indir = '../results/PSZ2G024.44+22.76/results/'
obsids = [26045]
z = 0.164700
R500 = 1104.6102507593205
center_x = 4197.5
center_y = 4147.4905792183945
N = 1

def runAphot(indir,obsids,z,R500,center_x,center_y,N):
    R500inpixels = getR500inpixels(z,R500)
    for i in range(N):
        expcorrect(indir,obsids)
        limit = getlimit(indir,obsids,center_x,center_y)
        ## Have user create background region
        if not os.path.exists(indir + 'aphot/bkg.reg'):
            map_file = indir + 'aphot/corrected_counts.img'
            sp.call(['bash', 'shell/ds9_cluster.sh', map_file])
            time.sleep(6)
            sp.call(['bash', 'shell/ds9_visual_improve.sh'])
            input('Add a background region then press Enter to continue...')
            print('------------------------------------------------------------')
            reg_file = indir + 'aphot/bkg.reg'
            sp.call(['bash', 'shell/ds9_save_regions.sh', reg_file])
            time.sleep(3)
            sp.call('killall ds9', shell=True)        
        bkg_rate = bkgrate(indir)
        calc(indir,center_x,center_y,R500inpixels,bkg_rate,limit)

runAphot(indir,obsids,z,R500,center_x,center_y,N)
