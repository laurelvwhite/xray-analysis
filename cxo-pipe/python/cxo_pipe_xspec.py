import subprocess as sp
from termcolor import colored
from astropy.io import fits
import os
import logging
import time
import numpy as np
import copy
import glob
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import ascii
from gdpyc import GasMap
import warnings
from tqdm import tqdm
import scipy.optimize as optimization
import shutil
import time
from sherpa.astro import ui
#from sherpa_contrib.utils import *
from xspec import *

def fit_spec_pyxspec(res_dir, obsids, z):

    mer_dir = os.path.join(res_dir, 'results') #change directory name
    cl_dir  = os.path.join(mer_dir, 'cluster')
    bkg_dir = os.path.join(mer_dir, 'background')
    fig_dir = mer_dir + "/figures/"

    if not os.path.exists(fig_dir):
        sp.call("mkdir " + fig_dir, shell=True)

    cl_spectra = sorted(glob.glob(os.path.join(cl_dir, 'cl_spectrum*')))
    bkg_spectra = sorted(glob.glob(os.path.join(bkg_dir, 'bkg_spectrum*')))

#    tmpdir = os.path.join(mer_dir, 'xspec')
    tmpdir = '/Users/laurelwhite/temp/{}/'.format(np.random.randint(1,100000000))

    try:
        os.mkdir(tmpdir)
    except FileExistsError:
        print('Directory already exists. Continuing...')
    cwd = os.getcwd()
    
    for filename in cl_spectra+bkg_spectra:
        shutil.copy(filename,tmpdir)

    file_T_prof = os.path.join(cl_dir, "spectro_T_prof.npz")

    if os.path.exists(file_T_prof):

        print(colored("Spectra already fitted", "white", None, ["bold"]))
        print("------------------------------------------------------------")

    else:  

        efile = os.path.join(mer_dir, 'efile_repro_raw_clean_nopts.fits')
        reg_file_cl = os.path.join(cl_dir, "spec_annuli.reg")
        efile_in = efile + "[bin sky=@" + reg_file_cl + "]"
        file_out = os.path.join(mer_dir, "spec_ann_stat.fits")
        sp.call(["bash", "shell/extract_content.sh", efile_in, file_out])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hdu = fits.open(file_out)
            tab_area_cl = hdu[1].data["AREA"]

        tab_obsid = obsids.split(",")
        
        logger = logging.getLogger('xspec')
        logger.setLevel(logging.WARN)
        t0 = time.time()

        N_spec_reg = tab_area_cl.size
        ICM_T_tab = np.zeros(N_spec_reg)
        ICM_T_tab_err = np.zeros(N_spec_reg)
        tab_area_bkg = [0]

        for ind_ann in tqdm(range(1, N_spec_reg + 1)):

            buff_T_err = 1.0
            buff_T = 1.0
            SNR_bin = 3
            T_max = 100.0

            fit_ind = 1
            for obsid in tab_obsid:
                file_area_bkg = os.path.join(mer_dir,  "bkg_area_" + obsid + ".txt")
                with open(file_area_bkg) as f:
                    content = f.readlines()
                area_bkg = float(content[0])
                tab_area_bkg.append(area_bkg)
                ui.load_data(
                    fit_ind,
                    os.path.join(cl_dir, "cl_spectrum_" + obsid + "_" + str(ind_ann) + ".pi"),
                )
                fit_ind += 1
#                tab_area_bkg.append(area_bkg)
#                ui.load_data(fit_ind, os.path.join(bkg_dir, "bkg_spectrum_" + obsid + ".pi"))
#                fit_ind += 1

            for i in range(1, fit_ind):
                full_spec = ui.get_data_plot(i)
                particle_spec = ui.get_bkg(i)
                wscale = np.where((full_spec.x > 9) & (full_spec.x < 12))
                newscale = np.trapz(
                    ui.get_counts(i)[wscale], full_spec.x[wscale]
                ) / np.trapz(particle_spec.counts[wscale], full_spec.x[wscale])
                if newscale != 0:
                    area_scale = ui.get_backscal(i)
                    bkg_scale = ui.get_bkg_scale(i)
                    ui.set_backscal(i, area_scale * newscale / bkg_scale)
                ui.subtract(i)
                ui.group_snr(i, SNR_bin)

            # Get/Set the console chatter level
            ch = Xset.chatter
            Xset.chatter = 5
            # Get/Set the log chatter level
            lch = Xset.logChatter
            Xset.logChatter = 10

            # Create and open a log file for XSPEC output
            # This returns a Python file object
            logFile = Xset.openLog(os.path.join(tmpdir, f"xspec_log_{ind_ann}.txt"))
            # Get the Python file object for the currently opened log
            logFile = Xset.log

            os.chdir(tmpdir)

            files = glob.glob('{tmpdir}/*')
            for f in files:
                os.remove(f)

            AllData.clear()
            AllData(' '.join([f'{i+1}:{i+1} {cl_dir}/cl_spectrum_{obs}_{ind_ann}.pi' for i, obs in enumerate(tab_obsid)]))

            AllData.ignore("bad")
            AllData.ignore("**-0.7")
            AllData.ignore("7.0-**")

            AllModels.setEnergies("0.01 100. 10000")

            Xset.cosmo = "70,,0.7"
            Xset.abund = 'angr'
            Xset.parallel.leven    = 8
            Xset.parallel.error    = 8
            Xset.parallel.walkers  = 8
            Xset.parallel.steppar  = 8
            Xset.parallel.goodness = 8

            Plot.device = "/null"
            Plot.xAxis = "keV"

            AllChains.clear()
            AllChains.defLength  = 10000
            AllChains.defBurn    = 5000
            AllChains.defWalkers = 8
            AllChains.defRescale = 0.125

            Fit.query = "yes"
            Fit.statMethod = 'cstat'
            Fit.nIterations = 100

            nH_val = float(np.load(os.path.join(mer_dir, "nH_value.npy")))

            m = Model('phabs*apec')
            m.componentNames

            m.phabs.nH = f'{nH_val} -1'
            m.apec.Redshift = f'{z} -1'
            m.apec.Abundanc = '0.3 -1'
            m.apec.kT = '6.0 0.01 0.01 1 15 20'
            m.apec.norm = 7e-4

            for i in range(1, fit_ind):
                tab_area_fact = tab_area_cl / tab_area_bkg[i]

#            bkg_dict = {}
#            for i, bkg in enumerate(tab_obsid):
#                sp_index  = 2*i + 2
#                bkg_index = f'm{sp_index}'
#                bkg_dict[bkg_index] = AllModels(sp_index)

#            for modname, modobj in bkg_dict.items():
#                modobj.apec.norm = '0 -1'
#                modobj.bremss.norm.link = f'p7 / {tab_area_fact[ind_ann - 1]}'
#                modobj.apec_4.norm.link = f'p11 / {tab_area_fact[ind_ann - 1]}'

            Fit.perform()

            mcmc_chain_file = f'{tmpdir}/chain_{ind_ann}.fits'
            mcmc_chain = Chain(mcmc_chain_file)

#            Fit.error("1. 2 5 7 11")

            with fits.open(mcmc_chain.fileName) as mcmc_result:
                mcmc_chain_array = np.array(mcmc_result[1].data.tolist())

            final_temp      = np.percentile(mcmc_chain_array[:,0], 50)
            final_temp_erru = np.percentile(mcmc_chain_array[:,0], 84) - final_temp
            final_temp_errd = final_temp - np.percentile(mcmc_chain_array[:,0], 16)

            ICM_T_tab[ind_ann - 1] = final_temp
            buff_T = final_temp
            T_max = final_temp

            if (final_temp_errd is None) & (final_temp_erru is None):
                ICM_T_tab_err[ind_ann - 1] = final_temp
                buff_T_err = final_temp
            elif (final_temp_errd is not None) & (final_temp_erru is None):
                ICM_T_tab_err[ind_ann - 1] = np.abs(final_temp_errd)
                buff_T_err = np.abs(final_temp_errd)
            else:
                ICM_T_tab_err[ind_ann - 1] = final_temp_erru
                buff_T_err = final_temp_erru

            AllChains.clear()

            # Close XSPEC's currently opened log file.
            Xset.closeLog()
            Plot.device = "/null"


        np.savez(file_T_prof, T=ICM_T_tab, Terr=ICM_T_tab_err)
        print(f'(xspec) total elapsed time: {(time.time()-t0)/60:.2f} min \n-----------------------------------')

        tmp_spec_files = glob.glob(tmpdir+'/*spectrum*')
        for tmp_specfile in tmp_spec_files:
            os.system('rm {}'.format(tmp_specfile))

        os.chdir(cwd)
