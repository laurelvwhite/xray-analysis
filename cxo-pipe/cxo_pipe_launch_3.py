import os
import time
from termcolor import colored
import subprocess as sp


import python.cxo_pipe_preproc as prepro
import python.cxo_pipe_spec as spec
import python.cxo_pipe_icm as icm
import python.cxo_pipe_plot as plt
import python.concentration as conc
import python.centroid as cent
import python.image as image

from param import *

ti = time.time()
res_dir = os.environ["CXO_RES_DIR"] + source_name + "/"

# Get center
Xdepro, Ydepro = prepro.find_peak_cent(res_dir, z, R500, use_peak, fixed_coord)
# Get background area
bkg_area = prepro.bkg_region(res_dir, z, R500, Xdepro, Ydepro, multiobs, obsids)
# Get number of annuli
N_ann = spec.find_spec_annuli(
    res_dir, Xdepro, Ydepro, bkg_area, z, R500, single_ann_spec, obsids
)
# Plot the fit results
plt.plot_spectra(res_dir, obsids)
# Fit the temperature profile with a Vikhlinin model
spec.fit_T_prof(res_dir, R500, z)
# Plot the fit results
plt.plot_T_prof(res_dir, R500)
# Compute the conversion factors from surface brightness to emission measure
spec.XSB_to_EM_coef(res_dir, obsids, z)

print("------------------------------------------------------------")
print(colored("Estimation of ICM profiles", "cyan", None, ["bold"]))
print("------------------------------------------------------------")
# Initialize log-spaced integration radius maps
Rproj, r_map, los_step_map = icm.init_integ_maps(z, R500)
# Run the MCMC analysis to fit the surface brightness profile
icm.mcmc_ne(res_dir, Rproj, r_map, los_step_map, z, R500)
# Clean the chains
icm.clean_chains(res_dir, "ne")
# Find the best-fit density model
icm.best_ne_model(res_dir, Rproj, r_map, los_step_map, z, R500)
if fit_kT_profile_directly:
    icm.mcmc_kT(res_dir, R500)
    icm.clean_chains(res_dir, 'kT')
    icm.best_icm_models_kTdirect(res_dir, z, R500, N_ann, Ysz)
else:
    if N_ann > 2:
        # Run the MCMC analysis to fit the temperature profile
        icm.mcmc_pe(res_dir)
        # Clean the chains
        icm.clean_chains(res_dir, "pe")
    # Find the best-fit ICM models
    icm.best_icm_models(res_dir, z, R500, N_ann, Ysz)
# Compute the cooling luminosity if requested
if compute_Lcool:
    icm.cooling_lum(
        res_dir, z, tcool_th, Xdepro, Ydepro, multiobs, obsids, input_XSZ_file, do_err
    )

print("------------------------------------------------------------")
print(colored("Analysis figures", "cyan", None, ["bold"]))
print("------------------------------------------------------------")

# Plot the ICM profiles
plt.plot_icm_profiles(res_dir, file_ACCEPT, z)
plt.plot_2D_posteriors(res_dir, N_ann, fit_kT_profile_directly)
plt.adaptive_map(res_dir, z, R500)
plt.compute_Aphot(res_dir, z, R500, bkg_area, obsids)
plt.cluster_id_card(res_dir, source_name, z)

conc.calc_concentration(source_name, z, R500)
cent.calc_centroid_shift(source_name, z, R500)
image.make_labeled_image(source_name, z, R500)

sp.call("cp param.py " + res_dir, shell=True)

te = time.time()
print(
    colored("=== *************************************** ===", "cyan", None, ["bold"])
)
print(colored("               End of program", "cyan", None, ["bold"]))
print(
    colored(
        "          Execution time: " + "{0:.1f}".format((te - ti) / 60.0) + " min",
        "cyan",
        None,
        ["bold"],
    )
)
print(
    colored("=== *************************************** ===", "cyan", None, ["bold"])
)
