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


print("------------------------------------------------------------")
print(colored("Pre-processing", "cyan", None, ["bold"]))
print("------------------------------------------------------------")
res_dir = os.environ["CXO_RES_DIR"] + source_name + "/"
# Import data
prepro.import_data(obsids, res_dir)
# Reprocess data
prepro.reprocess_data(obsids, res_dir)
# Apply energy cuts
prepro.energy_cut(obsids, res_dir)
# Remove flares
prepro.remove_flares(obsids, res_dir)
# Compute maps of the observations
if multiobs:
    prepro.reproj_obsids(obsids, res_dir)
else:
    prepro.make_images(obsids, res_dir)
    prepro.make_psf_map(res_dir)
# Find point sources
is_sources = prepro.find_sources(res_dir, multiobs)
# Check point source regions
prepro.check_sources(res_dir, is_sources)
# Compute a background region for each point source
prepro.bkg_point_sources(res_dir, is_sources)
# Subtract point sources
prepro.subtract_point_sources(res_dir, is_sources, multiobs, obsids)
# Find X-ray peak and centroid locations
Xdepro, Ydepro = prepro.find_peak_cent(res_dir, z, R500, use_peak, fixed_coord)
# Define background region for X-ray surface brightness and spectra
bkg_area = prepro.bkg_region(res_dir, z, R500, Xdepro, Ydepro, multiobs, obsids)
# Define annuli for the surface brightness profile
prepro.find_SB_annuli(res_dir, Xdepro, Ydepro, bkg_area, z, R500, fast_annuli, obsids)
# Compute weights to take vignetting into account
prepro.vignetting_prof(res_dir, obsids)
# Compute vignetting-corrected surface brightness profile
prepro.X_ray_SB_profile(res_dir, obsids, z)

print("------------------------------------------------------------")
print(colored("Spectral analysis", "cyan", None, ["bold"]))
print("------------------------------------------------------------")
## Extract background spectra
#spec.bkg_spectrum(res_dir, multiobs, obsids)
# Define the annuli to be used for cluster spectrum extraction
N_ann = spec.find_spec_annuli(
    res_dir, Xdepro, Ydepro, bkg_area, z, R500, single_ann_spec, obsids
)
# Extract a spectrum in each annulus for each obsid
spec.extract_cl_spectra(res_dir, multiobs, obsids)
# Compute the ARF and RMF in each annulis of the SB profile
spec.arf_per_bin(res_dir, multiobs, obsids)
# Compute the hydrogen column density in the cluster direction
spec.hydrogen_column(res_dir)
