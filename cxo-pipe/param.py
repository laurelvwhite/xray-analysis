source_name = "PSZ2G024.44+22.76"
obsids = "26045"
z = 0.164700
R500 = 1104.6102507593205
use_peak = True
fixed_coord = None
fast_annuli = True
Ysz = [6e-05, 1e-05]
single_ann_spec = False
fixed_spec_annuli = False
fit_kT_profile_directly = True
file_ACCEPT = None
compute_Lcool = False
tcool_th = 7.7
input_XSZ_file = None
do_err = True
tab_obsid = obsids.split(",")
if len(tab_obsid) > 1:
    multiobs = True
else:
    multiobs = False
