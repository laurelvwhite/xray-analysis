source_name = "PSZ2G002.82+39.23"
obsids = "26041"
z = 0.153300
R500 = 1201.8868335007987
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
