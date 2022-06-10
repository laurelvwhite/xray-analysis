source_name = "0152−1358"
obsids = "913,21703,22856"
z = 0.8325
R500 = 813.9644434605176
use_peak = True
fixed_coord = None
fast_annuli = True
Ysz = [6e-05, 1e-05]
single_ann_spec = False
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
