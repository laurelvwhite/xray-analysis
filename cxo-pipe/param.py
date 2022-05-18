source_name = "PSZ2G065.32-64.84"
obsids = "922,6934,7329,19596,19597,19598,20626,20627,20628,20629,20805,20806,20811,20817"
z = 0.085200
R500 = 930.2477408897814
use_peak = True
fixed_coord = None
fast_annuli = True
Ysz = [6e-05, 1e-05]
single_ann_spec = False
file_ACCEPT = None
compute_Lcool = True
tcool_th = 7.7
input_XSZ_file = None
do_err = True
tab_obsid = obsids.split(",")
if len(tab_obsid) > 1:
    multiobs = True
else:
    multiobs = False
