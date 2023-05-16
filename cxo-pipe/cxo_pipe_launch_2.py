import os
import time
from termcolor import colored
import subprocess as sp

import python.cxo_pipe_xspec as xspec

from param import *

res_dir = os.environ["CXO_RES_DIR"] + source_name + "/"

# Simultaneously fit the cluster and background spectra
xspec.fit_spec_pyxspec(res_dir, obsids, z)
