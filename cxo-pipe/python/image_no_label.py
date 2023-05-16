import os
import numpy as np
import subprocess as sp

from ciao_contrib.runtool import *
from cosmocalc import cosmocalc
from termcolor import colored

def make_unlabeled_image(clus):
    base_dir = os.environ['CXO_RES_DIR']
    res_dir = base_dir + clus + '/results/'

    ## Make label region
    with open(res_dir + 'label.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\n')
        f.write('global color=white dashlist=8 3 width=1 font="helvetica 16 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('image\n')

    sp.call(['bash', 'shell/labeled_image.sh', res_dir + 'figures/cropped_smoothed.img', res_dir + 'label.reg', res_dir + 'figures/unlabeled.png'])
