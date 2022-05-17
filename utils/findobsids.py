#!/usr/bin/env python

import argparse

from ciao_contrib.runtool import *

parser = argparse.ArgumentParser()
parser.add_argument('--ra')
parser.add_argument('--dec')
args = parser.parse_args()

ra = args.ra
dec = args.dec

def listobsids(ra, dec):
    obsidoutput = find_chandra_obsid(ra,dec,instrument='acisi',detail='obsid')
    obsidlist = str(obsidoutput).split()
    obsidlist = obsidlist[2:len(obsidlist)]
    return obsidlist

print(listobsids(ra,dec))
