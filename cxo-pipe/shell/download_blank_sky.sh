#! /bin/bash
mapfile=$1
outfile=$2
asolfile=$3

punlearn blanksky
pset blanksky evtfile=$mapfile
pset blanksky outfile=$outfile
pset blanksky asolfile=$asolfile
pset blanksky clobber="yes"
blanksky
