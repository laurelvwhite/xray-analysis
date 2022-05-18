#! /bin/bash
mapfile=$1
outfile=$2

blanksky.punlearn()
            blanksky.evtfile = '{}/prep/lightcurve/{}_evt2_lc.fits'.format(datadir,obsid)
            blanksky.outfile = '{}/prep/background/{}_bkg.evt'.format(datadir,obsid)
            blanksky.tmpdir = '{}/prep/background/'.format(datadir)
            blanksky.asolfile = '{}/{}/repro/*asol1.fits'.format(datadir,obsid)
            blanksky.clobber = 'yes'
            blanksky()

punlearn blanksky
pset blanksky evtfile=$mapfile
pset blanksky outfile=$outfile
pset blanksky clobber="yes"
blanksky
