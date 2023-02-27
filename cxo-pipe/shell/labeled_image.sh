#! /bin/bash
infile=$1
regfile=$2
imgfile=$3

ds9 -file $infile -log -scalelims 0.000000005 0.00001 -scale log exp 100000 -cmap magma -cmap 2.5 0.5 -colorbar no -region $regfile -zoom TO FIT -saveimage png $imgfile -exit
