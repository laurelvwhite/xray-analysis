#! /bin/bash
infile=$1
regfile=$2
bkg=$3
bkgreg=$4
outfile=$5

punlearn specextract
pset specextract infile=$infile'[sky=region('$regfile')]'
pset specextract bkgfile=$bkg'[sky=region('$bkgreg')]'
pset specextract outroot=$outfile
specextract mode=h binarfwmap=4 bkgresp=no clob+
