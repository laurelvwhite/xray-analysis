#! /bin/bash
infile=$1
regfile=$2
outfile=$3

punlearn specextract
pset specextract infile=$infile'[sky='$regfile']'
pset specextract outroot=$outfile
pset specextract binspec=2
specextract mode=h binarfwmap=4 bkgresp=no clob+
