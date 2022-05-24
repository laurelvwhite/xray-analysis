#! /bin/bash
mapfile=$1
outfile=$2
asolfile=$3

punlearn blanksky
pset blanksky evtfile=$mapfile
pset blanksky outfile=$outfile
pset blanksky asolfile=$asolfile
{
  expect << EOD
  set timeout -1
  spawn blanksky clobber=yes
  expect {
    "Source event file*" { 
      send "\n" 
      exp_continue 
    }
    "Input event file*" { 
      send "\n" 
      exp_continue 
    }
    "Enter output file*" { 
      send "\n" 
      exp_continue 
    }
    "Output directory path*" { 
      send "\n" 
      exp_continue 
    }
  }
EOD
}
