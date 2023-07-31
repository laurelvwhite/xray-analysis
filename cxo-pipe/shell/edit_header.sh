#! /bin/bash
infile=$1
key=$2
value=$3

punlearn dmhedit
pset dmhedit infile=$infile
pset dmhedit filelist="none"
pset dmhedit operation="add"
pset dmhedit key=$key
pset dmhedit value=$value
{
  expect << EOD
  set timeout -1
  spawn dmhedit
  expect {
    "Input dataset*" { 
      send "\n" 
      exp_continue 
    }
    "Edit list file name*" { 
      send "\n" 
      exp_continue 
    }
    "Operation*" { 
      send "\n" 
      exp_continue 
    }
    "Keyword*" { 
      send "\n" 
      exp_continue 
    }
    "Value*" { 
      send "\n" 
      exp_continue 
    }
  }
EOD
}
