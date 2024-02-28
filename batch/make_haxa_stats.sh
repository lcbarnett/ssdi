#!/bin/bash

# base name of this script (without extension)
scriptname=$(basename $0 .sh)

echo -e "\n*** Running batch script '"$scriptname"' ***\n"

# path to Matlab code directory
codedir=~/git/ssdi

# current directory (the directory this script is run in - log files will go here)
currdir=$(pwd -P)

logfile=$currdir/$scriptname.log

# Parameters

N=1000000
nmax=100

matcmds="clear; make_haxa_stats($nmax,$N,'$currdir'); quit"

# run Matlab

cd $codedir && nohup nice matlab -nojvm -nodisplay -r "$matcmds" > $logfile < /dev/null 2>&1 &
