#!/bin/bash

# base name of this script (without extension)
scriptname=$(basename $0 .sh)

echo -e "\n*** Running batch script '"$scriptname"' ***\n"

# path to Matlab code directory
codedir=~/git/ssdi

# current directory (the directory this script is run in - log files will go here)
currdir=$(pwd -P)

# Matlab invocation
runmatlab="nohup nice matlab -nojvm -nodisplay"

# Parameters

N=1000000
C=10

# run multiple concurrent Matlab sessions
for n in 10 20 30 40 50 60 70 80 90; do

	# the log file
	logfile=$currdir/$scriptname\_$n.log

	# Matlab commands
	matcmds="clear; gen_haxa_stats($n,$N,$C,'$currdir'); quit"

	# run Matlab
	cd $codedir && $runmatlab -r "$matcmds" > $logfile < /dev/null 2>&1 &

done
