#!/bin/bash

# base name of this script (without extension)
scriptname=$(basename $0 .sh)

echo -e "\n*** Running batch script '"$scriptname"' ***\n"

# path to Matlab code directory
codedir=$LOCALREPO/matlab/DI_opt_02

# current directory (the directory this script is run in - log files will go here)
currdir=$(pwd -P)

# Restrict cores?
if [ -n "$1" ]; then
	tset="taskset -c $1"
fi

# Matlab invocation
runmatlab="nohup nice $tset matlab -nojvm -nodisplay"

# Parameters
net=tnet9d

# run multiple concurrent Matlab sessions with different parameters
for m in {1..8}
do
    # the log file
    logfile=$currdir/$scriptname\_m\_$m.log

    # Matlab commands
    matcmds="\
        resdir = '$currdir';\
        rid    = '_"$net"_$m';\
        G      = $net;\
        r      = 16;\
        m      = $m;\
        niters = 10000;\
        nruns  = 100;\
        hist   = false;\
        gvdisp = false;\
        gpplot = false;\
        mseed  = 1234;\
        iseed  = 5678;\
        oseed  = 9012;\
        sim_opt_es_ssdd;\
        quit"

    # run Matlab as background job
    cd $codedir && $runmatlab -r "$matcmds" > $logfile < /dev/null 2>&1 &
done
