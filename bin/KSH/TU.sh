#!/bin/bash
# Support old and new TU format
TU=$1 && [ ! -f $TU ] && echo "$TU not found!" && exit -1
[ "$2" = "" ] && echo "Specify an option" && echo "TU.sh file [-dt] [-its]" && exit -1
if [ $2 = "-dt" ]
then
   # Time step
   awk '/Secondes/ && /pas de temps/ {val=$NF} /Average time per time step/ {val=$NF} END {print val}' $TU
elif [ $2 = "-its" ]
then
   # Iterations solver
   awk '/Iterations/ && /solveur/ {val=$NF} /Average number of iteration of the linear solver per call/ {val=$NF} END {print val}' $TU 
else
   echo "Option $2 unknown"
   exit -1
fi
