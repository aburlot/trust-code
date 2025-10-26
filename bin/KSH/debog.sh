#!/bin/bash
[ "$cpu" = "" ] && echo "cpu empty." && exit -1
jdd=`pwd`
jdd=`basename $jdd`
cp $jdd.data cpu.data
cp $jdd.data gpu.data
sed -i "1,$ s?Solve?Debog pb seq faces 1.e-6 0 Solve?g" cpu.data
sed -i "1,$ s?Solve?Debog pb seq faces 1.e-6 1 Solve?g" gpu.data
(source $cpu/env_TRUST.sh;$exec  cpu 2>&1 | tee cpu.out_err)
$exec gpu 2>&1 | tee gpu.out_err
compare_lata cpu.lml gpu.lml
