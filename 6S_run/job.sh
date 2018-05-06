#!/bin/bash

start=`date +%s`

s="d0"
s0="_rc="
out="_"

for ref in $(seq -1.0 0.05 0.0)
do
   for rc in $(seq 2.0 .05 3.5)
   do
      ##echo $ref $rc
      refs=$ref$s
      rcs=$rc$s
      log=$ref$out$rc
      ##sed -e "s/\-1.5d0/${refs}/g;s/2.0d0/${rcs}/g" in7.dat.temp > in7.dat_ref=$ref$s0$rc
      sed -e "s/\-5.1d0/${refs}/g;s/2.6d0/${rcs}/g" in7.dat.false_true > in7.dat
      ./atm > ./log_files/atm_false_true_out_$log

      sed -e "s/\-5.1d0/${refs}/g;s/2.6d0/${rcs}/g" in7.dat.true_false > in7.dat
      ./atm > ./log_files/atm_true_false_out_$log

      sed -e "s/\-5.1/${ref}/g;s/2.6/${rc}/g" plot.py.temp > plot.py
      python3 plot.py

   done
done

end=`date +%s`
runtime=$((end-start))
echo "------------- Run time: $runtime s ------------"
