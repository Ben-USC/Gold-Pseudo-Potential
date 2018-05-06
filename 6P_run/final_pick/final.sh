#!/bin/bash

start=`date +%s`

ref="0.05"
rc="2.9"

s="d0"
s0="_rc="
out="_"
refs=$ref$s
rcs=$rc$s
log=$ref$out$rc

rm AU_*
rm S_*
rm P_*
rm D_*
rm vloc.dat
rm ele.dat
rm out_AU

sed -e "s/\-5.1d0/${refs}/g;s/2.9d0/${rcs}/g" in7.dat.false_true > in7.dat
./atm > ./log_files/atm_false_true_out_$log

sed -e "s/\-5.1d0/${refs}/g;s/2.9d0/${rcs}/g" in7.dat.true_false > in7.dat
./atm > ./log_files/atm_true_false_out_$log

sed -e "s/\-5.1/${ref}/g;s/2.9/${rc}/g" plot.py.temp > plot.py
python3 plot.py

end=`date +%s`
runtime=$((end-start))
echo "------------- Run time: $runtime s ------------"
