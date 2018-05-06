#!/bin/bash

start=`date +%s`

ref="-0.3"
rc="2.70"

s="d0"
s0="_rc="
out="_"
refs=$ref$s
rcs=$rc$s
log=$ref$out$rc

sed -e "s/\-1.5d0/${refs}/g;s/2.0d0/${rcs}/g" in7.dat.false_true > in7.dat
./atm > ./log_files/atm_false_true_out_$log

sed -e "s/\-1.5d0/${refs}/g;s/2.0d0/${rcs}/g" in7.dat.true_false > in7.dat
./atm > ./log_files/atm_true_false_out_$log

sed -e "s/\-1.5/${ref}/g;s/2.0/${rc}/g" plot.py.temp > plot.py
python3 plot.py

mv AU_* ./output_files/
mv S_* P_* D_* ./output_files/
mv out_AU ./output_files/
mv vloc.dat ./output_files/
mv ele.dat ./output_files/

end=`date +%s`
runtime=$((end-start))
echo "------------- Run time: $runtime s ------------"
