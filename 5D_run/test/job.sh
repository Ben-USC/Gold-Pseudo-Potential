#!/bin/bash

start=`date +%s`

cp in7.dat.false_true in7.dat
./atm > ./atm_false_true_out

cp in7.dat.true_false in7.dat
./atm > ./atm_true_false_out

python3 plot.py

end=`date +%s`
runtime=$((end-start))
echo "------------- Run time: $runtime s ------------"
