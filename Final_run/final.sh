#!/bin/bash

start=`date +%s`

rm AU_*
rm S_*
rm P_*
rm D_*
rm vloc.dat
rm ele.dat
rm out_AU

cp ./in7.dat.false_true ./in7.dat

./atm > ./log_files/atm_false_true_out

cp ./in7.dat.true_false ./in7.dat

./atm > ./log_files/atm_true_false_out

python3 plot_5D.py

python3 plot_6S.py

python3 plot_6P.py

mv AU_* ./pseudo_potential_files
mv S_* P_* D_* ./output_files
mv ele.dat ./output_files
mv vloc.dat ./output_files

end=`date +%s`
runtime=$((end-start))
echo "------------- Run time: $runtime s ------------"
