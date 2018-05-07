#!/bin/bash

start=`date +%s`

s="d0"

ref5d="-0.3"
rc5d="2.7"
ref5ds=$ref5d$s
rc5ds=$rc5d$s

ref6s="-0.2"
rc6s="2.6"
ref6ss=$ref6s$s
rc6ss=$rc6s$s

ref6p="0.05"
rc6p="2.9"
ref6ps=$ref6p$s
rc6ps=$rc6p$s

sed -e "s/\-5.1d0/${ref5ds}/g;s/5.5d0/${rc5ds}/g;
        s/\-5.2d0/${ref6ss}/g;s/5.6d0/${rc6ss}/g;
        s/\-5.3d0/${ref6ps}/g;s/5.7d0/${rc6ps}/g" in7.dat.false_true > in7.dat

./atm > ./log_files/atm_false_true_out

sed -e "s/\-5.1d0/${ref5ds}/g;s/5.5d0/${rc5ds}/g;
        s/\-5.2d0/${ref6ss}/g;s/5.6d0/${rc6ss}/g;
        s/\-5.3d0/${ref6ps}/g;s/5.7d0/${rc6ps}/g" in7.dat.true_false > in7.dat

./atm > ./log_files/atm_true_false_out

sed -e "s/\-5.1/${ref5d}/g;s/5.5/${rc5d}/g" plot_5D.py > plot.py
python3 plot.py

sed -e "s/\-5.2/${ref6s}/g;s/5.6/${rc6s}/g" plot_6S.py > plot.py
python3 plot.py

sed -e "s/\-5.3/${ref6p}/g;s/5.7/${rc6p}/g" plot_6P.py > plot.py
python3 plot.py

mv AU_* ./pseudo_potential_files
mv S_* P_* D_* ./output_files
mv ele.dat ./output_files
mv vloc.dat ./output_files

cp *.png ../

end=`date +%s`
runtime=$((end-start))
echo "------------- Run time: $runtime s ------------"
