#!/bin/bash

start=`date +%s`

s="5d_ref\="
s0="_rc\="
s1=".eps"
s2=".jpg"

##for ref in $(seq -1.0 .05 0.0)
##do
##   for rc in $(seq 2.0 .05 3.5)
##   do
##      ss1=$s$ref$s0$rc$s1
##      ss2=$s$ref$s0$rc$s2
##      convert ./5D_figures/$ss1 ./5D_figures_converted/$ss2
##   done
##done

for eps in ./5D_figures/*.eps; do
    filename=${eps%.*}
    convert "$filename.eps" "$filename.jpg"
done

mv ./5D_figures/*.jpg ./5D_figures_converted/

end=`date +%s`
runtime=$((end-start))
echo "------------- Run time: $runtime s ------------"
