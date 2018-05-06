#!/bin/bash

s="d0"
s0="_rc="

for ref in $(seq -2.0 .05 -1.0)
do
   for rc in $(seq 1.0 .05 3.0)
   do
      ##echo $ref $rc
      refs=$ref$s
      rcs=$rc$s
      ##sed -e "s/\-1.5d0/${refs}/g;s/2.0d0/${rcs}/g" in7.dat.temp > in7.dat_ref=$ref$s0$rc
      sed -e "s/\-1.5d0/${refs}/g;s/2.0d0/${rcs}/g" in7.dat.temp > in7.dat
      ./atm
   done
done
