#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

ymin, ymax = -0.0000015, 0.000004

r10, r11, r12 = [],[],[]
r20, r21, r22 = [],[],[]

# AE wavefunctions
ffile1 = open('./D_Pae.dat', 'r')
# PS wavefunctions
ffile2 = open('./D_Pus.dat', 'r')

lines1 = ffile1.readlines()
lines2 = ffile2.readlines()
for line in lines1:
    r10.append(float(line.split()[0]))
    r11.append(float(line.split()[1]))
    r12.append(float(line.split()[2]))
   
for line in lines2:
    r20.append(float(line.split()[0]))
    r21.append(float(line.split()[1]))
    r22.append(float(line.split()[2]))

fig = plt.gcf()
fig.set_size_inches(8, 6)
##fig.savefig('test5D-wf.png', dpi=100)

p1 = plt.plot(r20, r21, "blue")
p1 = plt.plot(r20, r22, "red")
p1 = plt.plot(r10, r11, "blue", linestyle='dashed')
p1 = plt.plot(r10, r12, "red", linestyle='dashed')

plt.ylim(ymin, ymax)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
plt.ylabel('Wavefunction', size=15)
plt.xlabel('Radius (bohr)', size=15)
plt.legend(('PP_E=ae', 'PP_E=ref', 'AE_E=ae', 'AE_E=ref'), loc='upper right')

ref=-1.00
rc=3.00

pl = plt.vlines([rc], ymin, ymax, linestyle='dotted')

plt.title("5D WF: ref = "+str(ref)+", rc = "+str(rc))

save="../figures/WF_5D_ref="+str(ref)+"_rc="+str(rc)+".png"
fig.savefig(save, dpi=100)
##plt.show()
