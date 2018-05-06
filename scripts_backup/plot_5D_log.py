import numpy as np
import matplotlib.pyplot as plt

ymin, ymax = -4, 4

r10, r11, r12 = [],[],[]
r20, r21, r22 = [],[],[]

# AE wavefunctions
ffile1 = open('D_chiae.dat', 'r')
# PS wavefunctions
ffile2 = open('D_chil.dat', 'r')

lines1 = ffile1.readlines()
lines2 = ffile2.readlines()
for line in lines1:
    r10.append(float(line.split()[0]))
    r11.append(float(line.split()[1]))
   
for line in lines2:
    r20.append(float(line.split()[0]))
    r21.append(float(line.split()[1]))

fig = plt.gcf()
fig.set_size_inches(6.5, 4.5)

p1 = plt.plot(r20, r21, "blue")
p1 = plt.plot(r10, r11, "red", linestyle='dotted')
plt.ylim(ymin, ymax)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
plt.ylabel('Logarithmic derivative', size=15)
plt.xlabel('Energy (Ry)', size=15)
plt.legend(('PP', 'AE'), loc='upper right')

ref=-1.00
rc=3.00

plt.title("5D LG: ref = "+str(ref)+", rc = "+str(rc))

save="../figures/log_derivative_5D_ref="+str(ref)+"_rc="+str(rc)+".png"
fig.savefig(save, dpi=100)
