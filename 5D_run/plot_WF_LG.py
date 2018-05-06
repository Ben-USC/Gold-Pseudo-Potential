import numpy as np
import matplotlib.pyplot as plt

ymin1,ymax1 = -0.0000015,0.000004
r10, r11, r12 = [],[],[]
r20, r21, r22 = [],[],[]

# AE wavefunctions
file1 = open('./D_Pae.dat', 'r')
# PS wavefunctions
file2 = open('./D_Pus.dat', 'r')

lines1 = file1.readlines()
lines2 = file2.readlines()
for line in lines1:
    r10.append(float(line.split()[0]))
    r11.append(float(line.split()[1]))
    r12.append(float(line.split()[2]))
   
for line in lines2:
    r20.append(float(line.split()[0]))
    r21.append(float(line.split()[1]))
    r22.append(float(line.split()[2]))

ref=-0.50
rc=3.00

fig,axes=plt.subplots(1,2)
fig.set_size_inches(10, 4)

axes[0].plot(r20, r21, "blue")
axes[0].plot(r20, r22, "red")
axes[0].plot(r10, r11, "blue", linestyle='dashed')
axes[0].plot(r10, r12, "red", linestyle='dashed')
axes[0].set_ylim(ymin1, ymax1)
axes[0].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
axes[0].set_ylabel('Wavefunction', size=12)
axes[0].set_xlabel('Radius (bohr)', size=12)
axes[0].legend(('PP_E=ae', 'PP_E=ref', 'AE_E=ae', 'AE_E=ref'), loc='upper right')
axes[0].vlines([rc], ymin1, ymax1, linestyle='dotted')
axes[0].set_title("5D WF: ref = "+str(ref)+", rc = "+str(rc))


ymin2, ymax2 = -4, 4
rr10, rr11 = [],[]
rr20, rr21 = [],[]

# AE wavefunctions
ffile1 = open('D_chiae.dat', 'r')
# PS wavefunctions
ffile2 = open('D_chil.dat', 'r')

lines1 = ffile1.readlines()
lines2 = ffile2.readlines()
for line in lines1:
    rr10.append(float(line.split()[0]))
    rr11.append(float(line.split()[1]))
   
for line in lines2:
    rr20.append(float(line.split()[0]))
    rr21.append(float(line.split()[1]))

axes[1].plot(rr20, rr21, "blue")
axes[1].plot(rr10, rr11, "red", linestyle='dotted')
axes[1].set_ylim(ymin2, ymax2)
axes[1].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
axes[1].set_ylabel('Logarithmic derivative', size=12)
axes[1].set_xlabel('Energy (Ry)', size=12)
axes[1].legend(('PP', 'AE'), loc='upper right')
axes[1].set_title("5D LG: ref = "+str(ref)+", rc = "+str(rc))

plt.tight_layout()
save="../5D_figures/5D_ref="+str(ref)+"_rc="+str(rc)+".png"
fig.savefig(save, dpi=100)
