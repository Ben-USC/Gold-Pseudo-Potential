import numpy as np
import matplotlib.pyplot as plt

ref=-5.1
rc=5.5

fig,axes=plt.subplots(2,2)
fig.set_size_inches(8,6)

## plot AE and PS wavefunctions:
#ymin1,ymax1 = -0.0000015,0.000004
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
file1.close()
file2.close()
axes[0,0].plot(r20, r21, "blue")
axes[0,0].plot(r20, r22, "red")
axes[0,0].plot(r10, r11, "blue", linestyle='dashed')
axes[0,0].plot(r10, r12, "red", linestyle='dashed')
##axes[0].set_ylim(ymin1, ymax1)
(ymin1,ymax1) = axes[0,0].get_ylim()
axes[0,0].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
axes[0,0].set_ylabel('Wavefunction', size=12)
axes[0,0].set_xlabel('Radius (bohr)', size=12)
axes[0,0].legend(('PP_E=ae', 'PP_E=ref', 'AE_E=ae', 'AE_E=ref'), loc='best')
axes[0,0].vlines([rc], ymin1, ymax1, linestyle='dotted')
#axes[0,0].set_title("5d WF: ref = "+str(ref)+", rc = "+str(rc))
axes[0,0].set_title("5d: Wave Function")

# log derivative of AE wavefunctions
# and log derivative of PS wavefunctions
ymin2, ymax2 = -4, 4
rr10, rr11 = [],[]
rr20, rr21 = [],[]
ffile1 = open('D_chiae.dat', 'r')
ffile2 = open('D_chil.dat', 'r')
lines1 = ffile1.readlines()
lines2 = ffile2.readlines()
for line in lines1:
    rr10.append(float(line.split()[0]))
    rr11.append(float(line.split()[1]))
for line in lines2:
    rr20.append(float(line.split()[0]))
    rr21.append(float(line.split()[1]))
ffile1.close()
ffile2.close()
axes[0,1].plot(rr20, rr21, "blue")
axes[0,1].plot(rr10, rr11, "red", linestyle='dotted')
axes[0,1].set_ylim(ymin2, ymax2)
axes[0,1].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
axes[0,1].set_ylabel('Logarithmic derivative', size=12)
axes[0,1].set_xlabel('Energy (Ry)', size=12)
axes[0,1].legend(('PP', 'AE'), loc='best')
#axes[0,1].set_title("5D LG: ref = "+str(ref)+", rc = "+str(rc))
axes[0,1].set_title("5d: Logarithmic Derivative")

## plot Error in energy associated with E_cut:
r10, r11, r12 = [],[],[]
r20, r21, r22 = [],[],[]
xmin, xmax = 0, 50
ymin, ymax = 0, 0.002
axes[1,0].hlines([0.001], xmin, xmax, linestyle='dotted')
ffile1 = open('./D_delE.dat', 'r')
lines1 = ffile1.readlines()
for line in lines1:
    r10.append(float(line.split()[0]))
    r11.append(float(line.split()[1]))
ffile1.close()
axes[1,0].plot(r10, r11)
axes[1,0].set_xlim(xmin, xmax)
axes[1,0].set_ylim(ymin, ymax)
axes[1,0].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
axes[1,0].set_ylabel('Relative energy (Ry)', size=12)
axes[1,0].set_xlabel('Energy (Ry)', size=12)
axes[1,0].legend(('d orbital',), loc='upper right')
axes[1,0].set_title("5d: Error in Energy")

## plot Fourier components:
r10, r11, r12, r13 = [], [], [], []
ffile1 = open('./D_Qbar.dat', 'r')
lines1 = ffile1.readlines()
for line in lines1:
    r10.append(float(line.split()[0]))
    r11.append(float(line.split()[1]))
    r12.append(float(line.split()[2]))
    r13.append(float(line.split()[3]))
ffile1.close()
axes[1,1].plot(r10, r11)
axes[1,1].plot(r10, r12)
axes[1,1].plot(r10, r13)
axes[1,1].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
axes[1,1].set_ylabel('Fourier componets', size=12)
axes[1,1].set_xlabel('Energy (Ry)', size=12)
axes[1,1].legend(('L = 0', 'L = 2', 'L = 4'), loc='upper right')
axes[1,1].set_title("5d: Fourier Components")

fig.suptitle(f'5d orbital, ref = {ref}, rc = {rc}', fontsize=12)
plt.tight_layout()
plt.subplots_adjust(top=0.90)
##save="./5d_ref="+str(ref)+"_rc="+str(rc)+".eps"
##fig.savefig(save, format='eps')
save="./5d_ref="+str(ref)+"_rc="+str(rc)+".png"
fig.savefig(save, dpi=300)
