import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import os
import time
from matplotlib import colors
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(figsize=(14,10))
diffuse = np.zeros([41,41])
for i in range(41):
    temp = "sqe_%.2f.txt" % (i/10.0)
    sqe = np.loadtxt(temp)
    diffuse[:,i] = np.sum(sqe[0:1251,:],axis=0)  # 0-25 meV
    #diffuse[:,i] = np.sum(sqe,axis=1)  # 0-25meV
np.savetxt('TDS_0_25meV.txt',diffuse)


#SQE = pl.log10(np.load('diffuse_HK1_0-2meV.npy'))
#norm=colors.Normalize(np.min(SQE),np.max(SQE))
vmin = 1.8
vmax = 3.8
tds = ax.imshow(np.log10(diffuse),aspect='auto',interpolation='bicubic',origin='lower',vmin=vmin,vmax=vmax,cmap='viridis')
#tds = ax.imshow(diffuse,aspect='auto',interpolation='bicubic',origin='lower',cmap='viridis')
plt.colorbar(tds)

## change xtick
xpos = np.linspace(0,40,5)
xl = np.linspace(0,4,5)
ax.set_xticks(xpos)
ax.set_xticklabels(xl)

## change ytick
ypos = np.linspace(0,40,5)
yl = np.linspace(0,4,5)
ax.set_yticks(ypos)
ax.set_yticklabels(yl)

ax.set_xlabel('H 0 0 (r.l.u)')
ax.set_ylabel('0 0 L (r.l.u)')
plt.title('NaCl thermal diffuse scattering, L=0')
plt.savefig('TDS_0-25meV.png')
