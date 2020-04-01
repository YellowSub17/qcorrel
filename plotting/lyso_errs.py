import plot_and_process_utils as ppu
from pathlib import Path
import os
import sys
import matplotlib.pyplot as plt
import numpy as np









nQ = 150
qmax=0.08
nTheta = 360
rmax =(nQ*1e-10)/(2*qmax)*1e9

r_space = np.linspace(0,rmax,nQ)
t_space = np.linspace(0,180, nTheta)

padf_path = Path(os.getcwd()).parent.parent /'py3padf02'/'padf'/'output'


#r1r2imshow
vol = ppu.read_dbin(str(padf_path/'254l-sf_res16_padf.dbin'), nQ, nTheta)
r1r2 = ppu.extract_r1r2(vol)


tscale, rscale = np.meshgrid(t_space, r_space)


ppu.plot_r1r2map(r1r2*rscale**2, title='254L - 16$\AA$ res cif (rscale **2)',
                 extent=[0,180,rmax,0],save='254lmap1')
ppu.plot_r1r2map(r1r2*rscale**4, title='254L - 16$\AA$ res cif (rscale **4)',
                 extent=[0,180,rmax,0],save='254lmap2')

padfs = [   '254l-sf_res16_padf.dbin',
            '254l_res16_err1_tag0_padf.dbin',
            '254l_res16_err1_tag1_padf.dbin',
            '254l_res16_err1_tag2_padf.dbin',
            '254l_res16_err1_tag3_padf.dbin',]


cols = ['r','b','b','b', 'b']

plt.figure()
plot_ind = 125
plt.title(f'254L - 16$\AA$ res cif\nPlot Through r1=r2={np.round(r_space[plot_ind],2)}$\AA$')
for i, padf in enumerate(padfs):
    vol = ppu.read_dbin(str(padf_path/padf), 150, 360)
    r1r2 = ppu.extract_r1r2(vol)
    plt.plot(np.linspace(0,180,360), r1r2[plot_ind,:]/np.max(r1r2[plot_ind,:]), f':{cols[i]}', label=padf[:20])

plt.legend()
plt.savefig('254lplot125.png')


plt.figure()
plot_ind = 20
plt.title(f'254L - 16$\AA$ res cif\nPlot Through r1=r2={np.round(r_space[plot_ind],2)}$\AA$')
for i, padf in enumerate(padfs):
    vol = ppu.read_dbin(str(padf_path/padf), 150, 360)
    r1r2 = ppu.extract_r1r2(vol)
    plt.plot(np.linspace(0,180,360), r1r2[plot_ind,:]/np.max(r1r2[plot_ind,:]), f':{cols[i]}', label=padf[:20])

plt.legend()
plt.savefig('254lplot20.png')






plt.show()



