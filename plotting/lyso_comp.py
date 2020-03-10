import plot_and_process_utils as ppu
from pathlib import Path
import os
import sys
import matplotlib.pyplot as plt
import numpy as np



nQ = 150
qmax=0.22
nTheta = 180
rmax =(nQ*1e-10)/(2*qmax)*1e9

r_space = np.linspace(0,rmax,nQ)
t_space = np.linspace(0,180, nTheta)

padf_path = Path(os.getcwd()).parent.parent /'py3padf02'/'padf'/'output'



nQ = 150
qmax=0.22
nTheta = 180
rmax =(nQ*1e-10)/(2*qmax)*1e9

r_space = np.linspace(0,rmax,nQ)
t_space = np.linspace(0,180, nTheta)
tscale, rscale = np.meshgrid(t_space, r_space)


padfs = [
            '254l_res8_err1_tag0_padf.dbin',
            '253l_res8_err1_tag0_padf.dbin'
            ]

for padf in padfs:
#r1r2imshow
    vol = ppu.read_dbin(str(padf_path/padf), nQ, nTheta)
    r1r2 = ppu.extract_r1r2(vol)
    ppu.plot_r1r2map(r1r2*rscale**2, title=f'{padf}(rscale **2)',
                     extent=[0,180,rmax,0],save=f'{padf[:4]}_res8rsc2')
    ppu.plot_r1r2map(r1r2*rscale**4, title=f'{padf}(rscale **4)',
                     extent=[0,180,rmax,0],save=f'{padf[:4]}_res8rsc4')



cols = ['r','b','g','y']


plot_ind = 10
plt.figure()
for i, padf in enumerate(padfs):

    vol = ppu.read_dbin(str(padf_path/padf), nQ, nTheta)
    r1r2 = ppu.extract_r1r2(vol)

    plt.plot(np.linspace(0,180,nTheta), r1r2[plot_ind,:]/np.max(r1r2[plot_ind,:]), c=cols[i], label=padf)

plt.title(f'254L/253L comparision \nPlot Through r1=r2={np.round(r_space[plot_ind],2)}$\AA$')
plt.legend()
plt.savefig('253l254l_comp_10.png')



plot_ind = 100
plt.figure()
for i, padf in enumerate(padfs):

    vol = ppu.read_dbin(str(padf_path/padf), nQ, nTheta)
    r1r2 = ppu.extract_r1r2(vol)

    plt.plot(np.linspace(0,180,nTheta), r1r2[plot_ind,:]/np.max(r1r2[plot_ind,:]), c=cols[i], label=padf)

plt.title(f'254L/253L comparision \nPlot Through r1=r2={np.round(r_space[plot_ind],2)}$\AA$')
plt.legend()
plt.savefig('253l254l_comp_100.png')




plot_ind = 50
plt.figure()
for i, padf in enumerate(padfs):

    vol = ppu.read_dbin(str(padf_path/padf), nQ, nTheta)
    r1r2 = ppu.extract_r1r2(vol)

    plt.plot(np.linspace(0,180,nTheta), r1r2[plot_ind,:]/np.max(r1r2[plot_ind,:]), c=cols[i], label=padf)

plt.title(f'254L/253L comparision \nPlot Through r1=r2={np.round(r_space[plot_ind],2)}$\AA$')
plt.legend()
plt.savefig('253l254l_comp_50.png')
