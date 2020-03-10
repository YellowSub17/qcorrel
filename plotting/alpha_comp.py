import plot_and_process_utils as ppu
from pathlib import Path
import os
import sys
import matplotlib.pyplot as plt
import numpy as np



padf_path = Path(os.getcwd()).parent.parent /'py3padf02'/'padf'/'output'



nQ = 150
qmax=0.2
nTheta = 180
rmax =(nQ*1e-10)/(2*qmax)*1e9

r_space = np.linspace(0,rmax,nQ)
t_space = np.linspace(0,180, nTheta)


#r1r2imshow
vol = ppu.read_dbin(str(padf_path/'4zry-sf_res8_qcorrel_padf.dbin'), nQ, nTheta)
r1r2 = ppu.extract_r1r2(vol)


tscale, rscale = np.meshgrid(t_space, r_space)


ppu.plot_r1r2map(r1r2*rscale**2, title='4ZRY - 8$\AA$ res cif (rscale **2)',
                 extent=[0,180,rmax,0], save='4zry_maprs2')
ppu.plot_r1r2map(r1r2*rscale**3, title='4ZRY - 8$\AA$ res cif (rscale **3)',
                 extent=[0,180,rmax,0], save='4zry_maprs3')


#
#
# padfs = [  '4zry-sf_res8_qcorrel_padf.dbin',
# ]
#
#
# cols = ['r','m','g','y']
#
#
# plot_ind = 50
# plt.figure()
# for i, padf in enumerate(padfs):
#
#     vol = ppu.read_dbin(str(padf_path/padf), nQ, nTheta)
#     r1r2 = ppu.extract_r1r2(vol)
#     plt.plot(np.linspace(0,180,nTheta), r1r2[plot_ind,:]/np.max(r1r2[plot_ind,:]), c=cols[i])
#
#






nQ = 150
qmax=0.36
nTheta = 180
rmax =(nQ*1e-10)/(2*qmax)*1e9

r_space = np.linspace(0,rmax,nQ)
t_space = np.linspace(0,180, nTheta)


#r1r2imshow
vol = ppu.read_dbin(str(padf_path/'1cos-sf_res4_qcorrel_padf.dbin'), nQ, nTheta)
r1r2 = ppu.extract_r1r2(vol)


tscale, rscale = np.meshgrid(t_space, r_space)


ppu.plot_r1r2map(r1r2*rscale**2, title='1COS - 4$\AA$ res cif (rscale **2)',
                 extent=[0,180,rmax,0], save='1cos_maprs2')


ppu.plot_r1r2map(r1r2*rscale**3, title='1COS - 4$\AA$ res cif (rscale **3)',
                 extent=[0,180,rmax,0], save='1cos_maprs3')




# padfs = [  '1cos-sf_res4_qcorrel_padf.dbin',
#            ]
#
#
# cols = ['r','m','g','y']
#
#
# plot_ind = 50
# plt.figure()
# for i, padf in enumerate(padfs):
#
#     vol = ppu.read_dbin(str(padf_path/padf), nQ, nTheta)
#     r1r2 = ppu.extract_r1r2(vol)
#     plt.plot(np.linspace(0,180,nTheta), r1r2[plot_ind,:]/np.max(r1r2[plot_ind,:]), c=cols[i])
#
