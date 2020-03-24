from pathlib import Path
import os
import sys
sys.path.append(str(Path(os.getcwd()).parent))



import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np







padf_path = Path(os.getcwd()).parent.parent/'py3padf02'/'padf'/ 'output'

padf_fnames = ['253l-sf_ave_qcorrel_padf','254l-sf_ave_qcorrel_padf', '253l-sf_ave_comp254l-sf_ave_qcorrel_padf','254l-sf_ave_comp253l-sf_ave_qcorrel_padf']



nQ = 256
qmax=0.14
nTheta =360
tmax =180
rmax = 50e-10


q_space = np.linspace(0, qmax, nQ)
t_space = np.linspace(0,tmax, nTheta)
r_space = np.linspace(0,rmax, nQ)

tscale, rscale = np.meshgrid(t_space, r_space)

for padf_fname in padf_fnames:

    rvol = ppu.read_dbin(str(padf_path/f'{padf_fname}.dbin'), nQ, nTheta)



    r1r2 = ppu.extract_r1r2(rvol*rscale**2)
    ppu.plot_map(r1r2[:,:180], title=f'{padf_fname} r1r2',
                 extent=[0,tmax,rmax,0], save=f'{padf_fname}_r1r2')



    sumax0 =np.sum(rvol*rscale**2, axis=0)
    ppu.plot_map(sumax0[:,:180], title=f'{padf_fname} sumax0',
                 extent=[0,tmax,rmax,0], save=f'{padf_fname}_sumax0')



    rsx, rsy = np.meshgrid(r_space,r_space)
    sumax2 =np.sum(rvol, axis=2)*(rsx**(1.5))*(rsy**(1.5))
    ppu.plot_map(sumax2, title=f'{padf_fname} sumax2',
                 extent=[0,rmax,rmax,0], save=f'{padf_fname}_sumax2')


plt.show()

