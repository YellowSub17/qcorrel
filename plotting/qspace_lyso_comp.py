from pathlib import Path
import os
import sys
sys.path.append(str(Path(os.getcwd()).parent))



import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np



qcorrel_path = Path(os.getcwd()).parent /'dbins'

qcorrel_fnames = ['253l-sf_ave_qcorrel','254l-sf_ave_qcorrel', '253l-sf_ave_comp254l-sf_ave_qcorrel','254l-sf_ave_comp253l-sf_ave_qcorrel']



nQ = 256
qmax=0.14
nTheta =360
tmax =180


q_space = np.linspace(0, qmax, nQ)
t_space = np.linspace(0,tmax, nTheta)

tscale, qscale = np.meshgrid(t_space, q_space)


for qcorrel_fname in qcorrel_fnames:

    qvol = ppu.read_dbin(str(qcorrel_path/f'{qcorrel_fname}.dbin'), nQ, nTheta)

    r1r2 = ppu.extract_r1r2(qvol)
    r1r2 = ppu.convolve_gaussian(r1r2,2,2)
    ppu.plot_map(r1r2, title=f'{qcorrel_fname} r1r2',
                    extent=[0,tmax,qmax,0], save=f'{qcorrel_fname}_r1r2')


    sumax0 =np.sum(qvol, axis=0)
    sumax0 = ppu.convolve_gaussian(sumax0,2,2)
    ppu.plot_map(sumax0, title=f'{qcorrel_fname} sumax0',
                 extent=[0,tmax,qmax,0], save=f'{qcorrel_fname}_sumax0')

    sumax2 =np.sum(qvol, axis=2)
    sumax2 = ppu.convolve_gaussian(sumax2,2,2)
    ppu.plot_map(sumax2, title=f'{qcorrel_fname} sumax2',
                 extent=[0,qmax,qmax,0], save=f'{qcorrel_fname}_sumax2')



plt.show()

