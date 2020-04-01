from pathlib import Path
import os
import sys
sys.path.append(str(Path(os.getcwd()).parent))



import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np



qcorrel_path = Path(os.getcwd()).parent /'dbins'/'convol_test'

qcorrel_fnames = ['253l-sf_qcorrel','254l-sf_qcorrel']

for i in range(4,25,4):
    qcorrel_fnames.append(f'253l-sf_qcorrel_filter{i}')
    qcorrel_fnames.append(f'254l-sf_qcorrel_filter{i}')


nQ = 150
qmax=0.14
nTheta =360
tmax =180


q_space = np.linspace(0, qmax, nQ)
t_space = np.linspace(0,tmax, nTheta)

tscale, qscale = np.meshgrid(t_space, q_space)
for qcorrel_fname in qcorrel_fnames:

    qvol = ppu.read_dbin(str(qcorrel_path/f'{qcorrel_fname}.dbin'), nQ, nTheta)

    r1r2 = ppu.extract_r1r2(qvol)

    #r1r2 = np.sum(qvol, axis=0)


#    r1r2 = ppu.convolve_gaussian(r1r2,3,3)
    ppu.plot_map(r1r2, title=f'{qcorrel_fname} (scaled)',
                    extent=[0,tmax,qmax,0], save=f'{qcorrel_fname}')



plt.show()

