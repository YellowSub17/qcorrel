import os
from pathlib import Path
import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np


padf_path = Path(os.getcwd()).parent.parent.parent/'padf'/ 'output'/'paper'


padf_fnames = [
    
    '1cos-sf_paper_qcorrel_wonl0_padf',
]

rmaxs = [60]
cropr=20
rmin =0 

res = 0.2

ncropr = round(cropr/res)
nrmin = round(rmin/res)

nQ = 256
ntheta =360
tmax =360
cropt = 180
theta_lins = [0, 20, 40, 60, 90]





#####
for rmax, padf_fname in zip(rmaxs,padf_fnames):
    nrmax= round(rmax/res)
    t_space = np.linspace(0,tmax, ntheta)
    r_space = np.linspace(rmin,cropr, ncropr-nrmin)


    tscale, rscale = np.meshgrid(t_space, r_space)
    rxx, ryy = np.meshgrid(r_space, r_space)


    try:
        rvol = ppu.read_dbin(str(padf_path/f'{padf_fname}'), nrmax, ntheta)
    except FileNotFoundError:
        print(f'File Not Found: {padf_fname}')
        continue

    rvol = rvol[nrmin:ncropr,nrmin:ncropr, :cropt]
    print('rvol', rvol.shape)
    print('rxx', rxx.shape)
    print('ryy', ryy.shape)
    
    
    for theta in theta_lins:

        rrtheta = rvol[:,:,theta]
        rrtheta = np.clip(rrtheta, 0, np.max(rrtheta))
        print(rrtheta.shape)
        rrtheta = rrtheta*(rxx**1)
        rrtheta = rrtheta*(ryy**1)
        ppu.plot_map(rrtheta, extent=[0,cropr,0,cropr],
                     title=f'{padf_fname[:4].upper()} - $\\theta$={theta}$^\circ$',
                     save=f'{padf_fname[:4]}_rrtheta{theta}',cmap='viridis',
                     xlabel='$r$ / $\\AA$', ylabel='$r`$ / $\\AA$')


plt.show()
