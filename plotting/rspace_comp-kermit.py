from pathlib import Path
import os
import sys
sys.path.append(str(Path(os.getcwd()).parent))



import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np







padf_path = Path(os.getcwd()).parent.parent/'py3padf02'/'padf'/ 'output'

padf_fnames = [
    
    '1al1-sf_qcorrel',
    '1mft-sf_qcorrel',
    '1cos-sf_qcorrel',
    '4ggr-sf_qcorrel',
    '4lqt-sf_qcorrel'

]




nQ = 256
nTheta =360
tmax =180

rmaxs=[50,45,40,35,30,25,20]
res = round(256/50)

nRs = [round(rmax*res) for rmax in rmaxs]




for padf_fname in padf_fnames:
    for rmax, nR in zip(rmaxs, nRs):
        t_space = np.linspace(0,tmax, nTheta)
        r_space = np.linspace(0,rmax, nR)

        tscale, rscale = np.meshgrid(t_space, r_space)

        rvol = ppu.read_dbin(str(padf_path/f'{padf_fname}_rmax{rmax}_padf'), nR, nTheta)

        print(padf_fname, rvol.shape)



        r1r2 = ppu.extract_r1r2(rvol*rscale**2)
        ppu.plot_map(r1r2[:,:180], title=f'{padf_fname} rmax{rmax} r1r2',
                         extent=[0,tmax,rmax,0],
                         save=f'{padf_fname}_rmax{rmax}_r1r2')



        sumax0 =np.sum(rvol*rscale**2, axis=0)
        ppu.plot_map(sumax0[:,:180], title=f'{padf_fname} rmax{rmax} sumax0',
                     extent=[0,tmax,rmax,0], save=f'{padf_fname}_rmax{rmax}_sumax0')



        rsx, rsy = np.meshgrid(r_space,r_space)
        sumax2 =np.sum(rvol, axis=2)*(rsx**(1.5))*(rsy**(1.5))
        ppu.plot_map(sumax2, title=f'{padf_fname} rmax{rmax} sumax2',
                     extent=[0,rmax,rmax,0], save=f'{padf_fname}_rmax{rmax}_sumax2')

        plt.close('all')

plt.show()

