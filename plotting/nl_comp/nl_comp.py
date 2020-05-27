from pathlib import Path
import os

#append Qcorrel directory to get plot_and_process_utils
import plot_and_process_utils as ppu


import matplotlib.pyplot as plt
import numpy as np







padf_path = Path(os.getcwd()).parent.parent.parent/'py3padf02'/'padf'/ 'output'/'nl_comp'

padf_fnames = [
    
    '1al1-sf_qcorrel',
    '1mft-sf_qcorrel',
    '1cos-sf_qcorrel',


]




nQ = 256
nTheta =360
tmax =180

#nls=[36, 32,28,24, 20, 16,12]
nls=[56,60, 64]
res = 2 #pix/A
rmaxs = [30, 40, 50, 60, 70,80,90,100,110,120]

nRs = [round(rmax*res) for rmax in rmaxs]



for rmax, nR in zip(rmaxs, nRs):
    for padf_fname in padf_fnames:
        for nl in nls:
            t_space = np.linspace(0,tmax, nTheta)
            r_space = np.linspace(0,rmax, nR)

            tscale, rscale = np.meshgrid(t_space, r_space)
            try:
                rvol = ppu.read_dbin(str(padf_path/f'{padf_fname}_rmax{rmax}_nl{nl}_padf'), nR, nTheta)
            except FileNotFoundError:
                print('File Not Found: Skipping')
                continue


            print(padf_fname, rvol.shape)



            r1r2 = ppu.extract_r1r2(rvol*rscale**2)
            
            ppu.plot_map(r1r2[:,:180], title=f'{padf_fname} nl {nl} r1r2',
                  extent=[0,tmax,rmax,0],vmin=0, cmap='gist_stern', save=f'{padf_fname}_rmax{rmax}_nl{nl}_r1r2')



            sumax0 =np.sum(rvol*rscale**2, axis=0)
            
            ppu.plot_map(sumax0[:,:180], title=f'{padf_fname} nl {nl} sumax0',
                  extent=[0,tmax,rmax,0], cmap='gist_stern', vmin=0, save=f'{padf_fname}_rmax{rmax}_nl{nl}_sumax0')



            rsx, rsy = np.meshgrid(r_space,r_space)
            
            sumax2 =np.sum(rvol, axis=2)*(rsx**(1.5))*(rsy**(1.5))
            
            ppu.plot_map(sumax2, title=f'{padf_fname} nl {nl} sumax2', cmap='gist_stern', 
                  extent=[0,rmax,rmax,0],vmin=0, save=f'{padf_fname}_rmax{rmax}_nl{nl}_sumax2')

            plt.close('all')

plt.show()

