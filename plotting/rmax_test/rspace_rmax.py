from pathlib import Path
import os
import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np







padf_path = Path(os.getcwd()).parent.parent.parent/'padf'/ 'output'/ 'rmaxtest'

padf_fnames = [
    
    '1al1-sf_rmaxtest_qcorrel_rmax',
#    '1cos-sf_rmaxtest_qcorrel_rmax',
]

rmaxs = [40, 50, 60]

res = 0.2

nRs= [round (rmax/res) for rmax in rmaxs]


nQ = 256
nTheta =360
tmax =180

multi = np.ones((nRs[0]-1, 180))

ave = np.zeros((nRs[0]-1, 180))


#####
for rmax, nR in zip(rmaxs, nRs):
    for padf_fname in padf_fnames:
        t_space = np.linspace(0,tmax, nTheta)
        r_space = np.linspace(0,rmax, nR)

        tscale, rscale = np.meshgrid(t_space, r_space)
        try:
            rvol = ppu.read_dbin(str(padf_path/f'{padf_fname}{rmax}_padf'), nR, nTheta)
        except FileNotFoundError:
            print('File Not Found: Skipping')
            continue





        r1r2 = ppu.extract_r1r2(rvol)*rscale**2
        
        r1r2 = np.clip(r1r2, 0, np.max(r1r2))**(0.125)


        ppu.plot_map(r1r2[:nRs[0]-1,:180], title=f'{padf_fname}{rmax}  r1r2',
              extent=[0,tmax,0,r_space[nRs[0]-1]], cmap='viridis', save=f'{padf_fname}{rmax}r1r2')

        multi *= r1r2[:nRs[0]-1,:180]

        ave += r1r2[:nRs[0]-1,:180]/3
#
#
#        sumax0 =np.sum(rvol*rscale**2, axis=0)
#        
#        ppu.plot_map(sumax0[:,:180], title=f'{padf_fname}  sumax0',
#              extent=[0,tmax,0,rmax], cmap='gist_stern', vmin=0)
#
#
#
#        rsx, rsy = np.meshgrid(r_space,r_space)
#        
#        sumax2 =np.sum(rvol, axis=2)*(rsx**(1.5))*(rsy**(1.5))
#        
#        ppu.plot_map(sumax2, title=f'{padf_fname}  sumax2', cmap='gist_stern', 
#              extent=[0,rmax,0,rmax],vmin=0)
#
#i

ppu.plot_map(multi,extent=[0,tmax,0,r_space[nRs[0]-1]], cmap='viridis', save=f'{padf_fname}{rmax}r1r2_multi', title=f'{padf_fname}_multi')

ppu.plot_map(ave,extent=[0,tmax,0,r_space[nRs[0]-1]], cmap='viridis', save=f'{padf_fname}{rmax}r1r2_ave', title=f'{padf_fname}_ave')
plt.show()

