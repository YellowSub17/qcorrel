from pathlib import Path
import os
import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np







padf_path = Path(os.getcwd()).parent.parent.parent/'padf'/ 'output'/'res_comp'

padf_fnames = [
    
    '1al1-sf_paper_qcorrel',
    '1cos-sf_paper_qcorrel',
]

res_names = [
    'res_lower',
    'res_equalqmax',
    'res_higher',

]
nRs  = [
    28,
    30,
    48
]

nQ = 256
nTheta =360
tmax =180
rmax=40




for res_name, nR in zip(res_names, nRs):
    for padf_fname in padf_fnames:
        t_space = np.linspace(0,tmax, nTheta)
        r_space = np.linspace(0,rmax, nR)

        tscale, rscale = np.meshgrid(t_space, r_space)
        print(str(padf_path/f'{padf_fname}_{res_name}_padf'), nR, nTheta)
        try:
            rvol = ppu.read_dbin(str(padf_path/f'{padf_fname}_{res_name}_padf'), nR, nTheta)
        except FileNotFoundError:
            print('File Not Found: Skipping')
            continue


#        rvol = rvol**(2)


        r1r2 = ppu.extract_r1r2(rvol*rscale**2)
        
        ppu.plot_map(r1r2[:,:180], title=f'{padf_fname} r1r2 rbins nR{nR}',
              extent=[0,tmax,0,rmax],vmin=0, cmap='gist_stern', save=f'{padf_fname}_r1r2_nR{nR}')

#
#
        sumax0 =np.sum(rvol*rscale**2, axis=0)
        
        ppu.plot_map(sumax0[:,:180], title=f'{padf_fname}  sumax0 nR{nR}',
              extent=[0,tmax,0,rmax],vmin=0, cmap='gist_stern', save=f'{padf_fname}_sumax0_nR{nR}')
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
#
plt.show()

