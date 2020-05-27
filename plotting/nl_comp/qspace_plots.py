from pathlib import Path
import plot_and_process_utils as ppu
import os
import matplotlib.pyplot as plt
import numpy as np





dbin_path = Path(str(os.getcwd())).parent.parent / 'dbins'

print(dbin_path.absolute())

dbin_fnames = [
    '1al1-sf_highres_qcorrel',
    '1cos-sf_highres_qcorrel',
]

qmaxs = [0.3699, 0.5]


nQ = 256
nTheta =360

tmax =180



for dbin_fname, qmax in zip(dbin_fnames, qmaxs):
    t_space = np.linspace(0,tmax, nTheta)
    q_space = np.linspace(0,qmax, nQ)

    tscale, qscale = np.meshgrid(t_space, q_space)
    try:
        qvol = ppu.read_dbin(str(dbin_path/f'{dbin_fname}'), nQ, nTheta)
    except FileNotFoundError:
        print(f'File {dbin_fname} Not Found: Skipping')
        continue



    r1r2 = ppu.extract_r1r2(qvol)#*qscale**2)

    r1r2 = ppu.convolve_gaussian(r1r2)
    
    ppu.plot_map(r1r2, title=f'{dbin_fname} r1r2',
          extent=[0,tmax,0,qmax], cmap='viridis', xlabel='Correlation Angle $\\theta$ [degrees]',ylabel='Correlation length $q$ [1/$\AA$]', save=f'{dbin_fname}_r1r2' )



#    sumax0 =np.sum(qvol*qscale**2, axis=0)
#    
#    ppu.plot_map(sumax0, title=f'{dbin_fname} sumax0',
#          extent=[0,tmax,qmax,0], cmap='gist_stern', vmin=0)
#
#
#
#    qsx, qsy = np.meshgrid(q_space,q_space)
#    
#    sumax2 =np.sum(qvol, axis=2)*(qsx**2)*(qsy**2)
#    
#    ppu.plot_map(sumax2, title=f'{dbin_fname} sumax2', cmap='gist_stern', 
#          extent=[0,qmax,qmax,0],vmin=0)
#

plt.show()

