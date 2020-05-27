from pathlib import Path
import plot_and_process_utils as ppu
import os
import matplotlib.pyplot as plt
import numpy as np





dbin_path = Path(str(os.getcwd())).parent.parent / 'dbins'

print(dbin_path.absolute())

dbin_fnames = [
    '1al1-sf_paper_qcorrel',
    '1cos-sf_paper_qcorrel',
#    '4osd-sf_paper_qcorrel'
]

qmaxs = [0.3699, 0.3699]


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


    


    r1r2 = ppu.extract_r1r2(qvol)


    r1r2 = r1r2**(0.25)
    
    r1r2 = ppu.convolve_gaussian(r1r2, 2,2)

    r1r2 -= np.min(r1r2)
    r1r2 /= np.max(r1r2)

    ppu.plot_map(r1r2,aspect='auto', title=f'{dbin_fname[:4].upper()} - Q space correlation',extent=[0,tmax,0,qmax], cmap='magma', xlabel='$\\theta$ / $ ^{\circ}$',ylabel='$q_1=q_2$ / $\AA^{-1}$', save=f'{dbin_fname[:4]}_q1q2' )

#    plt.plot(t_space, 3.5175*t_space**(-1.038291), 'r,')

    
    plt.xlim(0,180)
    plt.ylim(0, qmax)


plt.show()

