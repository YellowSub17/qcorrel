from pathlib import Path
import plot_and_process_utils as ppu
import os
import matplotlib.pyplot as plt
import numpy as np


import matplotlib.pylab as pylab
print(pylab.rcParams.keys())
params = { 'axes.labelsize':16,
         'axes.titlesize':16,
        'figure.titlesize':20,
          'xtick.labelsize':14,
         'ytick.labelsize':14}
pylab.rcParams.update(params)


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

    # ppu.plot_map(r1r2,aspect='auto', title=f'{dbin_fname[:4].upper()} - $q$ Space Correlation',
                 # extent=[0,tmax,0,qmax], xlabel='$\\theta$ / $ ^{\circ}$',fig_size=(8,5.5), 
                 # ylabel='$q_1=q_2$ / $\AA^{-1}$', save=f'{dbin_fname[:4]}_q1q2_rv' , cmap='viridis')
    plt.figure(figsize=(8,5.5), dpi=100)
    plt.suptitle(f'{dbin_fname[:4].upper()} - $q$ Space Correlation')
    ax1 = plt.subplot(2,2,3)

    plt.imshow(r1r2, cmap='viridis', aspect='auto', extent=[0,tmax,0,qmax], origin='lower')
    plt.xlim(0,180)
    plt.ylim(0, qmax)

    # plt.title(f'{dbin_fname[:4].upper()} - $q$ Space Correlation')
    plt.xlabel('$\\theta$ / $ ^{\circ}$')
    plt.ylabel('$q_1=q_2$ / $\AA^{-1}$')

    ax2 = plt.subplot(2,2,4)
    # ax2.set_title('Mean Inten.($q$)')
    plt.plot(np.average(r1r2,axis=1), q_space)
    plt.ylim([0, qmax])
    plt.xlim([0, np.max(np.average(r1r2,axis=1))])
    # plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()

    # plt.xlabel('$q_1=q_2$ / $\AA^{-1}$')
    plt.xticks([])
    plt.yticks([])
    plt.xlabel('Mean Inten.($q$)')


    ax2 = plt.subplot(2,2,1)
    # ax2.set_title('Mean Inten.($\\theta$)')
    plt.plot(t_space, np.average(r1r2,axis=0))
    plt.xlim([0, 180])
    plt.ylim([0, np.max(np.average(r1r2,axis=0))])

    # plt.xlabel('$\\theta$ / $ ^{\circ}$')
    plt.xticks([])
    plt.yticks([])
    plt.subplots_adjust(wspace=0,
                    hspace=0.0)
    plt.ylabel('Mean Inten.($\\theta$)')



    # plt.colorbar()
    # plt.savefig(f'{dbin_fname[:4]}_q1q2_rv' )


#    plt.plot(t_space, 3.5175*t_space**(-1.038291), 'r,')

    

plt.show()

