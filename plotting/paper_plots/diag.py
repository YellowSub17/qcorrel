import numpy as np
import os
import sys
from pathlib import Path
import plot_and_process_utils as ppu
import matplotlib.pyplot as plt




if __name__ == '__main__':
    

    nQ = 256
    nTheta = 360
    rmax=120
    nR = round(rmax/0.2)

    rcrop = 30
    nrcrop = round(rcrop/0.2)
    
    rmin =5
    nrmin = round(rmin/0.2)


    output_path = Path(os.getcwd()).parent.parent.parent/'padf'/'output'/'paper'
#    output_path  = Path('~/rmit-onedrive/phd/python_projects/padf/output/paper/')

#    fname = output_path / '1cos-sf_paper_qcorrel_wnl0blrr_0.dbin'
#    fname = output_path / '1al1-sf_paper_qcorrel_wnl0blrr_0.dbin'
    fname = output_path / '4osd-sf_paper_qcorrel_wnl0blrr_0.dbin'
    
    blrr_0 = ppu.read_blrr(str(fname), nR)

#    blrr_0 = np.clip(blrr_0, 0, np.max(blrr_0))

    blrr_0 = blrr_0[nrmin:nrcrop, nrmin:nrcrop]
    

    
    rspace = np.linspace(rmin,rcrop, nrcrop-nrmin)

    rx,ry = np.meshgrid(rspace, rspace)


    ppu.plot_map(blrr_0, extent=[rmin,rcrop, rmin, rcrop],title=f'{fname.stem[:4]} blrr_0 rmin {rmin}', xlabel='Position [A]', ylabel='Postion [A]',save=f'{fname.stem[:4]}_blrr_0_rmin{rmin}')

    pdf = np.zeros(nrcrop-nrmin)

    for i in range(nrcrop-nrmin):
        pdf[i] = blrr_0[i,i]

    plt.figure()
#    plt.plot(rspace,pdf*rspace**3)
    plt.plot(rspace,pdf)
    plt.title(f'{fname.stem[:4]} blrr_0 diag rmin {rmin}')
    plt.xlabel('Position [A]')
    plt.ylabel('Intensity')
    plt.savefig(f'{fname.stem[:4]}blrr_0_diag_rmin{rmin}.png')


    plt.show()
