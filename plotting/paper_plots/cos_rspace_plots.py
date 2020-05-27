from pathlib import Path
import os
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
theta_lins = [i for i in range(55,66, 1)]



#####
for rmax, padf_fname in zip(rmaxs,padf_fnames):
    nrmax= round(rmax/res)
    t_space = np.linspace(0,tmax, ntheta)
    r_space = np.linspace(0,rmax, nrmax)

    tscale, rscale = np.meshgrid(t_space, r_space)
    try:
        rvol = ppu.read_dbin(str(padf_path/f'{padf_fname}'), nrmax, ntheta)
    except FileNotFoundError:
        print(f'File Not Found: {padf_fname}')
        continue





    r1r2 = ppu.extract_r1r2(rvol)
    np.savetxt('1cos_padf_r1r2.dat', r1r2, header=f'{r1r2.shape}')

#    r1r2 *= rscale**2
    r1r2 = r1r2[nrmin:ncropr, :cropt]
    


    r1r2 = np.clip(r1r2, 0, np.max(r1r2))

    r1r2 = r1r2**(0.125)


    r1r2 =r1r2 - np.min(r1r2)
    r1r2 = r1r2/np.max(r1r2)
    
    ppu.plot_map(r1r2, title=f'{padf_fname[:4].upper()} - PADF',
          extent=[0, cropt,rmin,cropr], cmap='Blues', save=f'{padf_fname[:4]}_r1r2', xlabel='$\\theta$ / $ ^{\circ}$', ylabel='$r_1 = r_2$ / $\\AA$')

    dists = [2.2, 6]
    cols = ['purple', 'red']
    for col, dist in zip(cols,dists):
        plt.plot(t_space, (dist/2)/np.sin(np.radians(t_space)/2),color=col, label=f'{dist} $\AA$')
        plt.plot(t_space, (dist/2)/np.sin(np.radians(180-t_space)/2),color=col )
    plt.xlim(0,180)
    plt.ylim(rmin, cropr)
    plt.legend(loc='upper right', ncol=1)
    plt.savefig(f'{padf_fname[:4]}_r1r2')

#    iso_peaks = [16.6, 10.1, 16.2, 15.2, 10.3, 9.8, 10.3, 15.4, 17.3, 20.4, 17.8, 11.7, 19.0 ]
#
#    plt.figure()
#    plt.title(f'{padf_fname}_r1r2_ave_theta_{theta_lins[-1]}-{theta_lins[0]}')
#    for peak in iso_peaks:
#
#        plt.axvline(peak)
#    theta_ave_plot = np.mean(r1r2[:,theta_lins[0]:theta_lins[-1]], axis=-1)
#    plt.plot(r_space[nrmin:ncropr].T,theta_ave_plot,color='red')
#    plt.savefig(f'{padf_fname}_r1r2_ave_theta_{theta_lins[-1]}-{theta_lins[0]}')    

plt.show()

