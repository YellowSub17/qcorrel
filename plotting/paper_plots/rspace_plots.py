from pathlib import Path
import os
import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np







padf_path = Path(os.getcwd()).parent.parent.parent/'padf'/ 'output'/'paper'


padf_fnames = [
    
#    '1al1-sf_paper_qcorrel_wnl0_padf',
    '1al1-sf_paper_qcorrel_wonl0_padf',
#    '1cos-sf_paper_qcorrel_wnl0_padf',
#    '1cos-sf_paper_qcorrel_wonl0_padf',
#     '4osd-sf_paper_qcorrel_wonl0_padf',
]

#rmaxs = [60,60,120]
rmaxs = [60,60]
cropr=25
rmin =9 

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
#    r1r2 *= rscale**2
    r1r2 = r1r2[nrmin:ncropr, :cropt]
    


    r1r2 = np.clip(r1r2, 0, np.max(r1r2))
    


    r1r2 =r1r2 - np.min(r1r2)
    r1r2 = r1r2/np.max(r1r2)
    
    ppu.plot_map(r1r2, title=f'{padf_fname}  r1r2',
          extent=[0, cropt,rmin,cropr], cmap='viridis', save=f'{padf_fname}_r1r2', xlabel='Correlation Angle $\\theta$ [degrees]', ylabel='Correlation length $r$ [$\\AA$]')

    iso_peaks = [16.6, 10.1, 16.2, 15.2, 10.3, 9.8, 10.3, 15.4, 17.3, 20.4, 17.8, 11.7, 19.0 ]

#    colors = plt.cm.jet(np.linspace(0,1,len(theta_lins)))
#
#    
#    plt.figure()
#    plt.title(f'{padf_fname}_r1r2_theta_lin')
#    print(len(colors), len(theta_lins))
#    for col, theta_lin in zip(colors, theta_lins):
#        theta_lin_plot = r1r2[:,theta_lin]
#        plt.plot(r_space[nrmin:ncropr].T,theta_lin_plot,color=col, label=str(theta_lin))
#  #16 A peak falls on arc of negative correlation  
#    iso_peaks = [16.6, 10.1, 16.2, 15.2, 10.3, 9.8, 10.3, 15.4, 17.3, 20.4, 17.8, 11.7, 19.0 ]
##    iso_peaks = []
#    for peak in iso_peaks:
#
#        plt.axvline(peak)
#
#    plt.legend()
#

    plt.figure()
    plt.title(f'{padf_fname}_r1r2_ave_theta_{theta_lins[-1]}-{theta_lins[0]}')
    for peak in iso_peaks:

        plt.axvline(peak)
    theta_ave_plot = np.mean(r1r2[:,theta_lins[0]:theta_lins[-1]], axis=-1)
    plt.plot(r_space[nrmin:ncropr].T,theta_ave_plot,color='red')
    plt.savefig(f'{padf_fname}_r1r2_ave_theta_{theta_lins[-1]}-{theta_lins[0]}')    

plt.show()

