from pathlib import Path
import os
import sys
sys.path.append(str(Path(os.getcwd()).parent))



import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np



padf_path = Path(os.getcwd()).parent.parent /'py3padf02'/'padf'/'output'

qcorrel_path = Path(os.getcwd()).parent /'dbins'

qcorrel_fnames =   [
                    '1al1-sf_qcorrel',
                    '1mft-sf_qcorrel',
                    '1cos-sf_qcorrel',
                    # '4ggr-sf_qcorrel',
                    # '5z6y-sf_qcorrel',
                    # '4lqt-sf_qcorrel',
                    # '2b3p-sf_qcorrel',
                    # '2b3p-sf_highres_qcorrel',
                    # '5z6y-sf_highres_qcorrel',
                    ]


nQ = 256
nTheta =360
tmax = 180

rmax= 50*1e-10

r_space = np.linspace(0,rmax,nQ)
t_space = np.linspace(0,tmax, nTheta)

tscale, rscale = np.meshgrid(t_space, r_space)


for qcorrel_fname in qcorrel_fnames:



    rvol = ppu.read_dbin(str(padf_path/f'{qcorrel_fname}_padf.dbin'), nQ, nTheta)



    r1r2 = ppu.extract_r1r2(rvol*rscale**2)
    ppu.plot_map(r1r2[:,:180], title=f'{qcorrel_fname} r1r2',
                extent=[0,tmax,rmax,0], save=f'{qcorrel_fname}_r1r2')



    sumax0 =np.sum(rvol*rscale**2, axis=0)
    ppu.plot_map(sumax0[:,:180], title=f'{qcorrel_fname} sumax0',
                 extent=[0,tmax,rmax,0], save=f'{qcorrel_fname}_sumax0')



    rsx, rsy = np.meshgrid(r_space,r_space)
    sumax2 =np.sum(rvol, axis=2)*(rsx**(1.5))*(rsy**(1.5))
    ppu.plot_map(sumax2, title=f'{qcorrel_fname} sumax2',
                 extent=[0,rmax,rmax,0], save=f'{qcorrel_fname}_sumax2')




#   plt.figure()

#   theta_lin_sum = np.zeros(nTheta)

#   for i in range(0, nQ):
#       theta_lin_sum += r1r2[i,:]
#  # plt.plot(t_space, r1r2[i,:], label=f'$r_1={np.round(r_space[i])}')

#   plt.plot(t_space, theta_lin_sum)
#   plt.title(f'{qcorrel_fname[:4]} real space theta sum')




#vlines = [24,33, 42, 47,60, 72, 83, 90,120]  #1cos

#   vlines = [30,60, 45,100 ]  #1cos

#   for vline in vlines:

#       plt.axvline(linewidth=1, x=vline)
#       plt.text(vline,np.max(theta_lin_sum),f'{vline}',rotation=90)

#plt.legend()





#plt.figure()


#for i in range(0, int(nTheta/2), 18):
##   plt.plot(r_space, r1r2[:,i], label=f'$theta={np.round(t_space[i])}')

#plt.legend()


plt.show()







