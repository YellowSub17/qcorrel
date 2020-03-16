from pathlib import Path
import os
import sys
sys.path.append(str(Path(os.getcwd()).parent))



import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
import numpy as np



padf_path = Path(os.getcwd()).parent.parent /'py3padf02'/'padf'/'output'

qcorrel_path = Path(os.getcwd()).parent /'dbins'

qcorrel_fnames =[ '1cos-sf_qcorrelnl20', '1mft-sf_qcorrelnl20']





nQ = 150
qmax=0.3
nTheta =360
tmax = 360

#rmax =(nQ*1e-10)/(2*qmax)
rmax= 50*1e-10

q_space = np.linspace(0, qmax, nQ)
r_space = np.linspace(0,rmax,nQ)
t_space = np.linspace(0,tmax, nTheta)

tscale, rscale = np.meshgrid(t_space, r_space)


for qcorrel_fname in qcorrel_fnames:
## Real Space
#r1r2imshow
    rvol = ppu.read_dbin(str(padf_path/f'{qcorrel_fname}_padf.dbin'), nQ, nTheta)

    rvol = rvol*rvol
    #r1r2 = ppu.extract_r1r2(rvol)*rscale**3
    r1r2 = np.sum(rvol, axis=0)*rscale**2




    ppu.plot_r1r2map(r1r2[:80,:], title=f'{qcorrel_fname[:4]} padf (scaled)',
                extent=[0,tmax,rmax,0])



    plt.figure()

    theta_lin_sum = np.zeros(nTheta)

    for i in range(0, nQ):
        theta_lin_sum += r1r2[i,:]
   # plt.plot(t_space, r1r2[i,:], label=f'$r_1={np.round(r_space[i])}')

    plt.plot(t_space, theta_lin_sum)
    plt.title(f'{qcorrel_fname[:4]} real space theta sum')




#vlines = [24,33, 42, 47,60, 72, 83, 90,120]  #1cos

    vlines = [30,60, 45,100 ]  #1cos

    for vline in vlines:

        plt.axvline(linewidth=1, x=vline)
        plt.text(vline,np.max(theta_lin_sum),f'{vline}',rotation=90)

#plt.legend()





#plt.figure()


#for i in range(0, int(nTheta/2), 18):
 #   plt.plot(r_space, r1r2[:,i], label=f'$theta={np.round(t_space[i])}')

#plt.legend()


plt.show()







