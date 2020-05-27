import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import plot_and_process_utils as ppu



##JACK NORMAL
j_data = np.loadtxt('1al1_ex_rcor.dat')
j_data = np.flip(j_data.T, axis=1)
j_data = np.clip(j_data, 0 , np.max(j_data))

t_space = np.linspace(0,180, j_data.shape[1])
r_space = np.linspace(0,20, j_data.shape[0])

rmax = 20
rmin = 2.5
res = j_data.shape[0]/rmax

nrmin = round(rmin*res)


ppu.plot_map(j_data[nrmin:,:], extent=[0,180,rmin,20], vmax=15, xlabel='$\\theta$ / $ ^{\circ}$', ylabel ='$r_1 = r_2$ /$\AA$', cmap='viridis')



arc_dists = [1.7, 2.7, 4.1,6]#, 10]
cols = plt.cm.autumn(np.linspace(0,1, len(arc_dists)))


for col, arc_dist in zip(cols, arc_dists):

    plt.plot(t_space, (arc_dist/2) /np.sin(np.radians(t_space/2)), '--',label=arc_dist, color=col)
    plt.plot(t_space, (arc_dist/2) /np.sin(np.radians((180-t_space)/2)), '--', color=col)

plt.xlim(0,180)
plt.ylim(rmin, 20)
plt.legend()


##JACK BLURRED
t_mesh, r_mesh = np.meshgrid(t_space, r_space)
j_data_blur = np.loadtxt('1al1_ex_rcor.dat')
j_data_blur = np.flip(j_data_blur.T, axis=1)


j_data_blur = np.clip(j_data_blur, 0 , np.max(j_data_blur))

j_data_blur = ppu.convolve_gaussian(j_data_blur, 4,4)


j_data_blur =j_data_blur/(r_mesh**2)
j_data_blur[0,:] =np.zeros((1, j_data_blur.shape[1]))






ppu.plot_map(j_data_blur[nrmin:,:],cmap='viridis',extent=[0,180,rmin,20], vmax=0.03)

xpnt = [65, 65, 72, 65, 90, 90, 27, 45, 90, 23, 50, 90, 90, 75, 55, 80]
ypnt = [8, 10, 11.8, 5.8, 4.5, 6.8, 5, 9, 9, 7, 7, 15, 13, 13.5, 12, 10.75]

cols = plt.cm.bone(np.linspace(0,1, len(xpnt)))
for x,y,col in zip(xpnt, ypnt, cols):
    plt.plot(x,y,'x', label=f'{x}, {y}',color = col )

plt.legend()






##PAT PADF
padf_fname = Path(os.getcwd()).parent.parent.parent/'padf'/ 'output'/'paper'/'1al1-sf_paper_qcorrel_wonl0_padf'

rmax = 60
cropr=20
res = 0.2
ncropr = round(cropr/res)
nrmax= round(rmax/res)
nrmin = round(rmin/res)

ntheta =360
tmax =360
cropt = 180

t_space = np.linspace(0,tmax, ntheta)
r_space = np.linspace(0,rmax, nrmax)

tscale, rscale = np.meshgrid(t_space, r_space)
rvol = ppu.read_dbin(f'{padf_fname}', nrmax, ntheta)

r1r2 = ppu.extract_r1r2(rvol)


r1r2 *= rscale**2
r1r2 = r1r2[nrmin:ncropr, :cropt]



r1r2 = np.clip(r1r2, 0, np.max(r1r2))
r1r2 = r1r2**(0.125)


r1r2 =r1r2 - np.min(r1r2)
r1r2 = r1r2/np.max(r1r2)

ppu.plot_map(r1r2,
      extent=[0, cropt,rmin,cropr], cmap='viridis', xlabel='$\\theta$ / $ ^{\circ}$', ylabel='$r_1 = r_2$ / $\\AA$')
arc_dists = [1.7, 2.7, 4.1,6]#, 10]
cols = plt.cm.autumn(np.linspace(0,1, len(arc_dists)))


#for col, arc_dist in zip(cols, arc_dists):
#
#    plt.plot(t_space, (arc_dist/2) /np.sin(np.radians(t_space/2)), '--',label=arc_dist, color=col)
#    plt.plot(t_space, (arc_dist/2) /np.sin(np.radians((180-t_space)/2)), '--', color=col)
#
cols = plt.cm.bone(np.linspace(0,1, len(xpnt)))
for x,y,col in zip(xpnt, ypnt, cols):
    plt.plot(x,y,'x', label=f'{x}, {y}',color = col )



plt.xlim(0,180)
plt.ylim(rmin, 20)
plt.legend()




plt.show()

