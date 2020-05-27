import numpy as np
import matplotlib.pyplot as plt
import plot_and_process_utils as ppu




j_data = np.loadtxt('1al1_ex_rcor.dat')

j_data = np.flip(j_data.T, axis=1)

j_data = np.clip(j_data, 0, np.max(j_data))



j_data = j_data**(0.125)

#j_data = ppu.convolve_gaussian(j_data, 1,1)
ppu.plot_map(j_data, extent=[0, 180, 0, 20])

t_space = np.linspace(0,180,j_data.shape[1])

plt.plot(t_space, 1.4/np.sin(np.radians(t_space-1)),label=1.4)

plt.plot(t_space, 2.5/np.sin(np.radians(t_space-1)), label=2.5)
plt.plot(t_space, 4.8/np.sin(np.radians(t_space-1)), label=4.8)
plt.plot(t_space, 12.5/np.sin(np.radians(t_space-1)), label=12.5)
plt.plot(t_space, 15.4/np.sin(np.radians(t_space-1)), label=15.4)
plt.legend()
plt.xlim(0,180)
plt.ylim(0,20)

plt.show()
