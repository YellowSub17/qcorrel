import matplotlib.pyplot as plt
import numpy as np
import imageio
from pathlib import Path

import struct



def plot_map(map, title='', save=None, cmap='plasma', extent=None,
        xlabel = 'Correlation Angle $\Delta$ [Degrees]', ylabel= 'Scattering Magnitude $q$ [1/$\AA$]'):

    plt.figure()
    plt.title(title)
    plt.imshow(map, cmap=cmap, aspect='auto', extent=extent)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar()
    if type(save)==str:
        save_path = Path('saved_plots') / f'{save}.png'
        plt.savefig(str(save_path))
        plt.close(plt.gcf().number)

def save_tiff(im, name="IMAGE"):

    if name[-4:]=='tiff':
        save_path = Path('tiffs') / f'{name}'
    else:
        save_path = Path('tiffs') / f'{name}.tiff'
    imageio.imsave(save_path, im)


def save_dbin(data, name='DBIN'):

    if name[-4:]=='dbin':
        save_path = Path('dbins') / f'{name}'
    else:
        save_path = Path('dbins') / f'{name}.dbin'
    
    f = open(save_path, "wb")
    fmt = '<' + 'd' * data.size
    bin_in = struct.pack(fmt, *data.flatten()[:])
    f.write(bin_in)
    f.close()

def read_dbin(path, nQ, nTheta):

    if path[-4:]=='dbin':
        save_path = Path( f'{path}')
    else:
        save_path = Path(f'{path}.dbin')
    arr = np.fromfile(str(save_path))
    arr = arr.reshape(nQ,nQ,nTheta)

    return arr




def array_shift(array, xshift=0, yshift=0):
    ###### Gaussian Convolution Functions by AM
    # shift - a 2D version of numpy's roll
    array = np.roll(array, xshift, 0)
    array = np.roll(array, yshift, 1)
    return array


def make_gaussian(nx, ny, rad=None, rady=-1., cenx=None, ceny=None, invert=0, norm=False):
    ###### Gaussian Convolution Functions by AM
    ## make a 2D array with a gaussian
    # set defaults
    if rad is None: rad = np.min(nx, ny) / 2
    if cenx is None: cenx = nx / 2
    if ceny is None: ceny = ny / 2
    radsq = rad ** 2
    if rady == -1.:
        radysq = radsq
    else:
        radysq = rady ** 2

    # define the circle
    x = np.outer(np.arange(0 - nx / 2, nx - nx / 2, 1), np.ones(ny))
    # print x.size, x.shape
    y = np.outer(np.ones(nx), np.arange(0 - ny / 2, ny - ny / 2, 1))
    # print y.size, y.shape

    a = np.zeros([nx, ny])
    a = np.exp(-(x ** 2) / radsq - (y ** 2) / radysq)
    a[int(nx / 2), int(ny / 2)] = 1.0

    a = array_shift(a, int(cenx - nx / 2), int(ceny - ny / 2))

    # normalise if required
    if norm == True: a *= 1. / np.sum(a)

    return a


def convolve_gaussian(image, rad=3, rady=1):
    ###### Gaussian Convolution Functions by AM
    # convolve an image with a gaussian (blur)
    c = make_gaussian(image.shape[0], image.shape[1], rad, rady, cenx=0, ceny=0, norm=True)
    fc = np.fft.fft2(c)
    fimage = np.fft.fft2(image)
    output = np.real(np.fft.ifft2(np.conjugate(fc) * fimage))
    return output

def extract_r1r2(vol):
    
    nQ = vol.shape[0]
    nTheta = vol.shape[-1]

    r1r2map = np.zeros((nQ, nTheta))

    for i in range(nQ):
        r1r2 = vol[i,i,:]
        r1r2map[i,:]=r1r2

    return r1r2map




def convolve_3D_gaussian(vol, wx, wy, wz, filter_size = 9):
    from scipy import ndimage
    xx,yy,zz = np.meshgrid(np.linspace(-wx,wx,filter_size), np.linspace(-wy,wy,filter_size),np.linspace(-wz,wz,filter_size))
    filter = np.exp(-xx**2 - yy**2 - zz**2)/(np.sqrt(2*np.pi)**3)
    filter -=np.min(filter)
    filter /= np.max(filter)

    new_vol = ndimage.convolve(vol, filter, mode='constant', cval=0.0)

    return new_vol


if __name__ =='__main__':
    








    pdb_code = '4yug'
    nQ =150 
    nTheta = 360
    qmax = 1.25
    
    rmax=(nQ*1e-10)/(2*qmax)*1e9
    fname = f'{pdb_code}-nq{nQ}-nt{nTheta}-qm{qmax}_padf.dbin'
    
    path_to_file = Path(f'/home/pat/rmit-onedrive/phd/python_projects/py3padf02/padf/output/') /fname
   
    arr = read_dbin(str(path_to_file), nQ,nTheta)
    


    r1r2map1 = extract_r1r2(arr)

    
    plot_map(r1r2map1, extent=[0,180, rmax, 0], ylabel='Correlation Distance $r_1$ [nm]',title='4yug')
    
    pdb_code = '4yum'
    nQ =150 
    nTheta = 360
    qmax = 1.25
    
    rmax=(nQ*1e-10)/(2*qmax)*1e9
    fname = f'{pdb_code}-nq{nQ}-nt{nTheta}-qm{qmax}_padf.dbin'
    
    path_to_file = Path(f'/home/pat/rmit-onedrive/phd/python_projects/py3padf02/padf/output/') /fname
   
    arr = read_dbin(str(path_to_file), nQ,nTheta)
    


    r1r2map2 = extract_r1r2(arr)

    
    plot_map(r1r2map2, extent=[0,180, rmax, 0], ylabel='Correlation Distance $r_1$ [nm]',title='4YUM')
    
    
    #plot_map(arr[35,:,:], extent=[0, 180, rmax, 0], ylabel= 'Correlation Distance $r_1$ [nm]', title='$r_1 = index 35$')
    #plot_map(arr[36,:,:],  extent=[0, 180, rmax, 0], ylabel = 'Correlation Distance $r_1$ [nm]', title='$r_1 = index 36$')

    #plot_map(arr[:,:,75], extent=[0, rmax, rmax, 0],ylabel= 'Correlation Distance $r_1$ [nm]', xlabel= 'Correlation Distance $r_1$ [nm]',title='$\Delta = index 75$')
    #plot_map(arr[:,:,76], extent=[0, rmax, rmax, 0],ylabel= 'Correlation Distance $r_1$ [nm]', xlabel= 'Correlation Distance $r_1$ [nm]',title='$\Delta = index 76$')

    
    #plt.figure()
    #plt.plot(np.linspace(0, rmax, nQ),arr[:,45,6])
    #plt.xlabel('Correlation Distance $r_2$ [nm]')
    #plt.ylabel('Correlation Intensity')
    #plt.title('$\Delta =90^o$, $r_1 = 5 [nm]$')

    
    #plt.figure()
    #plt.plot(np.linspace(0, 180, nTheta), arr[2,2,:])
    #plt.xlabel('Correlation Angle $\Delta$ [Degrees]')
    #plt.ylabel('Correlation Intensity')

    #plt.title('$r_1 =2.5 [nm]$, $r_2 = 2.5 [nm]$')



    diff = r1r2map2 - r1r2map1
    plot_map(diff, extent=[0,180, rmax, 0], ylabel='Correlation Distance $r_1$ [nm]',title='diff')
    


    degree= 80
    plt.figure()
    plt.title(f'R1=R2, linplot through $\Delta={degree}$')
    plt.plot(np.linspace(0, rmax, nQ),r1r2map1[:,degree], label='4yug')
    plt.plot(np.linspace(0, rmax, nQ),r1r2map2[:,degree], label=pdb_code)
    plt.plot(np.linspace(0, rmax, nQ),diff[:,degree], label='diff')
    plt.xlabel('Correlation distance R1=R2 [nm]')
    plt.legend()
    


    rind =9 
    plt.figure()
    plt.plot(np.linspace(0, 180, nTheta),r1r2map1[rind,:], label='4yug')
    plt.plot(np.linspace(0, 180, nTheta),r1r2map2[rind,:], label=pdb_code)
    plt.plot(np.linspace(0, 180, nTheta),diff[rind,:], label='diff')
    plt.xlabel('Corelation angle [degrees]')
    plt.legend()



    
    plt.show()



