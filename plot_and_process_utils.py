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



def plot_r1r2map(map, title='', save=None, cmap='plasma', extent=None,
        xlabel = 'Correlation Angle $\Delta$ [Degrees]', ylabel= 'Correlation distance $r_1=r_2$ [$\AA$]'):

    plt.figure()
    plt.title(title)
    plt.imshow(map, cmap=cmap, aspect='auto', extent=extent)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar()
    if type(save)==str:
        save_path = Path('saved_plots') / f'{save}.png'
        plt.savefig(str(save_path))


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
    


    padf_out_path = Path('/home/pat/rmit-onedrive/phd/python_projects/py3padf02/padf/output')
    dbins = []
    
    dbins.append(padf_out_path / '1cos-sf_res4_qcorrel_padf.dbin')


    
    
    nQ = 150
    qmax = 0.08
    nTheta = 180 
    rmax =(nQ*1e-10)/(2*qmax)*1e9


    r_scale = np.linspace(0, rmax, nQ)**2
    theta_scale = np.linspace(0, 180, nTheta)

    tt,rr = np.meshgrid(theta_scale, r_scale)

    

    rcrop_ind =0
    for dbin in dbins:
        arr = read_dbin(str(dbin), nQ, nTheta)

        r1r2map = extract_r1r2(arr)
        
        plot_map(r1r2map[rcrop_ind:,:]*rr[rcrop_ind:,:], extent=[0,nTheta,rmax, r_scale[rcrop_ind]])
        plt.figure()
        plt.plot(np.linspace(0,180,nTheta), r1r2map[10, :]/np.max(r1r2map[10,:]))
    plt.show()
    


