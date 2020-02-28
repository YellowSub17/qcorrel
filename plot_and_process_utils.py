import matplotlib.pyplot as plt
import numpy as np
import imageio
from pathlib import Path

import struct



def plot_map(map, title='', save=None, cmap='plasma', extent=None, xlabel = 'Correlation Angle $\Delta$ [Degrees]',
             ylabel= 'Scattering Magnitude $q$ [1/$\AA$]'):

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
    plt.show()

def save_tiff(im, name="IMAGE"):

    if name[-4:]=='tiff':
        save_path = Path('tiffs') / f'{name}'
    else:
        save_path = Path('tiffs') / f'{name}.tiff'
    imageio.imsave(save_path, im)


def save_dbin(data, name='DBIN'):

    if name[-4:]=='dbin':
        save_path = Path('tiffs') / f'{name}'
    else:
        save_path = Path('tiffs') / f'{name}.dbin'
    
    f = open(save_path, "wb")
    fmt = '<' + 'd' * data.size
    bin_in = struct.pack(fmt, *data.flatten()[:])
    f.write(bin_in)
    f.close()



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




def convolve_3D_gaussian(vol, wx, wy, wz, filter_size = 9):
    from scipy import ndimage
    xx,yy,zz = np.meshgrid(np.linspace(-wx,wx,filter_size), np.linspace(-wy,wy,filter_size),np.linspace(-wz,wz,filter_size))
    filter = np.exp(-xx**2 - yy**2 - zz**2)/(np.sqrt(2*np.pi)**3)
    filter -=np.min(filter)
    filter /= np.max(filter)

    new_vol = ndimage.convolve(vol, filter, mode='constant', cval=0.0)

    return new_vol


if __name__ =='__main__':
    filter = convolve_3D_gaussian(1,1,1,1, filter_size=10)

    print(filter.shape)
    plot_map(filter[:,0,:])



