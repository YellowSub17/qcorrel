import matplotlib.pyplot as plt
import numpy as np
import imageio
from pathlib import Path

import struct



def plot_map(map, title='',min_max = False, cb=True, save=None, cmap='plasma', extent=None,dpi=100,
        xlabel = '', ylabel= '', aspect='auto', vmin=None, vmax=None, origin='lower', fig_size=(6.4,4.8)):

    if min_max:
        map = map-np.min(map)
        map = map/np.max(map)
    plt.figure(figsize=fig_size, dpi=dpi)
    plt.title(title)
    plt.imshow(map, cmap=cmap, aspect=aspect, extent=extent, vmin=vmin, vmax=vmax, origin=origin)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if cb:
        plt.colorbar()
    if type(save)==str:
        save_path = f'{save}.png'
        plt.savefig(str(save_path))






def save_tiff(im, name="IMAGE"):

    if name[-4:]=='tiff':
        save_path =  f'{name}'
    else:
        save_path =  f'{name}.tiff'
    imageio.imsave(save_path, im.astype(np.float32))


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


def read_blrr(path, nQ):

    if path[-4:]=='dbin':
        save_path = Path( f'{path}')
    else:
        save_path = Path(f'{path}.dbin')
    arr = np.fromfile(str(save_path))
    arr = arr.reshape(nQ,nQ)

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


def convolve_gaussian(image, rad=1, rady=1):
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




def convolve_3D_gaussian(vol, wx, wy, wz, filter_size = 9, mode='constant', cval=0.0):
    from scipy import ndimage
    xx,yy,zz = np.meshgrid(np.linspace(-wx/2,wx/2,filter_size), np.linspace(-wy/2,wy/2,filter_size),np.linspace(-wz/2,wz/2,filter_size))
    g_filter = np.exp(-xx**2 - yy**2 - zz**2)#/(np.sqrt(2*np.pi)**3)
    #g_filter -=np.min(g_filter)
    #g_filter /= np.max(g_filter)

    new_vol = ndimage.convolve(vol, g_filter, mode=mode, cval=cval)

    return new_vol#, g_filter


def write_log(path, **kwargs):

    f=open(path, 'w')
    f.write('## Q space correlation log\n\n')
    for kw in kwargs:
        f.write(f'{kw} = {kwargs[kw]}\n\n' )
    f.close()

#if __name__ =='__main__':
#    import os
#    import matplotlib
#
#    def animate_vol(vol, interval=100, view='Q', start=0,end=-1 ):
#
#        if view=='Q':
#            pass
#        elif view=='T':
#            vol = np.rot90(vol, axes=(0,-1))
#
#        fig, ax = plt.subplots()
#
#        vol = vol[start:end, :,:]
#        ims = []
#        for i in range(vol.shape[0]):
#
#            im = ax.imshow(vol[i,:,:], animated=True)
#            title = ax.text(0.5,1.05,f"Frame {i}",
#                            size=plt.rcParams["axes.titlesize"],
#                            ha="center", transform=ax.transAxes, )
#
#            ims.append([im, title])
#
#        ani = animation.ArtistAnimation(fig, ims, interval=interval, blit=False)
#        return ani
#
#
#
#    padf_path = Path(os.getcwd()).parent /'py3padf02'/'padf'/'output'
#
#
#
#
#    vol = read_dbin(str(padf_path/'4ggr-sf_qcorrel_padf.dbin'),256, 360, half_theta=True)
#    ani = animate_vol(vol, view='Q')
#    plt.show()
#
 
