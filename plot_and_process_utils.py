import matplotlib.pyplot as plt
import numpy as np
import imageio
from pathlib import Path

import struct



def plot_map(map, title='', save=None, cmap='plasma', extent=None,
        xlabel = '', ylabel= ''):

    plt.figure()
    plt.title(title)
    plt.imshow(map, cmap=cmap, aspect='auto', extent=extent)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
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




def convolve_3D_gaussian(vol, wx, wy, wz, filter_size = 9, mode='constant', cval=0.0):
    from scipy import ndimage
    xx,yy,zz = np.meshgrid(np.linspace(-wx/2,wx/2,filter_size), np.linspace(-wy/2,wy/2,filter_size),np.linspace(-wz/2,wz/2,filter_size))
    g_filter = np.exp(-xx**2 - yy**2 - zz**2)#/(np.sqrt(2*np.pi)**3)
    #g_filter -=np.min(g_filter)
    #g_filter /= np.max(g_filter)

    new_vol = ndimage.convolve(vol, g_filter, mode=mode, cval=cval)

    return new_vol, g_filter


def write_log(path, **kwargs):

    f=open(path, 'w')
    f.write('## Q space correlation log\n\n')
    for kw in kwargs:
        f.write(f'{kw} = {kwargs[kw]}\n\n' )
    f.close()

if __name__ =='__main__':


    dbin_path = Path('dbins/convol_test/')

    nQ = 256
    nTheta = 360
    qmax = 0.14


    dbins = []

    dbins.append(dbin_path / '253l-sf_ave_qcorrel.dbin')

    dbins.append(dbin_path / '254l-sf_ave_qcorrel.dbin')


    convols = [4,8,12,16,20,24]
    for dbin in dbins:
        print(str(dbin))
        qcorrel_vol = read_dbin(str(dbin), nQ, nTheta)

        for filter_size in convols:
            print(filter_size)
            convolved_qcorrel_vol, x = convolve_3D_gaussian(qcorrel_vol, 6,6,6, filter_size, mode='reflect')
            save_dbin(convolved_qcorrel_vol,f'convol_test/{dbin.stem}_filter{filter_size}')





#   # Convol test
#   vol = np.random.random((100,100,100))
#   plt.figure()
#   plt.imshow(vol[:,:,50])


#   for filter_size in [1]:# range(4,25 ,4):
#       vol1, g_filter1 = convolve_3D_gaussian(vol,6,6,6,16, mode='nearest')

#       plt.figure()
#       plt.imshow(vol1[:,:,50])
#       plt.figure()
#       plt.imshow(g_filter1[:,:,int(filter_size/2)])

#       vol1, g_filter1 = convolve_3D_gaussian(vol,6,6,6,16, mode='mirror')

#       plt.figure()
#       plt.imshow(vol1[:,:,50])
#       plt.figure()
#       plt.imshow(g_filter1[:,:,int(filter_size/2)])


#       vol1, g_filter1 = convolve_3D_gaussian(vol,6,6,6,16, mode='wrap')

#       plt.figure()
#       plt.imshow(vol1[:,:,50])
#       plt.figure()
#       plt.imshow(g_filter1[:,:,int(filter_size/2)])


#       vol1, g_filter1 = convolve_3D_gaussian(vol,6,6,6,16, mode='reflect')

#       plt.figure()
#       plt.imshow(vol1[:,:,50])
#       plt.figure()
#       plt.imshow(g_filter1[:,:,int(filter_size/2)])


#       vol1, g_filter1 = convolve_3D_gaussian(vol,6,6,6,16, mode='constant', cval=1.0)

#       plt.figure()
#       plt.imshow(vol1[:,:,50])
#       plt.figure()
#       plt.imshow(g_filter1[:,:,int(filter_size/2)])



#       vol1, g_filter1 = convolve_3D_gaussian(vol,6,6,6,16, mode='constant', cval=0.0)

#       plt.figure()
#       plt.imshow(vol1[:,:,50])
#       plt.figure()
#       plt.imshow(g_filter1[:,:,int(filter_size/2)])



#   plt.show()

