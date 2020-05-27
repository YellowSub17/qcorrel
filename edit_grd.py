



import numpy as np
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

import matplotlib.animation as animation


def read_grd(fname):
    grd_path = Path(f'grds/{fname}.grd')

    grd_file = open(str(grd_path), 'r')




    huge_flat_list = []
    grd_file_slices = grd_file.read().split('\n\n')

    arr_shape = grd_file_slices[0].split('\n')[2].split()


    arr_shape = (int(arr_shape[0]), int(arr_shape[1]),int(arr_shape[2]))


    for slice in grd_file_slices:
        slice_lines = slice.split('\n')
        #print(slice_lines)
        for slice_line in slice_lines:
            vals = slice_line.split()
            for val in vals:
                huge_flat_list.append(val)
    huge_flat_arr = np.array([ float(i) for i in huge_flat_list[10:]])

    vol = huge_flat_arr.reshape(arr_shape)
    return vol


def normalize1(arr):
    arr-=np.min(arr)
    arr/=np.max(arr)

    return arr






if __name__=='__main__':



    vol1 = read_grd('1nqp-ed_res1,7_iso0,5')



    vol1= normalize1(vol1**2)

    print(np.where(vol1==np.max(vol1)))


    print(vol1.shape)
    plt.imshow(np.sum(vol1[-4:-2, :,:], axis=0))
    #
    plt.show()
    # #
    # #
    #

