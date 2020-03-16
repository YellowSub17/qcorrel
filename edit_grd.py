



import numpy as np
from pathlib import Path




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
    import matplotlib.pyplot as plt

    vol254l = read_grd('254l-ed_res0,7_iso1')
    vol253l = read_grd('253l-ed_res0,7_iso1')


    vol254l = normalize1(vol254l**2)
    vol253l = normalize1(vol253l**2)

    vol = vol254l-vol253l

    print(np.where(vol==np.max(vol)))



    plt.imshow(vol[45,:,:])

    plt.show()




