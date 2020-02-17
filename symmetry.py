import numpy as np
import spglib





HM_NUMBER_DICT = {
    'P 1': 1,
    'P 21 21 21': 115,
    'P 21 21 2': 112,
    'P 21 3': 492,
    'P 2 3': 489,
    'I 41 3 2': 510
}


def four_over_m(h,k,l):


    multiplicity = 8

    total_sym_mat = np.zeros( (multiplicity,3, 3) )


    total_sym_mat[0,...] = np.array([[1,0,0], [0,1,0], [0,0,1]])
    total_sym_mat[1,...] = np.array([[-1,0,0], [0,-1,0], [0,0,1]])
    total_sym_mat[2,...] = np.array([[0,-1,0], [1,0,0], [0,0,1]])
    total_sym_mat[3,...] = np.array([[-1,0,0], [0,-1,0], [0,0,-1]])




    ops_ind = 4

    while np.array_equal(total_sym_mat[-1,...], np.zeros((3,3))):

        for i in range(ops_ind):
            for j in range(ops_ind):
                new_sym_mat = np.matmul(total_sym_mat[i,...], total_sym_mat[j,...])


                is_new_sym=True

                for sym_mat in total_sym_mat[:ops_ind,...]:

                    if np.array_equal(sym_mat, new_sym_mat):
                        is_new_sym=False
                        break

                if is_new_sym:

                    total_sym_mat[ops_ind,...] = new_sym_mat
                    ops_ind +=1


    return total_sym_mat

#
# def apply_sym(reflections, spg_code):
#
#
#     # Look up the space group number of the cell
#     HM_number = HM_NUMBER_DICT[spg_code]
#
#     # get the point space gorup (laue group)
#
#
#     if HM_number ==1 or HM_number ==2:
#         return mmm(reflections)
#
#     elif laue_code =='C1' or laue_code == 'Ci':
#         return identity(reflections)
#
#     elif laue_code =='T' or laue_code =='Th':
#         return m3(reflections)
#
#     elif laue_code=='O' or laue_code=='Td' or laue_code=='Oh':
#         return m3m(reflections)
#
#
#     else:
#         print('WARNING: NO SYMETERY APPLIED')
#         print(f'spg: {spg_code}, pg: {laue_code}')
#         return identity(reflections)
#
#
#


ans = four_over_m(1,1,1)
print(*ans, sep='\n\n')

#
#
# test_hkl = np.random.random((100,4))
# test_hkl[0,:3] = [1,1,3]
# test_hkl[:, :3] *= 10
# test_hkl[:, :3] = np.round(test_hkl[:, :3])
# test_hkl[:,3] *= 100
#
#
#








