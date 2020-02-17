import numpy as np
import spglib





HM_NUMBER = {
    'P 1': 1,
    'P 21 21 21': 115,
    'P 21 21 2': 112,
    'P 21 3': 492,
    'P 2 3': 489,
    'I 41 3 2': 510
}



def apply_sym(reflections, spg_code):


    # Look up the space group number of the cell
    HM_number = HM_NUMBER[spg_code]

    # get the point space gorup (laue group)
    laue_code = spglib.get_spacegroup_type(HM_number)['pointgroup_international']


    if laue_code == 'D2h' or laue_code =='D2':
        return mmm(reflections)

    elif laue_code =='C1' or laue_code == 'Ci':
        return identity(reflections)

    elif laue_code =='T' or laue_code =='Th':
        return m3(reflections)

    elif laue_code=='O' or laue_code=='Td' or laue_code=='Oh':
        return m3m(reflections)


    else:
        print('WARNING: NO SYMETERY APPLIED')
        print(f'spg: {spg_code}, pg: {laue_code}')
        return identity(reflections)




def identity(reflections):
    print('Applying identity symmetry')
    return reflections

def mmm(reflections):
    print('Applying mmm symmetry')
    for i in range(3):
        reflections = mirror(i, reflections)
    return reflections

def m3(reflections):
    print('Applying m3 symmetry')
    reflections  = mirror(2, reflections)
    reflections  = rotate(3, 2, reflections)

    return reflections

def m3m(reflections):
    print('Applying m3m symmetry')
    reflections  = mirror(2, reflections)
    reflections  = rotate(3, 2, reflections)
    reflections  = mirror(2, reflections)
    return reflections


def mirror(m_ind, reflections):

    # create a new reflection array that is double the size of the old array
    new_refle = np.zeros((2 * reflections.shape[0] + 1, 4))

    # for each reflection point
    for i, reflection in enumerate(reflections):
        # add this reflection point to the new reflection list
        new_refle[i, :] = reflection


        #only negate hk or l if it is non zero
        if new_refle[i, m_ind] != 0:
            new_refle[reflections.shape[0] + i, :] = reflection
            new_refle[reflections.shape[0] + i, m_ind] *= -1

    return new_refle

def rotate(nfold, axis, reflections):

    # create a new reflection array that is nfold times the size of the old array
    new_refle = np.zeros(( nfold* reflections.shape[0] + 1, 4))

    # for each set of rotated reflections
    for i in range(nfold):
        #calcualte the angle of rotation
        # note, for i=0 we have theta=0 and so this is the original unrotated point.
        theta = i*2*np.pi/nfold

        # rotation matrices in each yz, xz, xy plane
        Rotx = np.array([[1,0,0],
                         [0, np.cos(theta), -np.sin(theta)],
                         [0, np.sin(theta), np.cos(theta)]])

        Roty = np.array([[np.cos(theta),0,np.sin(theta)],
                         [0, 1, 0],
                         [-np.sin(theta), 0, np.cos(theta)]])

        Rotz = np.array([[np.cos(theta), -np.sin(theta),0],
                         [np.sin(theta), np.cos(theta), 0],
                         [0, 0,1]])


        #axis of rotation is same index within rot
        rot = [Rotx, Roty, Rotz]



        # for each reflection
        for j, reflection in enumerate(reflections):
            # put the rotated hkl components in the new reflections array, round to int
            new_refle[i*reflections.shape[0] + j, :3] = np.round(np.matmul(rot[axis], reflection[:3]))
            # put the hkl intensity into new reflection list (unmodifed by reflection
            new_refle[i*reflections.shape[0] + j, 3] = reflection[3]


    return new_refle





test_hkl = np.random.random((100,4))
test_hkl[0,:3] = [1,1,3]
test_hkl[:, :3] *= 10
test_hkl[:, :3] = np.round(test_hkl[:, :3])
test_hkl[:,3] *= 100
rotated = rotate(8, 2, test_hkl)










