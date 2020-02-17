import time
import CifFile
import numpy as np
import matplotlib.pyplot as plt
import symmetry as sym

def calc_q(h, k, l, ast, bst, cst):
    '''
    Calculate the q vector from a hkl index and reciprocal lattice vector
    '''
    return h * ast + k * bst + l * cst


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))



def get_cif_reflections(cif):
    '''
    Given a cif file, make an array with coloumns h,k,l,intensity
    then apply the relevant symmetry operations on the hkl
    then return the list of all refelctions
    '''
    # read in reflections and convert str to int
    h_refl = cif[cif.visible_keys[0]]['_refln.index_h']
    h_refl = np.array([int(h) for h in h_refl])

    k_refl = cif[cif.visible_keys[0]]['_refln.index_k']
    k_refl = np.array([int(k) for k in k_refl])

    l_refl = cif[cif.visible_keys[0]]['_refln.index_l']
    l_refl = np.array([int(l) for l in l_refl])

    # read in inten, but also remove ? marks (assume low intensity, 0)
    try:
        inten_meas = cif[cif.visible_keys[0]]['_refln.intensity_meas']
    except KeyError:
        inten_meas = cif[cif.visible_keys[0]]['_refln.F_meas_au']

    inten_meas = np.array(inten_meas)
    inten_meas = np.where(inten_meas == '?', '0', inten_meas)
    inten_meas = np.array([float(inten) for inten in inten_meas])

    # init the list of reflections
    reflections = np.zeros((len(h_refl), 4))

    # input the reflections into the initalize array ( is there a better way?)
    reflections[:, 0] = h_refl
    reflections[:, 1] = k_refl
    reflections[:, 2] = l_refl

    reflections[:, 3] = inten_meas


    reflections = sym.apply_sym(reflections, cif[cif.visible_keys[0]]['_symmetry.space_group_name_H-M'])


    reflections=reflections[~np.all(reflections == 0, axis=1)]

    return reflections


def calc_cif_qs(cif):

    '''
    given a cif file, find the reflections and calculated the scattering vector q for each reflection
    return a an array with coloumns qx,qy,qz, |q|, 0,0, intens. (ordered my |q| lowest to highest down the col)
    note that coloumsn of 0 are for upper and lower bounds later, and that the cif is searched for a space group to
    apply reflections
    '''
    reflections = get_cif_reflections(cif)

    qs = np.zeros((reflections.shape[0], 7), dtype=np.float64)

    # get cell paraments to make a,b,c
    a_len = float(sf_cif[cif.visible_keys[0]]['_cell.length_a']) * np.array([1, 0, 0], dtype=np.float64)
    b_len = float(sf_cif[cif.visible_keys[0]]['_cell.length_b']) * np.array([0, 1, 0], dtype=np.float64)
    c_len = float(sf_cif[cif.visible_keys[0]]['_cell.length_c']) * np.array([0, 0, 1], dtype=np.float64)

    # calculate a*, b*, c*
    ast_len = np.cross(b_len, c_len) / np.dot(a_len, np.cross(b_len, c_len))
    bst_len = np.cross(a_len, c_len) / np.dot(a_len, np.cross(b_len, c_len))
    cst_len = np.cross(a_len, b_len) / np.dot(a_len, np.cross(b_len, c_len))

    # for every reflection
    for i, reflection in enumerate(reflections):
        # get the magnitude of the q vector
        q = calc_q(reflection[0], reflection[1], reflection[2], ast_len, bst_len, cst_len)

        qs[i, :3] = q
        qs[i, 3] = np.sqrt(q.dot(q))
        qs[i, 6] = reflection[3]


    return qs


def correlate1(cif, dQ_shell, dQ, dTheta, hist=True, out_q_max=False):

    '''
    given a cif file, get an array of q vectors from the reflections, ordered from lowest |q| to highest
    find the correlation shells of thickness dQ_shell for each vector
    correlate each vector with the vectors within it's shell, and map it to a correlation map interpolated to
    have pixels dQ by dTheta
    '''
    # get a list of sorted q vectors in the cif file
    qs_sort = calc_cif_qs(cif)
    print(f'Max q magnitude: {qs_sort[-1, 3]}')

    # number of angle and q indices
    nTheta = int(360 / dTheta)
    nQ = int(qs_sort[-1, 3] / dQ)

    # init histogram of theta values
    correl = np.zeros((nQ, nTheta))

    # for every q vector, we want to find how many indices above and below we need to go before the difference in |q| > dQ
    print('Finding upper and lower bounds of correlation shells')
    for i, q in enumerate(qs_sort):
        # init the upper and lower index values, how far above and below we look to correlate with
        upper_ind = 0
        lower_ind = 0

        # current difference in magnitude above and below is 0
        upper_mag = 0
        lower_mag = 0

        # continously check if the current difference in magnitude in vectors below is less then dQ
        while abs(lower_mag) < abs(dQ_shell):
            # if it is, then make sure that highest the lower_ind can go is the current position in the list
            if i - lower_ind == 0:
                break  # this avoids index errors
            # get the current comparison q vector
            q_prime = qs_sort[i - lower_ind]

            # get the difference in magnitude between the two q vectors, this gets checked in the next loop
            lower_mag = q[3] - q_prime[3]
            # increment how indices that we look below in the next iteration
            lower_ind += 1

        # repeat the process, but now calculate the upper bound
        while abs(upper_mag) < abs(dQ_shell):
            if i + upper_ind == qs_sort.shape[0]:
                break
            q_prime = qs_sort[i + upper_ind]
            upper_mag = q_prime[3] - q[3]

            upper_ind += 1

        # once we have the upper and lower bounds for our q shell, save these values in qs_sort

        qs_sort[i, 4] = lower_ind
        qs_sort[i, 5] = upper_ind

    # for every q vector
    for i, q in enumerate(qs_sort):
        print(f'Correlating: {i}/{qs_sort.shape[0]} (Q vector {q})', sep='\t')

        # extract the upper and lower bounds of correlation
        q_lower_ind = int(q[4])
        q_upper_ind = int(q[5])

        # for every vector within the bounds of +/- dQ (q[4]:q[5])
        for j, q_prime in enumerate(qs_sort[q_lower_ind:q_upper_ind]):

            # calculate the angle between the vectors, and the index for theta in the corell map
            theta = np.degrees(angle_between(q[:3], q_prime[:3]))
            theta_ind = int((theta / 360.0) * (nTheta - 1))

            # reflection about 180
            #theta_ind2 = int(((360 - theta) / 360.0) * (nTheta - 1))

            # calculate which index for q based in magnitude
            q_ind = int((q[3] / qs_sort[-1, 3]) * (nQ - 1))

            # increment the value of the histogram for this q magnitude and angle
            if hist:
                correl[q_ind, theta_ind] += 1
                #correl[q_ind, theta_ind2] += 1
            else:
                correl[q_ind, theta_ind] += q[6] * q_prime[6]
                #correl[q_ind, theta_ind2] += q[6] * q_prime[6]


    if out_q_max:
        return correl, qs_sort[-1, 3]
    else:
        return correl



def correlate(cif, dQ, dTheta):

    '''
    given a cif file, get an array of q vectors from the reflections, ordered from lowest |q| to highest
    find the correlation shells of thickness dQ_shell for each vector
    correlate each vector with the vectors within it's shell, and map it to a correlation map interpolated to
    have pixels dQ by dTheta
    '''

    # get a list of sorted q vectors in the cif file
    qs_sort = calc_cif_qs(cif)


    # number of angle and q indices
    nTheta = int(360 / dTheta)
    nQ = int(np.max(qs_sort[:, 3]) / dQ)
    #print(nQ)
    # init histogram of theta values
    correl = np.zeros((nQ, nTheta))
    hist = np.zeros((nQ, nTheta))



    for i, q in enumerate(qs_sort):
        iq = int(q[3]/dQ)
        if (iq>=0) and (iq < nQ):
            qs_sort[i,4] = iq




    # for every q vector
    for i, q in enumerate(qs_sort):

        if i%400== 0:
            print(f'Correlating: {i}/{qs_sort.shape[0]} (Q vector {q})', sep='\t')

        # extract the upper and lower bounds of correlation
        iq = int(q[4])

        # for every vector within the bounds of +/- dQ (q[4]:q[5])
        for j, q_prime in enumerate(qs_sort):
            iq_prime = int(q_prime[4])

            if i==j or iq != iq_prime:
                continue


            # calculate the angle between the vectors, and the index for theta in the corell map
            theta = np.degrees(angle_between(q[:3], q_prime[:3]))
            theta_ind = int((theta / 180.0) * (nTheta - 1))


            # increment the value of the histogram for this q magnitude and angle
            correl[iq, theta_ind] += q[6]*q_prime[6]

            hist[iq, theta_ind] += 1


    return correl, hist







protein_names = ['5dk8']

dqs = [0.001]

hkl_num = [10]

# Read in Cif file
for name in protein_names:
    for num in hkl_num:

        print('Reading Cif')
        sf_cif = CifFile.ReadCif(f'cifs\\{name}-sf_less{num}.cif')
        for dq in dqs:

            print(f'STARTING NEW - {name} {num} {dq}')

            correl, hist = correlate(sf_cif, dq,  0.5)


            plt.figure()
            plt.title(f'{name}_less{num} correl - dq={dq}')
            plt.imshow(correl, cmap='plasma', aspect='auto')
            plt.colorbar()
            # plt.savefig(f'saved_plots\\{name}_less{num}_correl_dq{int(dq*1000)}')
            # plt.close('all')

            plt.figure()
            plt.title(f'{name}_less{num} log10 correl - dq={dq}')
            plt.imshow(np.log10(correl), cmap='plasma', aspect='auto')
            plt.colorbar()
            # plt.savefig(f'saved_plots\\{name}_less{num}_logcorrel_dq{int(dq*1000)}')
            # plt.close('all')
            
            
            plt.figure()
            plt.title(f'{name}_less{num} hist - dq={dq}')
            plt.imshow(hist, cmap='plasma', aspect='auto')
            plt.colorbar()

            plt.figure()
            plt.title(f'{name}_less{num} log10 hist - dq={dq}')
            plt.imshow(np.log10(hist), cmap='plasma', aspect='auto')
            plt.colorbar()
            # plt.savefig(f'saved_plots\\{name}_less{num}_hist_dq{int(dq*1000)}')
            # plt.close('all')

            print('\n\n\n')




