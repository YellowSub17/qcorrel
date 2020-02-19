
import numpy as np
import symmetry as sym


def calc_q(h, k, l, ast, bst, cst):
    '''
    Calculate the q vector from a hkl index and reciprocal lattice vector
    '''
    return h * ast + k * bst + l * cst


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
   # print(vector, np.linalg.norm(vector))
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



    #removes the hkl=000
    reflections = reflections[np.where(np.all(reflections[...,:3]!=0, axis=1))[0]]

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
    a_len = float(cif[cif.visible_keys[0]]['_cell.length_a']) * np.array([1, 0, 0], dtype=np.float64)
    b_len = float(cif[cif.visible_keys[0]]['_cell.length_b']) * np.array([0, 1, 0], dtype=np.float64)
    c_len = float(cif[cif.visible_keys[0]]['_cell.length_c']) * np.array([0, 0, 1], dtype=np.float64)

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


def correlate(cif, dQ, dTheta, theta_bounds = None, q_bounds=None):

    '''
    given a cif file, get an array of q vectors from the reflections, ordered from lowest |q| to highest
    find the correlation shells of thickness dQ_shell for each vector
    correlate each vector with the vectors within it's shell, and map it to a correlation map interpolated to
    have pixels dQ by dTheta
    '''



    # get a list of sorted q vectors in the cif file
    qs_sort = calc_cif_qs(cif)




    if theta_bounds ==None:
        theta_min = 0
        theta_max = 185
    else:
        theta_min = theta_bounds[0]
        theta_max = theta_bounds[1]


    if q_bounds==None:
        q_min = 0
        q_max = np.max(qs_sort[:, 3])
    else:
        q_min = q_bounds[0]
        q_max = q_bounds[1]



    # number of angle and q indices
    nTheta = int(round( (theta_max-theta_min) / dTheta))+1
    nQ = int(round( (q_max-q_min)/ dQ))+1
    #print(nQ)
    # init histogram of theta values

    correl = np.zeros((nQ, nTheta))
    hist = np.zeros((nQ, nTheta))



    for i, q in enumerate(qs_sort):
        if q[3] > q_max or q[3] < q_min:
            qs_sort[i,4] = -1
            continue

        iq = int(round(q[3]/dQ))
        if (iq>=0) and (iq < nQ):
            qs_sort[i,4] = iq
        else:
            print('bing', q, iq)




    # for every q vector
    for i, q in enumerate(qs_sort):

        if q[4] == -1:
            continue

        if i%400== 0:
            print(f'Correlating: {i}/{qs_sort.shape[0]} (Q vector {q})', sep='\t')


        iq = int(q[4])

        q_primes = np.where(qs_sort[:,4]==iq)[0]


        for j in q_primes:

            if i==j:
                continue


            q_prime = qs_sort[j,...]


            # calculate the angle between the vectors, and the index for theta in the corell map
            theta = np.degrees(angle_between(q[:3], q_prime[:3]))


            if theta > theta_max or theta < theta_min:

                continue


            theta_ind = int(round(((theta-theta_min) / dTheta)))


            # increment the value of the histogram for this q magnitude and angle
            correl[iq, theta_ind] += q[6]*q_prime[6]

            hist[iq, theta_ind] += 1


    if q_bounds==None:
        return correl, hist, q_max


    return correl, hist






