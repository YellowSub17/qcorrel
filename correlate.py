import symmetry as sym
import numpy as np
from pathlib import Path
import concurrent.futures

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
        inten_meas = np.array(inten_meas)
        inten_meas = np.where(inten_meas == '?', '0', inten_meas)
        inten_meas = np.array([float(inten) for inten in inten_meas])
    except KeyError:
        inten_meas = cif[cif.visible_keys[0]]['_refln.F_meas_au']
        inten_meas = np.array(inten_meas)
        inten_meas = np.where(inten_meas == '?', '0', inten_meas)
        inten_meas = np.array([float(inten)**2 for inten in inten_meas])



    # init the list of reflections
    reflections = np.zeros((len(h_refl), 4))

    # input the reflections into the initalize array ( is there a better way?)
    reflections[:, 0] = h_refl
    reflections[:, 1] = k_refl
    reflections[:, 2] = l_refl

    reflections[:, 3] = inten_meas

    reflections = sym.apply_sym(reflections, cif[cif.visible_keys[0]]['_symmetry.space_group_name_H-M'])

    # removes the hkl=000
    reflections = reflections[np.where(np.all(reflections[..., :3] != 0, axis=1))[0]]

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



    alpha =np.radians(float(cif[cif.visible_keys[0]]['_cell.angle_alpha']))

    beta = np.radians(float(cif[cif.visible_keys[0]]['_cell.angle_beta']))
    gamma =np.radians(float(cif[cif.visible_keys[0]]['_cell.angle_gamma']))



    a_unit =   np.array([1, 0, 0], dtype=np.float64)

    b_unit =   np.array([np.cos(gamma),np.sin(gamma) , 0], dtype=np.float64)

    c_unit =   np.array([   np.cos(beta),
                           (np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma),
                            np.sqrt( 1 - np.cos(beta)**2 -  ( (np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma))**2 ) ]
        , dtype=np.float64)

    # get cell paraments to make a,b,c
    a = float(cif[cif.visible_keys[0]]['_cell.length_a'])*a_unit
    b = float(cif[cif.visible_keys[0]]['_cell.length_b'])*b_unit
    c = float(cif[cif.visible_keys[0]]['_cell.length_c'])*c_unit


    # calculate a*, b*, c*
    ast = np.cross(b, c) / np.dot(a, np.cross(b, c))
    bst = np.cross(a, c) / np.dot(a, np.cross(b, c))
    cst = np.cross(a, b) / np.dot(a, np.cross(b, c))

    # for every reflection
    for i, reflection in enumerate(reflections):
        # get the magnitude of the q vector
        q = calc_q(reflection[0], reflection[1], reflection[2], ast, bst, cst)

        qs[i, :3] = q
        qs[i, 3] = np.sqrt(q.dot(q))
        qs[i, 6] = reflection[3]

    return qs


def correlate_qq(cif, dQ, dTheta, theta_bounds=None, q_bounds=None):
    '''
    given a cif file, get an array of q vectors from the reflections, ordered from lowest |q| to highest
    find the correlation shells of thickness dQ_shell for each vector
    correlate each vector with the vectors within it's shell, and map it to a correlation map interpolated to
    have pixels dQ by dTheta
    '''

    # get a list of sorted q vectors in the cif file
    qs_sort = calc_cif_qs(cif)

    if theta_bounds == None:
        theta_min = 0
        theta_max = 185
    else:
        theta_min = theta_bounds[0]
        theta_max = theta_bounds[1]

    if q_bounds == None:
        q_min = 0
        q_max = np.max(qs_sort[:, 3])
    else:
        q_min = q_bounds[0]
        q_max = q_bounds[1]

    # number of angle and q indices
    nTheta = int(round((theta_max - theta_min) / dTheta)) + 1
    nQ = int(round((q_max - q_min) / dQ)) + 1

    # init histogram of theta values
    correl = np.zeros((nQ, nTheta))
    hist = np.zeros((nQ, nTheta))

    for i, q in enumerate(qs_sort):
        if q[3] > q_max or q[3] < q_min:
            qs_sort[i, 4] = -1
            continue

        iq = int(round(q[3] / dQ))
        if (iq >= 0) and (iq < nQ):
            qs_sort[i, 4] = iq
        else:
            print('bing', q, iq)

    # for every q vector
    for i, q in enumerate(qs_sort):

        if q[4] == -1:
            continue

        if i % 400 == 0:
            print(f'Correlating: {i}/{qs_sort.shape[0]} (Q vector {q})', sep='\t')

        iq = int(q[4])

        q_primes = np.where(qs_sort[:, 4] == iq)[0]

        for j in q_primes:

            if i == j:
                continue

            q_prime = qs_sort[j, ...]

            # calculate the angle between the vectors, and the index for theta in the corell map
            theta = np.degrees(angle_between(q[:3], q_prime[:3]))

            if theta > theta_max or theta < theta_min:
                continue

            theta_ind = int(round(((theta - theta_min) / dTheta)))

            # increment the value of the histogram for this q magnitude and angle
            correl[iq, theta_ind] += q[6] * q_prime[6]

            hist[iq, theta_ind] += 1

    if q_bounds == None:
        return correl, hist, q_max

    return correl, hist


def full_correlate(cif, nQ, nTheta, qmax=0.1):

    '''
    Produce a 3D correlation volume of the scattering vectors within a CIF file. nQ and qmax specify the resolution in
    qspace, ntheta defines the angular resolution (theta max locked at 180).
    '''
    # get the scattering vectors
    qs = calc_cif_qs(cif)

    # get indices where q is les then qmax
    correl_vec_indices = np.where(qs[:, 3] < qmax)[0]

    # remove scattering vectors with magnitude less then qmax
    qs = qs[correl_vec_indices]

    # init the correlation and historgrams
    correl = np.zeros((nQ, nQ, nTheta), dtype=np.float32)
   # hist = np.zeros((nQ, nQ, nTheta), dtype=np.float32)

    # for every scattering vector
    for i, q in enumerate(qs):

        if i%400==0:
            print(f'Correlating vector {i}/{len(qs)}. q={q[:3]}')
        # find the position of this vector in the correlation matrix
        q_ind = int(round((q[3] / float(qmax)) * (nQ - 1)))

        #for every other scattering vector
        for q_prime in qs:

            # find the position of this vector in the correlation matrix
            q_prime_ind = int(round((q_prime[3] / float(qmax)) * (nQ - 1)))

            # find the angle between the two vectors
            theta = np.degrees(angle_between(q[:3], q_prime[:3]))

            # find the position of this angle in the correlation matrix
            theta_ind = int(round((theta / 180.0) * (nTheta - 1)))

            # correlate the vectors intensity, add ths value into the matrix
            correl[q_ind, q_prime_ind, theta_ind] += q[6] * q[6]

            # add one to the histogram
    #        hist[q_ind, q_prime_ind, theta_ind] += 1

    return correl#, hist


def correlate_chunk_worker(inputs):
    '''
    worker function for each chunk in the threaded function
    inputs is a tuple of the whole list of q vectors, a list of q vectors just for this chunk,
    nQ, nTheta and qmax to place the calculated values in the matrix
    '''
    # get inputs from the tuple (processes doesn't like multiple inputs?)
    qs, q_primes, nQ, nTheta, qmax = inputs

    # init the volume for just this chunk
    chunk_correl = np.zeros((nQ, nQ, nTheta), dtype=np.float32)
    #chunk_hist = np.zeros((nQ, nQ, nTheta), dtype=np.float32)

    # correlat the vectors as in the other functions.
    for q_i, q in enumerate(qs):
        if q_i%400==0:
            print(f'Correlating vector {q_i}/{len(qs)}. q={q[:3]}')
        q_ind = int(round(((nQ - 1) * (q[3] / qmax))))

        for q_prime_i, q_prime in enumerate(q_primes):
            q_prime_ind = int(round(((nQ - 1) * (q_prime[3] / qmax))))
            theta = np.degrees(angle_between(q[:3], q_prime[:3]))
            theta_ind = int(round((theta / 180.0) * (nTheta - 1)))

            chunk_correl[q_ind, q_prime_ind, theta_ind] += q[6] * q_prime[6]
            #chunk_hist[q_ind, q_prime_ind, theta_ind] += 1

    return chunk_correl#, chunk_hist


#def calc_model_cif_qs(path_to_txt, ):
#    
#    f = open(str(path_to_txt), 'r')
#    cont = f.read()
#    lines = cont.split('\n')
#    for line in lines[:10]:
#    f = np.loadtxt(str(path_to_txt))
#    print(f)









def full_correlate_threaded(cif, nQ, nTheta, qmax, nChunks=4, model_cif=False):
    '''
    Produce a 3D correlation volume of the scattering vectors within a CIF file. nQ and qmax specify the resolution in
    qspace, ntheta defines the angular resolution (theta max locked at 180). Threaded for speed.
    nChunks defines the how the correlation volume is split up (defualt 4, for a quad core cpu)
    '''


    if model_cif:
        qs = calc_model_cif_qs(cif)
    else:
        # get the scattering vectors
        qs = calc_cif_qs(cif)

    # get indices where q is les then qmax
    correl_vec_indices = np.where(qs[:, 3] < qmax)[0]

    # remove scattering vectors with magnitude less then qmax
    qs = qs[correl_vec_indices]
    print(f'Max q: {np.max(qs[:,3])}')

    # init the correlation and historgrams
    correl = np.zeros((nQ, nQ, nTheta), dtype=np.float32)
    #hist = np.zeros((nQ, nQ, nTheta), dtype=np.float32)

    # if the number of correlating vectors is not divisbile by the number of
    # chunks, work out how many are left over. Also work out how big each normal chunk is

    remainder = len(qs) % nChunks
    chunk_size = (len(qs) - remainder) / nChunks

    # init a list where we will place the arguments for each function that will be called
    all_args = []

    # for each chunk (each call of the worker function)
    for chunk_i in range(nChunks):

        # init a list of arguments for this chunk (function)
        args = []

        # find the bounds of the chunk
        b1 = int(chunk_size * chunk_i)
        b2 = int(chunk_size * (chunk_i + 1))

        # if its the last chunck, we need a bit more for the remainder
        if chunk_i == nChunks - 1:
            b2 += remainder

        # append the rest of the arguments for the worker function
        args.append(qs)
        args.append(qs[b1:b2])
        args.append(nQ)
        args.append(nTheta)
        args.append(qmax)

        # append the arguments for this function to the list of all the function arguments
        all_args.append(tuple(args))

    # run a pool of processes for each chunk
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # get the outputs from each function
        results = [executor.submit(correlate_chunk_worker, args) for args in all_args]
        # for each result, sum them together
        for f in concurrent.futures.as_completed(results):
            correl += f.result()
            #hist += f.result()[1]

    return correl#, hist



def print_qmax(cif, qmaxs= []):
    '''
    Produce a 3D correlation volume of the scattering vectors within a CIF file. nQ and qmax specify the resolution in
    qspace, ntheta defines the angular resolution (theta max locked at 180). Threaded for speed.
    nChunks defines the how the correlation volume is split up (defualt 4, for a quad core cpu)
    '''


    # get the scattering vectors
    qs = calc_cif_qs(cif)


    print(f'Max q: {np.max(qs[:,3])}')
    for qmax in qmaxs:
        # get indices where q is les then qmax
        correl_vec_indices = np.where(qs[:, 3] < qmax)[0]

        # remove scattering vectors with magnitude less then qmax
        qs_of_qmax = qs[correl_vec_indices]


        # get indices where q is les then qmax
        print(f'qmax={qmax}, Number of Q vectors {len(qs_of_qmax)}')



if __name__ == '__main__':
    import CifFile
    import plot_and_process_utils as ppu
    import time
    import sys
    import os





    alpha_path = Path('cifs/alpha')
    GFP_path = Path('cifs/GFP')
    keratin_path = Path('cifs/alpha/keratin')
    macro_path = Path('cifs/macro')
    model_path = Path('cifs/sf_from_model')
    cif_file_names = []


#    # Normal size
#    cif_file_names.append(alpha_path / '1al1-sf.cif')
#    cif_file_names.append(alpha_path / '1mft-sf.cif')
#    cif_file_names.append(alpha_path / '1cos-sf.cif')
#    qmaxs=[0.3,0.39,0.35]



#   # High Res
#    cif_file_names.append(alpha_path / '1al1-sf.cif')
#    cif_file_names.append(alpha_path / '1cos-sf.cif')
#    qmaxs=[0.3699, 0.5]


#   # mid Res for comp
#    cif_file_names.append(alpha_path / '1cos-sf.cif')
#    qmaxs=[0.3699]



    # cif_file_names.append(alpha_path/'4omz-sf.cif')
    # cif_file_names.append(alpha_path/'4osd-sf.cif')
    # qmaxs = [0.3699/2,0.3699/2 ]


    cif_file_names.append(alpha_path/'4omz-sf.cif')
    cif_file_names.append(alpha_path/'1cos-sf.cif')
    cif_file_names.append(alpha_path/'1al1-sf.cif')
    qmaxs = [0.3699/2,0.3699,0.3699 ]




    nq=256
    nt=360



    for cif, qmax in zip(cif_file_names, qmaxs):

        print(f'Reading {cif}')







        sf_cif = CifFile.ReadCif(str(cif))

        print_qmax(sf_cif, [qmax, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6])
# #        print('\n\n')


        # start =time.strftime('%H:%M:%S, %A %d %B')
        # print(f'starting at {start}')
        # correl = full_correlate_threaded(sf_cif, nq,nt, qmax,nChunks=4)
        # end =time.strftime('%H:%M:%S, %A %d %B')
        # print(f'ended at {end}')


        # dbin_fname = str(cif.stem+'_paper_qcorrel')


        # ppu.save_dbin(correl,dbin_fname)

        # log_fname=str(Path('dbins')/  f'{dbin_fname}_log.txt')
        # ppu.write_log(log_fname, cif=cif, nQ=nq, nTheta=nt, qmax=qmax, Start_Time=start,End_Time=end)

