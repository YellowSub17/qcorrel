import time
import CifFile
import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt




def diff_cif(fname1, fname2, bins=50):
    '''
    anlayize the difference in reflections between two cif files.
    '''
    print('Reading Cif 1')
    cif1 = CifFile.ReadCif(str(fname1))

    print('Reading Cif 2')
    cif2 = CifFile.ReadCif(str(fname2))

    h1, k1, l1, I1, Isig1 = import_hklIF(cif1)
    h2, k2, l2, I2, Isig2 = import_hklIF(cif2)

    refl1 = [(int(h), int(k), int(l)) for h, k, l in zip(h1, k1, l1)]
    refl2 = [(int(h), int(k), int(l)) for h, k, l in zip(h2, k2, l2)]



    diffs = []
    for i, refl in enumerate(set(refl1) &set(refl2)):
        diffs.append(I1[refl1.index(refl)] - I2[refl2.index(refl)])

    diffs = np.array(diffs)

    print(f'Maximum difference: {np.max(diffs)}')

    print(f'Minimum difference: {np.min(diffs)}')

    print(f'Average difference: {np.mean(diffs)}')

    plt.figure()
    plt.hist(diffs, bins=bins)

    plt.title(f'Difference in intensity: {fname1}, {fname2}')
    plt.xlabel('Difference')
    plt.ylabel('Counts')



def comp_cif(fname1, fname2):
    '''
    Compare the reflections between two cif files.
    Save a new cif file with the reflections that are found in both, with intensity values of fname1
    '''
    print('Reading Cif 1')
    cif1 = CifFile.ReadCif(str(fname1))

    print('Reading Cif 2')
    cif2 = CifFile.ReadCif(str(fname2))

    h1, k1, l1, I1, Isig1 = import_hklIF(cif1)
    h2, k2, l2, I2, Isig2 = import_hklIF(cif2)

    refl1 = [(int(h), int(k), int(l)) for h, k, l in zip(h1, k1, l1)]
    refl2 = [(int(h), int(k), int(l)) for h, k, l in zip(h2, k2, l2)]

    intersection = list(set(refl1) & set(refl2))



    new_h = []
    new_k = []
    new_l = []
    new_I = []
    new_Isig = []


    for i, refl in enumerate(intersection):
        new_h.append(refl[0])
        new_k.append(refl[1])
        new_l.append(refl[2])
        #DEBUG PRINT STATEMENTS
        print(refl)
        # print(i, I1[i])
        print(refl1.index(refl), I1[refl1.index(refl)])
        print('\n\n')
        new_I.append(I1[refl1.index(refl)])
        new_Isig.append( Isig1[refl1.index(refl)])






    # convert the ints to str
    new_h = [str(h) for h in new_h]
    new_k = [str(k) for k in new_k]
    new_l = [str(l) for l in new_l]

    # save str list back into cif dict
    cif1[ cif1.visible_keys[0]]['_refln.index_h'] = new_h
    cif1[ cif1.visible_keys[0]]['_refln.index_k'] = new_k
    cif1[ cif1.visible_keys[0]]['_refln.index_l'] = new_l

    new_I = [str(I) for I in new_I]
    new_Isig = [str(Isig) for Isig in new_Isig]

    if '_refln.intensity_meas' in cif1[cif1.visible_keys[0]]:
        cif1[ cif1.visible_keys[0]]['_refln.intensity_meas'] = new_I
        cif1[ cif1.visible_keys[0]]['_refln.intensity_sigma'] = new_Isig

    elif '_refln.F_meas_au' in cif1[ cif1.visible_keys[0]]:
        cif1[ cif1.visible_keys[0]]['_refln.F_meas_au'] = new_I
        cif1[ cif1.visible_keys[0]]['_refln.F_meas_sigma_au'] = new_Isig
    else:
        print('CHECK CIF ENTRY NAME!!')
        return None







    # write the dict to file

    out_fname = fname1.stem + f'_comp{fname2.stem}.cif'

    outfile = fname1.parent / out_fname

    print(f'Saving to:{outfile}')
    out_file = open(str(outfile), 'w')
    out_file.write(cif1.WriteOut())
    out_file.close()

    # debug = np.array([new_h, new_k, new_l, new_I, new_Isig])
    # return debug


def cut_hkl(fname, resolutions):
    '''
    reduce the number of reflections in a CIF file fname
    if h,k or l are greater then bound, ignore this reflection
    remaining reflections are saved as outfname
    '''
    print(f'Reading Cif {fname}')
    cif = CifFile.ReadCif(str(fname))  # read the cif file

    cif_key = cif.visible_keys[0]

    a_len = float(cif[cif_key]['_cell.length_a'])
    b_len = float(cif[cif_key]['_cell.length_b'])
    c_len = float(cif[cif_key]['_cell.length_c'])



    h_refl_orig, k_refl_orig, l_refl_orig, I_refl_orig, I_sig_orig = import_hklIF(cif)

    for res in resolutions:

        print(f'Computing resolution: {res}')
        # int new lists
        new_h = []
        new_k = []
        new_l = []
        new_I = []
        new_Isig = []

        q_max = 1.0 / (2 * res)

        h_bound = int(round(a_len * q_max))
        k_bound = int(round(b_len * q_max))
        l_bound = int(round(c_len * q_max))

        print(f'HKL BOUNDS: {h_bound}, {k_bound}, {l_bound}')

        for h, k, l, I, Isig in zip(h_refl_orig, k_refl_orig, l_refl_orig, I_refl_orig, I_sig_orig):
            if abs(h) > h_bound or abs(k) > k_bound or abs(l) > l_bound:
                continue  # skip reflections greater than bounds

            # append indices less then bounds to the new lists
            new_h.append(h)
            new_k.append(k)
            new_l.append(l)
            new_I.append(I)
            new_Isig.append(Isig)

        # convert the ints to str
        new_h = [str(h) for h in new_h]
        new_k = [str(k) for k in new_k]
        new_l = [str(l) for l in new_l]

        # save str list back into cif dict
        cif[cif_key]['_refln.index_h'] = new_h
        cif[cif_key]['_refln.index_k'] = new_k
        cif[cif_key]['_refln.index_l'] = new_l

        new_I = [str(I) for I in new_I]
        new_Isig = [str(Isig) for Isig in new_Isig]

        if '_refln.intensity_meas' in cif[cif_key]:
            cif[cif_key]['_refln.intensity_meas'] = new_I
            cif[cif_key]['_refln.intensity_sigma'] = new_Isig

        elif '_refln.F_meas' in cif[cif_key]:
            cif[cif_key]['_refln.F_meas'] = new_I
            cif[cif_key]['_refln.F_sigma'] = new_Isig
        # write the dict to file

        out_fname = fname.stem + f'_res{res}.cif'

        outfile = fname.parent / out_fname

        print(f'Saving to:{outfile}')
        out_file = open(str(outfile), 'w')
        out_file.write(cif.WriteOut())
        out_file.close()


def import_hklIF(cif):
    cif_key = cif.visible_keys[0]

    # read in the reflection indices and convert to int
    h_refl_orig = cif[cif_key]['_refln.index_h']
    h_refl_orig = np.array([int(h) for h in h_refl_orig])

    k_refl_orig = cif[cif_key]['_refln.index_k']
    k_refl_orig = np.array([int(k) for k in k_refl_orig])

    l_refl_orig = cif[cif_key]['_refln.index_l']
    l_refl_orig = np.array([int(l) for l in l_refl_orig])

    I_refl_orig = np.zeros(len(l_refl_orig), dtype=float)
    I_sig_orig = np.zeros(len(l_refl_orig), dtype=float)

    F_refl_orig = np.zeros(len(l_refl_orig), dtype=float)
    F_sig_orig = np.zeros(len(l_refl_orig), dtype=float)

    if '_refln.intensity_meas' in cif[cif_key]:
        print('Reading INTENSITY')
        I_refl_orig = cif[cif_key]['_refln.intensity_meas']
        I_sig_orig = cif[cif_key]['_refln.intensity_sigma']

        for i, (I, sig) in enumerate(zip(I_refl_orig, I_sig_orig)):
            if I == '?' or sig == '?':
                I_refl_orig[i] = '0'
                I_sig_orig[i] = '0'
            else:
                continue

        I_refl_orig = np.array([float(I) for I in I_refl_orig])
        I_sig_orig = np.array([float(sig) for sig in I_sig_orig])

        return h_refl_orig, k_refl_orig, l_refl_orig, I_refl_orig, I_sig_orig



    elif '_refln.F_meas_au' in cif[cif_key]:
        print('Reading STRUCTURE FACTORS')
        F_refl_orig = cif[cif_key]['_refln.F_meas_au']
        F_sig_orig = cif[cif_key]['_refln.F_meas_sigma_au']

        for i, (F, sig) in enumerate(zip(F_refl_orig, F_sig_orig)):
            if F == '?' or sig == '?':
                F_refl_orig[i] = '0'
                F_sig_orig[i] = '0'
            else:
                continue

        # Conversion from F to I (check implementation?)
        F_refl_orig = np.array([float(F) for F in F_refl_orig])
        F_sig_orig = np.array([float(sig) for sig in F_sig_orig])

        return h_refl_orig, k_refl_orig, l_refl_orig, F_refl_orig, F_sig_orig

    else:
        print('NO INTENSITY INFO FOUND IN CIF')
        exit()


if __name__ == "__main__":
    import sys
    import os

    #
    sys.path.append(str(Path(os.getcwd()).parent))

    base_path = Path('cifs/cell_size')

    comp_cif(base_path / '253l-sf_ave.cif', base_path / '254l-sf_ave.cif')
