import time
import CifFile
import numpy as np
import matplotlib.pyplot as plt


def cut_hkl(fname, resolutions):
    '''
    reduce the number of reflections in a CIF file fname
    if h,k or l are greater then bound, ignore this reflection
    remaining reflections are saved as outfname
    '''
    print(f'Reading Cif {fname}')
    cif = CifFile.ReadCif(f'{fname}.cif') # read the cif file

    cif_key = cif.visible_keys[0]

    a_len = float(cif[cif_key]['_cell.length_a'])
    b_len = float(cif[cif_key]['_cell.length_b'])
    c_len = float(cif[cif_key]['_cell.length_c'])

    #read in the reflection indices and convert to int
    h_refl_orig = cif[cif_key]['_refln.index_h']
    h_refl_orig = np.array([int(h) for h in h_refl_orig])

    k_refl_orig = cif[cif_key]['_refln.index_k']
    k_refl_orig = np.array([int(k) for k in k_refl_orig])

    l_refl_orig = cif[cif_key]['_refln.index_l']
    l_refl_orig = np.array([int(l) for l in l_refl_orig])




    for res in resolutions:

        print(f'Computing resolution: {res}')
        #int new lists
        new_h = []
        new_k = []
        new_l = []

        q_max = 1.0/res

        h_bound = int(round(a_len*q_max))
        k_bound = int(round(b_len*q_max))
        l_bound = int(round(c_len*q_max))

        #print(h_bound, k_bound, l_bound)

        print('Looping indices')
        #remove indices higher then bound
        for h, k, l in zip(h_refl_orig, k_refl_orig, l_refl_orig):
            if h > h_bound or k > k_bound or l > l_bound:
                continue        #skip reflections greater than bounds
            else:
                #append indices less then bounds to the new lists
                new_h.append(h)
                new_k.append(k)
                new_l.append(l)
        print('Finished Looping Indices')

        #convert the ints to str
        new_h = [str(h) for h in new_h]
        new_k = [str(k) for k in new_k]
        new_l = [str(l) for l in new_l]

        #save str list back into cif dict
        cif[cif_key]['_refln.index_h'] = new_h
        cif[cif_key]['_refln.index_k'] = new_k
        cif[cif_key]['_refln.index_l'] = new_l

        #write the dict to file
        out_file = open(f'{fname}_res{res}.cif', 'w')
        out_file.write(cif.WriteOut())
        out_file.close()





## Utility to edit cif files for less reflections

pdb_codes = ['4yug', '4yum','4yuh', '4yui', '4yuj', '4yuk', '4yul']

resolutions = [5]

for name in pdb_codes:


    start = time.time()
    cut_hkl(f'cifs\\{name}-sf',resolutions)
    end = time.time()
    print(f'Run time: {end-start} seconds.\n\n\n')
