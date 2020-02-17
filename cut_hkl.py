import time
import CifFile
import numpy as np
import matplotlib.pyplot as plt


def cut_hkl(fname, bound, outfname):
    '''
    reduce the number of reflections in a CIF file fname
    if h,k or l are greater then bound, ignore this reflection
    remaining reflections are saved as outfname
    '''
    print('Reading Cif File...')
    sf_cif = CifFile.ReadCif(fname) # read the cif file
    print('Done')


    cif_key = sf_cif.visible_keys[0]

    #read in the reflection indices and convert to int
    h_refl = sf_cif[cif_key]['_refln.index_h']
    h_refl = np.array([int(h) for h in h_refl])

    k_refl = sf_cif[cif_key]['_refln.index_k']
    k_refl = np.array([int(k) for k in k_refl])

    l_refl = sf_cif[cif_key]['_refln.index_l']
    l_refl = np.array([int(l) for l in l_refl])


    #int new lists
    new_h = []
    new_k = []
    new_l = []


    print('Looping indices')
    #remove indices higher then bound
    for h, k, l in zip(h_refl, k_refl, l_refl):
        if h > bound or k > bound or l > bound:
            continue        #skip reflections greater than bounds
        else:
            #append indices less then bounds to the new lists
            new_h.append(h)
            new_k.append(k)
            new_l.append(l)
    print('Finished Looping Indices')

    #convert the ints to str
    h_refl = [str(h) for h in new_h]
    k_refl = [str(k) for k in new_k]
    l_refl = [str(l) for l in new_l]

    #save str list back into cif dict
    sf_cif[cif_key]['_refln.index_h'] = h_refl
    sf_cif[cif_key]['_refln.index_k'] = k_refl
    sf_cif[cif_key]['_refln.index_l'] = l_refl

    #write the dict to file
    out_file = open(outfname, 'w')
    out_file.write(sf_cif.WriteOut())
    out_file.close()





### Utility to edit cif files for less reflections

pdb_codes = ['9ins','5fya','5hul', '5mj0']

for name in pdb_codes:
    for num in range(5,11,5):

        start = time.time()
        print(f'Removing {num} and above hkl indices on {name}')
        cut_hkl(f'cifs\\{name}-sf.cif',num, f'cifs\\{name}-sf_less{num}.cif')
        end = time.time()
        print(f'Run time: {end-start} seconds.')
