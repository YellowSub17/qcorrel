import time
import CifFile
import numpy as np
import matplotlib.pyplot as plt



def cut_hkl(fname, resolutions, error=None):
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

    if '_refln.intensity_meas' in cif[cif_key]:
        I_refl_orig = cif[cif_key]['_refln.intensity_meas']
        I_sig_orig = cif[cif_key]['_refln.intensity_sigma']

        for i, (I, sig) in enumerate(zip(I_refl_orig, I_sig_orig)):
            if I=='?' or sig=='?':
                I_refl_orig[i] = '0'
                I_sig_orig[i] ='0'
            else:
                continue

        I_refl_orig = np.array([float(I) for I in I_refl_orig if I != '?'])
        I_sig_orig = np.array([float(sig) for sig in I_sig_orig])



    elif '_refln.F_meas_au' in cif[cif_key]:
        I_refl_orig = cif[cif_key]['_refln.F_meas_au']
        I_sig_orig = cif[cif_key]['_refln.F_meas_sigma_au']

        for i, (I, sig) in enumerate(zip(I_refl_orig, I_sig_orig)):
            if I=='?' or sig=='?':
                I_refl_orig[i] = '0'
                I_sig_orig[i] ='0'
            else:
                continue

        I_refl_orig = np.array([(float(I))**2 for I in I_refl_orig if I != '?'])
        I_sig_orig = np.array([(float(sig))**2 for sig in I_sig_orig if sig !='?'])



    else:
        print('NO INTENSITY OR F_MEAS FOUND IN CIF.')
        exit()




    for res in resolutions:

        print(f'Computing resolution: {res}')
        #int new lists
        new_h = []
        new_k = []
        new_l = []
        new_I = []
        new_sig = []

        q_max = 1.0/res

        h_bound = int(round(a_len*q_max))
        k_bound = int(round(b_len*q_max))
        l_bound = int(round(c_len*q_max))

        print(f'HKL BOUNDS: {h_bound}, {k_bound}, {l_bound}')

        for h, k, l, I, sig in zip(h_refl_orig, k_refl_orig, l_refl_orig,I_refl_orig,I_sig_orig):
            if abs(h) > h_bound or abs(k) > k_bound or abs(l) > l_bound:
                continue        #skip reflections greater than bounds

            # elif I=='?' or sig=='?':
            #     continue
            # else:



            if error==None:
                dI = 0


            elif error =='sig':
                dI = np.round(np.random.normal(0, float(sig)))



            #append indices less then bounds to the new lists
            new_h.append(h)
            new_k.append(k)
            new_l.append(l)
            new_I.append(I+dI)
            new_sig.append(sig)

        #convert the ints to str
        new_h = [str(h) for h in new_h]
        new_k = [str(k) for k in new_k]
        new_l = [str(l) for l in new_l]
        new_I = [str(I) for I in new_I]
        new_sig = [str(sig) for sig in new_sig]

        #save str list back into cif dict
        cif[cif_key]['_refln.index_h'] = new_h
        cif[cif_key]['_refln.index_k'] = new_k
        cif[cif_key]['_refln.index_l'] = new_l

        if '_refln.intensity_meas' in cif[cif_key]:
            cif[cif_key]['_refln.intensity_meas'] = new_I
            cif[cif_key]['_refln.intensity_sigma'] = new_sig

        elif '_refln.F_meas_au' in cif[cif_key]:
            cif[cif_key]['_refln.F_meas_au'] = new_I
            cif[cif_key]['_refln.F_meas_sigma_au'] = new_sig


        #write the dict to file
        out_file = open(f'{fname}_res{res}.cif', 'w')
        out_file.write(cif.WriteOut())
        out_file.close()


#
# def error_pm(fname, mag=None):
#



pdb_codes = ['CypA\\4yug','CypA\\4yuh', 'CypA\\4yui', 'CypA\\4yuj', 'CypA\\4yuk', 'CypA\\4yul', 'CypA\\4yum',
                'GFP\\2b3p', 'GFP\\5z6y', 'GFP\\6b9c', 'GFP\\2q6p','GFP\\4lqt']
pdb_codes = ['Lyso\\253l','Lyso\\254l']

#pdb_codes = ['Errors\\253l','Errors\\254l']





####INCREASEING ORDER!!
resolutions = [4,8,16]




## Utility to edit cif files for less reflections
for name in pdb_codes:
    start = time.time()
    cut_hkl(f'cifs\\{name}-sf',resolutions)
    end = time.time()
    print(f'Run time: {end-start} seconds.\n\n\n')
