import time
import CifFile
import numpy as np
from pathlib import Path


def add_err(fname, outfname, res=None, err_scale=1, tag=0):


    # fname = Path(f'cifs/{dir}/{pdb}-sf.cif')
    # if res!=None:
    #     fname = Path(f'cifs/{dir}/{pdb}-sf_res{res}.cif')

    print(f'Reading Cif {fname}')
    cif = CifFile.ReadCif(fname) # read the cif file

    cif_key = cif.visible_keys[0]

    num_ref = len(cif[cif_key]['_refln.index_h'])

    I_refl_orig= np.zeros(num_ref, dtype=float)
    I_sig_orig=np.zeros(num_ref, dtype=float)

    F_refl_orig = np.zeros(num_ref, dtype=float)
    F_sig_orig = np.zeros(num_ref, dtype=float)



    if '_refln.intensity_meas' in cif[cif_key]:
        I_refl_orig = cif[cif_key]['_refln.intensity_meas']
        I_sig_orig = cif[cif_key]['_refln.intensity_sigma']

        for i, (I, sig) in enumerate(zip(I_refl_orig, I_sig_orig)):
            if I=='?' or sig=='?':
                I_refl_orig[i] = '0'
                I_sig_orig[i] ='0'
            else:
                continue

        I_refl_orig = np.array([float(I) for I in I_refl_orig])
        I_sig_orig = np.array([float(sig) for sig in I_sig_orig])



    if '_refln.F_meas_au' in cif[cif_key]:
        F_refl_orig = cif[cif_key]['_refln.F_meas_au']
        F_sig_orig = cif[cif_key]['_refln.F_meas_sigma_au']

        for i, (F, sig) in enumerate(zip(F_refl_orig, F_sig_orig)):
            if F=='?' or sig=='?':
                F_refl_orig[i] = '0'
                F_sig_orig[i] ='0'
            else:
                continue

        # Conversion from F to I (check implementation?)
        # F_refl_orig = np.array([np.round((float(F))**2,3) for F in F_refl_orig])
        # F_sig_orig = np.array([np.round((float(sig))**2,3) for sig in F_sig_orig])
        F_refl_orig = np.array([float(F) for F in F_refl_orig])
        F_sig_orig = np.array([float(sig) for sig in F_sig_orig])




    #int new lists
    new_I= []
    new_Isig = []
    new_F = []
    new_Fsig = []



    for  I, Isig, F, Fsig in zip(I_refl_orig,I_sig_orig,F_refl_orig,F_sig_orig):


        #append indices less then bounds to the new lists

        new_I.append(I+np.random.normal(0, err_scale*Isig))
        new_Isig.append(Isig)
        new_F.append(F+np.random.normal(0, err_scale*Fsig))
        new_Fsig.append(Fsig)


    if np.max(I_refl_orig)!=0:
        new_I = [str(I) for I in new_I]
        new_Isig = [str(Isig) for Isig in new_Isig]
        cif[cif_key]['_refln.intensity_meas'] = new_I
        cif[cif_key]['_refln.intensity_sigma'] = new_Isig

    if np.max(F_refl_orig)!=0:
        new_F = [str(F) for F in new_F]
        new_Fsig = [str(Fsig) for Fsig in new_Fsig]
        cif[cif_key]['_refln.F_meas_au'] = new_F
        cif[cif_key]['_refln.F_meas_sigma_au'] = new_Fsig


    #write the dict to file

    out_file = open(f'{outfname}-sf_res{res}_err{err_scale}_tag{tag}.cif', 'w')
    out_file.write(cif.WriteOut())
    out_file.close()




def cut_hkl(fname, resolutions, outfname):
    '''
    reduce the number of reflections in a CIF file fname
    if h,k or l are greater then bound, ignore this reflection
    remaining reflections are saved as outfname
    '''
    print(f'Reading Cif {fname}')
    cif = CifFile.ReadCif(fname) # read the cif file

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

    I_refl_orig= np.zeros(len(l_refl_orig), dtype=float)
    I_sig_orig=np.zeros(len(l_refl_orig), dtype=float)

    F_refl_orig = np.zeros(len(l_refl_orig), dtype=float)
    F_sig_orig = np.zeros(len(l_refl_orig), dtype=float)


    if '_refln.intensity_meas' in cif[cif_key]:
        I_refl_orig = cif[cif_key]['_refln.intensity_meas']
        I_sig_orig = cif[cif_key]['_refln.intensity_sigma']

        for i, (I, sig) in enumerate(zip(I_refl_orig, I_sig_orig)):
            if I=='?' or sig=='?':
                I_refl_orig[i] = '0'
                I_sig_orig[i] ='0'
            else:
                continue

        I_refl_orig = np.array([float(I) for I in I_refl_orig])
        I_sig_orig = np.array([float(sig) for sig in I_sig_orig])



    if '_refln.F_meas_au' in cif[cif_key]:
        F_refl_orig = cif[cif_key]['_refln.F_meas_au']
        F_sig_orig = cif[cif_key]['_refln.F_meas_sigma_au']

        for i, (F, sig) in enumerate(zip(F_refl_orig, F_sig_orig)):
            if F=='?' or sig=='?':
                F_refl_orig[i] = '0'
                F_sig_orig[i] ='0'
            else:
                continue

        # Conversion from F to I (check implementation?)
        # F_refl_orig = np.array([np.round((float(F))**2,3) for F in F_refl_orig])
        # F_sig_orig = np.array([np.round((float(sig))**2,3) for sig in F_sig_orig])
        F_refl_orig = np.array([float(F) for F in F_refl_orig])
        F_sig_orig = np.array([float(sig) for sig in F_sig_orig])






    for res in resolutions:

        print(f'Computing resolution: {res}')
        #int new lists
        new_h = []
        new_k = []
        new_l = []
        new_I = []
        new_Isig = []
        new_F = []
        new_Fsig = []

        q_max = 1.0/res

        h_bound = int(round(a_len*q_max))
        k_bound = int(round(b_len*q_max))
        l_bound = int(round(c_len*q_max))

        print(f'HKL BOUNDS: {h_bound}, {k_bound}, {l_bound}')

        for h, k, l, I, Isig, F, Fsig in zip(h_refl_orig, k_refl_orig, l_refl_orig,I_refl_orig,I_sig_orig,F_refl_orig,F_sig_orig):
            if abs(h) > h_bound or abs(k) > k_bound or abs(l) > l_bound:
                continue        #skip reflections greater than bounds


            #append indices less then bounds to the new lists
            new_h.append(h)
            new_k.append(k)
            new_l.append(l)
            new_I.append(I)
            new_Isig.append(Isig)
            new_F.append(F)
            new_Fsig.append(Fsig)

        #convert the ints to str
        new_h = [str(h) for h in new_h]
        new_k = [str(k) for k in new_k]
        new_l = [str(l) for l in new_l]



        #save str list back into cif dict
        cif[cif_key]['_refln.index_h'] = new_h
        cif[cif_key]['_refln.index_k'] = new_k
        cif[cif_key]['_refln.index_l'] = new_l

        if np.max(I_refl_orig)!=0:
            new_I = [str(I) for I in new_I]
            new_Isig = [str(Isig) for Isig in new_Isig]
            cif[cif_key]['_refln.intensity_meas'] = new_I
            cif[cif_key]['_refln.intensity_sigma'] = new_Isig

        if np.max(F_refl_orig)!=0:
            new_F = [str(F) for F in new_F]
            new_Fsig = [str(Fsig) for Fsig in new_Fsig]
            cif[cif_key]['_refln.F_meas_au'] = new_F
            cif[cif_key]['_refln.F_meas_sigma_au'] = new_Fsig


        #write the dict to file
        out_file = open(f'{outfname}-sf_res{res}.cif', 'w')
        out_file.write(cif.WriteOut())
        out_file.close()



if __name__=="__main__":

    import sys
    import os
    #
    sys.path.append(str(Path(os.getcwd()).parent))

    from email_alert.alert import alert



    cif_path = Path('cifs/alpha/keratin')
    proteins = ['3tnu',
                '4zry',
                '6e2j',
                '6jfv',
                '6uui'
                ]





#   proteins = ['1al1',
#               '1cos',
#               '1mft',
#               '3hf0',
#               '102l',
#               '115l'
#               ]

    # proteins = [('CypA', '4yug'),
    #             ('CypA', '4yuh'),
    #             ('CypA', '4yui'),
    #             ('CypA', '4yuj'),
    #             ('CypA', '4yuk'),
    #             ('CypA', '4yul'),
    #             ('CypA', '4yum'),
    #
    #             ('GFP', '2b3p'),
    #             ('GFP', '2q6p'),
    #             ('GFP', '4lqt'),
    #             ('GFP', '5z6y'),
    #             ('GFP', '6b9c'),
    #
    #             ('Lyso', '253l'),
    #             ('Lyso', '254l')
    #             ]
    
    ####INCREASEING ORDER!!
    resolutions = [16,8,4]
    #
    #
    ## Utility to edit cif files for less reflections
    for protein in proteins:
        start = time.time()
        protien_path = cif_path / f'{protein}-sf.cif'
        out_path = cif_path / protein
        cut_hkl(str(protien_path),resolutions, str(out_path))
        end = time.time()
        print(f'Run time: {end-start} seconds.\n\n\n')

    #alert(sub='CUT HKL DONE')

    # for protein in proteins:
    #     for res in resolutions:
    #         for tag in range(2, 10):
    #             start = time.time()
    #             protien_path = cif_path / protein [0] / f'{protein[1]}-sf_res{res}.cif'
    #             out_path = cif_path / 'Errors' / protein [0]/ protein[1]
    #             add_err(str(protien_path), str(out_path), res, tag=tag)
    #             end = time.time()
    #             print(f'Run time: {end-start} seconds.\n\n\n')

    #alert(sub='ADD HKL DONE')
