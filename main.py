from correlate import *
import CifFile
import plot_and_process_utils as ppu
import matplotlib.pyplot as plt

#protein_names = ['CypA\\4yug', 'CypA\\4yum','CypA\\4yuh', 'CypA\\4yui', 'CypA\\4yuj', 'CypA\\4yuk', 'CypA\\4yul']

#protein_names = ['GFP\\5z6y', 'GFP\\6b9c','GFP\\2q6p','GFP\\2b3p',]

protein_names = ['Lyso\\253l','Lyso\\254l']

# protein_names = ['Errors\\253l','Lyso\\253l']
# protein_names = ['GFP\\2q6p']


dqs = [0.001]

res_num = [16, 8]

# Read in Cif file
for name in protein_names:
    for res in res_num:

        print('Reading Cif')
        sf_cif = CifFile.ReadCif(f'cifs\\{name}-sf_res{res}.cif')
        for dq in dqs:

            print(f'STARTING NEW - {name} {res} {dq}')

            qmax=0.22
            correl, hist = correlate(sf_cif, dq,  0.5, q_bounds=[0,qmax])

            #correl, hist, qmax = correlate(sf_cif, dq,  0.5)


            correl = ppu.convolve_gaussian(correl, 3,3)



            #plt.plot(correl[-5,:], label=f'{res}')
            ppu.save_tiff(correl.astype(np.int32), f'correl_{name[:3]}{name[-4:]}_res{res}')




            ppu.plot_map(correl,title=f'correl {name}, res {res}, qmax {qmax}', extent=[0,180, qmax,0], save=f'correl_{name[:3]}{name[-4:]}_res{res}')




            print('\n\n\n')

