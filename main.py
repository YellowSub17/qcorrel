from correlate import *
import CifFile
import plot_and_process_utils as ppu

protein_names = ['4yug','4yuh', '4yui', '4yuj', '4yuk', '4yul', '4yum']

dqs = [0.001]

res_num = [5]

# Read in Cif file
for name in protein_names:
    for res in res_num:

        print('Reading Cif')
        sf_cif = CifFile.ReadCif(f'cifs\\{name}-sf_res{res}.cif')
        for dq in dqs:

            print(f'STARTING NEW - {name} {res} {dq}')


            correl, hist, qmax = correlate(sf_cif, dq,  0.5)

            correl = ppu.convolve_gaussian(correl, 10,10)

            ppu.save_tiff(correl.astype(np.int32), f'correl_{name}_res{res}')

            ppu.plot_map(correl,title=f'correl {name}, res {res}', save=f'correl_{name}_res{res}', extent=[0,180, qmax,0])




            print('\n\n\n')

