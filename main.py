from correlate import *
import CifFile
import plot_and_process_utils as ppu
import matplotlib.pyplot as plt
from pathlib import Path


cif_path = Path('cifs')


proteins=[('CypA', '4yug'),]



res_num = [16, 8]

# Read in Cif file
for protein in proteins:
    for res in res_num:

        fname = f'{protein[1]}-sf_res{res}.cif'
        protein_path = cif_path /protein[0] /fname

        print(protein_path)


        print('Reading Cif')
        sf_cif = CifFile.ReadCif(str(protein_path))

        
        print(f'STARTING NEW - {fname} {res}')

        qmax=0.22
        correl, hist = correlate(sf_cif, 0.001,  0.5, q_bounds=[0,qmax])

        #correl, hist, qmax = correlate(sf_cif, dq,  0.5)

        correl = ppu.convolve_gaussian(correl, 3,3)

        ppu.save_tiff(correl.astype(np.int32),f'{protein_path.stem}.tiff')


        ppu.plot_map(correl,title=f'correl {fname}, res {res}, qmax {qmax}', extent=[0,180, qmax,0], save=f'{protein_path.stem}.tiff')



        print('\n\n\n')

