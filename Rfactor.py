
import imageio
import CifFile
import numpy as np
import matplotlib.pyplot as plt


def calc_R_from_cif(cif1, cif2,c=1):

    cif_key1 = cif1.visible_keys[0]
    h_refl_orig1 = cif1[cif_key1]['_refln.index_h']
    k_refl_orig1 = cif1[cif_key1]['_refln.index_k']
    l_refl_orig1 = cif1[cif_key1]['_refln.index_l']
    try:
        inten_meas_orig1 = cif1[cif_key1]['_refln.intensity_meas']
    except KeyError:
        inten_meas_orig1 = cif1[cif_key1]['_refln.F_meas_au']


    cif_key2 = cif2.visible_keys[0]
    h_refl_orig2 = cif2[cif_key2]['_refln.index_h']
    k_refl_orig2 = cif2[cif_key2]['_refln.index_k']
    l_refl_orig2 = cif2[cif_key2]['_refln.index_l']
    try:
        inten_meas_orig2 = cif2[cif_key2]['_refln.intensity_meas']
    except KeyError:
        inten_meas_orig2 = cif2[cif_key2]['_refln.F_meas_au']


    hkls1 = np.array([h_refl_orig1, k_refl_orig1, l_refl_orig1, inten_meas_orig1]).T
    hkls2 = np.array([h_refl_orig2, k_refl_orig2, l_refl_orig2, inten_meas_orig2]).T



    I1 = []
    I2 = []
    R=0
    missing= []
    for hkl1 in hkls1:
        match_found=False
        for hkl2 in hkls2:

            if np.equal(hkl1[:3,...].astype(int), hkl2[:3,...].astype(int)).all():
                match_found=True
                #print(hkl1, hkl2)
                I1.append(float(hkl1[3]))
                I2.append(float(hkl2[3]))


                R += (float(hkl1[3]) - float(c*hkl2[3]))/(float(hkl1[3])+1)



                continue

        if not match_found:
            print(f'No matching vector found for {hkl1}')
            missing.append(hkl1)


    return I1, I2,hkls1, hkls2, missing, R




def calc_R_from_correl(correl1, correl2, c=1):


    R = 0
    for q in range(correl1.shape[0]):
        print(q)
        for theta in range(correl1.shape[1]):
            if correl1[q,theta] ==0:
                continue
            else:

                R += (correl1[q,theta] - correl2[q,theta])/correl1[q,theta]


    return R




if __name__=='__main__':

    from pathlib import Path



    cif1_path = Path('cifs/Errors/253l-sf_res16.cif')
    cif2_path = Path('cifs/Lyso/253l-sf_res16.cif')

    cif1 = CifFile.ReadCif(str(cif1_path)) # read the cif file
    cif2 = CifFile.ReadCif(str(cif2_path)) # read the cif file



    I1, I2,hkls1, hkls2,missing, R2=  calc_R_from_cif(cif1, cif2)






    tiff1_path = Path('tiffs/4yug-sf_res8.tiff')

    tiff2_path = Path('tiffs/4yuh-sf_res8.tiff')

    tiff1 = imageio.imread(str(tiff1_path))
    tiff2 = imageio.imread(str(tiff2_path))


    R = calc_R_from_correl(tiff1, tiff2)

    plt.figure()
    plt.imshow(tiff1)

    plt.figure()
    plt.imshow(tiff2)


