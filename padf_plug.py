from pathlib import Path


if __name__== '__main__':
    import os
    import sys

    sys.path.append('/home/patrick/rmit-onedrive/phd/python_projects')



    qcorrel_dbins_path = Path('/home/patrick/rmit-onedrive/phd/python_projects/Qcorrel/dbins/')
    padf_path = Path('/home/patrick/rmit-onedrive/phd/python_projects/py3padf02/padf/')




    dbins = []

#    dbins.append(qcorrel_dbins_path / '253l-sf_ave_qcorrel.dbin')
#    dbins.append(qcorrel_dbins_path / '254l-sf_ave_qcorrel.dbin')
#    dbins.append(qcorrel_dbins_path / '253l-sf_ave_comp254l-sf_ave_qcorrel.dbin')
#    dbins.append(qcorrel_dbins_path / '254l-sf_ave_comp253l-sf_ave_qcorrel.dbin')
#    qmax = [0.14,0.14,0.14,0.14]


#    dbins.append(qcorrel_dbins_path / '1al1-sf_qcorrel.dbin')
#    dbins.append(qcorrel_dbins_path / '1mft-sf_qcorrel.dbin')
#    dbins.append(qcorrel_dbins_path / '1cos-sf_qcorrel.dbin')
#    dbins.append(qcorrel_dbins_path / '4ggr-sf_qcorrel.dbin')
#    dbins.append(qcorrel_dbins_path / '5z6y-sf_qcorrel.dbin')
#    dbins.append(qcorrel_dbins_path / '4lqt-sf_qcorrel.dbin')
#    dbins.append(qcorrel_dbins_path / '2b3p-sf_qcorrel.dbin')
#    qmaxs = [0.3,0.39,0.35,0.35,0.15,0.3,0.15]


    dbins.append(qcorrel_dbins_path / '1al1-sf_qcorrel.dbin')
    dbins.append(qcorrel_dbins_path / '1mft-sf_qcorrel.dbin')
    dbins.append(qcorrel_dbins_path / '1cos-sf_qcorrel.dbin')
    dbins.append(qcorrel_dbins_path / '4ggr-sf_qcorrel.dbin')
    dbins.append(qcorrel_dbins_path / '4lqt-sf_qcorrel.dbin')
    qmaxs = [0.3,0.39,0.35,0.35,0.3]


    
    nQ =256
    nTheta = 360
    
    
    rmaxs = [50, 45, 40, 35,30,25,20]

    res = round(256/50)   #pix/A

    nRs = [round(rmax*res) for rmax in rmaxs]


    nl = 20
    

    for rmax, nR in zip(rmaxs, nRs):
        for dbin,qmax in zip(dbins, qmaxs):



            fname = dbin.stem


            config_file = open(str(padf_path /"config.txt"), 'w')

            config_file.write(f'correlationfile = {str(dbin)}\n\n')
            config_file.write(f'outpath = {str(padf_path/"output")}\n\n')

            config_file.write(f'wavelength = 1e-10\n\n')

            config_file.write(f'tag = {fname}_rmax{rmax}\n\n')


            config_file.write(f'nthq = {nTheta}\n\n')

            config_file.write(f'nq = {nQ}\n\n')

            config_file.write(f'nthr = {nTheta}\n\n')

            config_file.write(f'nr = {nR}\n\n')


            config_file.write(f'nl = {nl}\n\n')


            config_file.write(f'qmax = {float(qmax)/1e-10}\n\n')

            config_file.write(f'rmax = {rmax*1e-10}\n\n')



            config_file.close()

            os.system(f'cd {padf_path}')
            os.system(f'cd {padf_path}')
            cmd = f'{padf_path}/padf {padf_path/"config.txt"}'

            os.system(cmd)
            stream1 = os.popen(f'rm {padf_path/"output"/"*bl*"}')
            stream2 = os.popen(f'rm {padf_path/"output"/"*r_vs_l*"}')







