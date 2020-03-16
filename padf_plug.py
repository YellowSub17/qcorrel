from pathlib import Path


if __name__== '__main__':
    import os
    import sys

    sys.path.append('/home/pat/rmit-onedrive/phd/python_projects')
    from email_alert.alert import alert



    qcorrel_dbins_path = Path('/home/pat/rmit-onedrive/phd/python_projects/Qcorrel/dbins/')
    padf_path = Path('/home/pat/rmit-onedrive/phd/python_projects/py3padf02/padf/')





    dbins = []
    dbins.append(qcorrel_dbins_path / '1cos-sf_qcorrel.dbin')


    dbins.append(qcorrel_dbins_path / '1mft-sf_qcorrel.dbin')


    #dbins.append(qcorrel_dbins_path / '4osd-sf_res4_qcorrel.dbin')



    #dbins.append(qcorrel_dbins_path /'6jfv-sf_res4_qcorrel.dbin')


    #dbins.append(qcorrel_dbins_path/ '6q5j-sf_res4_qcorrel.dbin')

    #dbins.append(qcorrel_dbins_path / '254l-sf_qcorrel.dbin')

    #dbins.append(qcorrel_dbins_path / '253l-sf_qcorrel.dbin')

    nQ = 150
    nTheta = 360
    qmaxs = [0.3, 0.3]#,  0.14,0.14,0.14]
    rmaxs = [20]*2


    nls = [20]
    try:
        for nl in nls:
            for dbin, qmax, rmax in zip(dbins, qmaxs,rmaxs):



                fname = dbin.stem


                config_file = open(str(padf_path /"config.txt"), 'w')

                config_file.write(f'correlationfile = {str(dbin)}\n\n')
                config_file.write(f'outpath = {str(padf_path/"output")}\n\n')

                config_file.write(f'wavelength = 1e-10\n\n')

                config_file.write(f'tag = {fname}nl{nl}\n\n')


                config_file.write(f'nthq = {nTheta}\n\n')

                config_file.write(f'nq = {nQ}')

                config_file.write(f'nthr = {nTheta}\n\n')

                config_file.write(f'nr = {nQ}\n\n')


                config_file.write(f'nl = {nl}\n\n')


                config_file.write(f'qmax = {float(qmax)/1e-10}\n\n')

                config_file.write(f'rmax = {rmax*1e-10}\n\n')

                #config_file.write(f'rmax = {float(nQ)*(1e-10)/(2*float(qmax))}\n\n')


                config_file.close()

                os.system(f'cd {padf_path}')
                os.system(f'cd {padf_path}')
                cmd = f'{padf_path}/padf {padf_path/"config.txt"}'

                os.system(cmd)
                stream1 = os.popen(f'rm {padf_path/"output"/"*bl*"}')
                stream2 = os.popen(f'rm {padf_path/"output"/"*r_vs_l*"}')

                #alert(msg=f'Finished padf calc: {str(dbin.stem)}')




            alert(msg='PADF correlation finished')

    except:
        print('err')
        #alert(msg='ERROR')

