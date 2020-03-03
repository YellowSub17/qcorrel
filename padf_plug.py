from pathlib import Path


if __name__== '__main__':
    import os
    import sys

    sys.path.append('/home/pat/rmit-onedrive/phd/python_projects')
    from email_alert.alert import alert

    

    qcorrel_dbins_path = Path('/home/pat/rmit-onedrive/phd/python_projects/Qcorrel/dbins/')
    padf_path = Path('/home/pat/rmit-onedrive/phd/python_projects/py3padf02/padf/')

    dbins = os.listdir(str(qcorrel_dbins_path))

    for dbin in dbins:
    #for i in range(1):
        #dbin = dbins[i]
        print(dbin)
        qcorrelation_path = qcorrel_dbins_path /dbin
    
        fname = qcorrelation_path.stem
    
        
        pdb_code, nQ, nTheta, qmax = fname.split('-')
    
        nQ = nQ[2:]
        nTheta = nTheta[2:]
        qmax = qmax[2:]
    
    
        config_file = open(str(padf_path /"config.txt"), 'w')
    
        config_file.write(f'correlationfile = {str(qcorrelation_path)}\n\n')
        config_file.write(f'outpath = {str(padf_path/"output")}\n\n')
    
        config_file.write(f'wavelength = 1e-10\n\n')
        
        config_file.write(f'tag = {fname}\n\n')
    
        
        config_file.write(f'nthq = {nTheta}\n\n')
    
        config_file.write(f'nq = {nQ}')
    
        config_file.write(f'nthr = {nTheta}\n\n')
    
        config_file.write(f'nr = {nQ}\n\n')
    
    
        config_file.write(f'nl = 20\n\n')
    
    
        config_file.write(f'qmax = {float(qmax)/1e-10}\n\n')
        
        config_file.write(f'rmax = {float(nQ)*(1e-10)/(2*float(qmax))}\n\n')
        
    
        config_file.close()
    
        #os.system(f'cd {padf_path}')
        os.system(f'cd {padf_path}')
        cmd = f'{padf_path}/padf {padf_path/"config.txt"}'
    
        os.system(cmd)
        stream = os.popen(f'rm {padf_path/"output"/"*bl*"}')
        stream = os.popen(f'rm {padf_path/"output"/"*r_vs_l*"}')
        
   
    alert(sub='CODE FINISHED', msg='THe code has finished running.')



