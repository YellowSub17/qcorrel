import numpy as np
import symmetry as sym
from correlate import angle_between, calc_cif_qs, calc_q, get_cif_reflections, unit_vector
import multiprocessing as mp
import concurrent.futures
import threading
from pathlib import Path

def correlate_worker(inputs):

    qs, q_primes, nQ, nTheta, qmax = inputs
    chunk_correl = np.zeros(( nQ,nQ,nTheta ), dtype=np.float32)
    chunk_hist = np.zeros(( nQ,nQ,nTheta ), dtype=np.float32)
    for q_i, q in enumerate(qs):
        print(f'Correlating vector {q_i}/{len(qs)}. q={q[:3]}')
        q_ind = int(round(( (nQ-1) * (q[3] / qmax ))))

        for q_prime_i, q_prime in enumerate(q_primes):

            q_prime_ind = int(round(( (nQ-1) * (q_prime[3] / qmax ))))
            theta = np.degrees(angle_between(q[:3], q_prime[:3]))
            theta_ind = int(round((theta / 180.0) * (nTheta-1)))

            chunk_correl[q_ind, q_prime_ind, theta_ind] += q[6]*q_prime[6]
            chunk_hist[q_ind, q_prime_ind, theta_ind] += 1

    return chunk_correl, chunk_hist

    # np.savetxt(f'thread_tmp\\chunk_{chunk_id}_th.txt', hist.flatten())





def full_correlate_threaded(cif, nQ, nTheta, qmax, nChunks):
    qs = calc_cif_qs(cif)

    correl_vec_indices = np.where(qs[:, 3] < qmax)[0]

    qs = qs[correl_vec_indices]

    hist = np.zeros((nQ, nQ, nTheta), dtype=np.float32)
    correl = np.zeros((nQ, nQ, nTheta), dtype=np.float32)

    remainder = len(qs) % nChunks
    chunk_size = (len(qs) - remainder) / nChunks



    all_args = []

    for chunk_i in range(nChunks):

        args =[]

        b1 = int(chunk_size * chunk_i)
        b2 = int(chunk_size * (chunk_i + 1))
        if chunk_i == nChunks-1:
            b2+=remainder
        args.append(qs)
        args.append(qs[b1:b2])
        args.append(nQ)
        args.append(nTheta)
        args.append(qmax)
        all_args.append(tuple(args))



    with concurrent.futures.ProcessPoolExecutor() as executor:

        results = [executor.submit(correlate_worker, args) for args in all_args]

        for f in concurrent.futures.as_completed(results):

            correl +=f.result()[0]
            hist +=f.result()[1]

    return correl, hist







if __name__ == '__main__':
    import CifFile
    import plot_and_process_utils as ppu
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import time

    prot_path = Path('cifs/Lyso/253l-sf_res8.cif')

    sf_cif = CifFile.ReadCif(str(prot_path))

    start = time.time()
    correl, hist = full_correlate_threaded(sf_cif, 150,360, 0.05, 4)
    print(f'time taken {time.time() - start}')


    ppu.plot_map(np.sum(hist, axis=0))
    ppu.save_tiff(np.sum(hist, axis=0), 'hist_no_conv_thread')
    convol_hist = ppu.convolve_3D_gaussian(hist, 0.5,0.5,0.5, 9)
    ppu.plot_map(np.sum(convol_hist, axis=0))
    ppu.save_tiff(np.sum(convol_hist, axis=0), 'hist_conv_thread')

    ppu.plot_map(np.sum(correl, axis=0))
    ppu.save_tiff(np.sum(correl, axis=0), 'correl_no_conv_thread')
    convol_correl = ppu.convolve_3D_gaussian(correl, 0.5,0.5,0.5, 9)
    ppu.plot_map(np.sum(convol_correl, axis=0))
    ppu.save_tiff(np.sum(convol_correl, axis=0), 'correl_conv_thread')


    # ppu.plot_map(np.sum(correl, axis=1))
    # ppu.plot_map(np.sum(correl, axis=2))
    #
    # x,y,z = np.where(correl!=0)
    #
    #
    # fig = plt.figure()
    # Axes3D(fig).scatter(x,y,z)
    # plt.show()
    #
