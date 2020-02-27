import numpy as np
import symmetry as sym
from correlate import angle_between, calc_cif_qs, calc_q, get_cif_reflections, unit_vector
import multiprocessing as mp
import threading
from pathlib import Path

def correlate_worker(qs, q_primes, hist, nQ, nTheta, qmax, chunk_id):

    # chunk_hist = np.zeros(( nQ,nQ,nTheta ), dtype=int)
    for q_i, q in enumerate(qs):
        print(f'Correlating vector {q_i}/{len(qs)}. q={q[:3]}')
        q_ind = int(round(( (nQ-1) * (q[3] / qmax ))))

        for q_prime_i, q_prime in enumerate(q_primes):

            q_prime_ind = int(round(( (nQ-1) * (q_prime[3] / qmax ))))
            theta = np.degrees(angle_between(q[:3], q_prime[:3]))
            theta_ind = int(round((theta / 180.0) * (nTheta-1)))

            hist[q_ind, q_prime_ind, theta_ind] += 1
            # chunk_hist[q_ind, q_prime_ind, theta_ind] += 1

    # np.savetxt(f'thread_tmp\\chunk_{chunk_id}_th.txt', hist.flatten())





def full_correlate_threaded(cif, nQ, nTheta, qmax, nChunks):
    qs = calc_cif_qs(cif)

    correl_vec_indices = np.where(qs[:, 3] < qmax)[0]

    qs = qs[correl_vec_indices]

    hist = np.zeros((nQ, nQ, nTheta))

    remainder = len(qs) % nChunks
    chunk_size = (len(qs) - remainder) / nChunks

    chunk_bounds = [None] * nChunks

    for chunk_i in range(nChunks):
        b1 = chunk_size * chunk_i
        b2 = chunk_size * (chunk_i + 1)
        chunk_bounds[chunk_i]=[int(b1), int(b2)]

    chunk_bounds[-1][-1] += remainder


    threads = [None] *nChunks

    for chunk_i in range(nChunks):
        print(chunk_i)
        q_primes = qs[chunk_bounds[chunk_i][0]:chunk_bounds[chunk_i][1]]
        threads[chunk_i] = threading.Thread(target=correlate_worker, args=(qs,q_primes, hist, nQ,nTheta, qmax, chunk_i))
        threads[chunk_i].start()


    for chunk_i in range(nChunks):
        threads[chunk_i].join()


    return hist




if __name__ == '__main__':
    import CifFile
    import plot_and_process_utils as ppu
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import time

    prot_path = Path('cifs/Lyso/253l-sf_res8.cif')

    sf_cif = CifFile.ReadCif(str(prot_path))

    start = time.time()
    hist = full_correlate_threaded(sf_cif, 150,150, 0.2, 8)
    print(f'time taken {time.time() - start}')

    ppu.plot_map(np.sum(hist, axis=0))
    ppu.plot_map(np.sum(hist, axis=1))
    ppu.plot_map(np.sum(hist, axis=2))

    x,y,z = np.where(hist!=0)


    fig = plt.figure()
    Axes3D(fig).scatter(x,y,z)
    plt.show()

