from multiprocessing.pool import ThreadPool
import numpy as np

nThreads = 2
x = np.zeros((5,3))


def thread_fn(*args):


    x +=n


inputs = ((x, 0), (x,3), (x, 5))


pool = ThreadPool(processes=nThreads)

async_res = pool.apply_async(thread_fn, inputs)

val = async_res.get()


print(val)
