
import threading
import numpy as np

arr = np.zeros( (10, 5))


def fill_thread(x, i):
    x[:,i] +=1


threadObj = threading.Thread(target=fill_thread, args=[arr, 0])
threadObj2 = threading.Thread(target=fill_thread, args=[arr, 1])

threadObj.start()
threadObj2.start()

print(arr)