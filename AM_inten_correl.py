import numpy as np
import matplotlib.pyplot as plt

## make a 2D array with a gaussian
def make_gaussian(nx, ny, rad=None, rady=-1., cenx=None, ceny=None, invert=0, norm=False ):
    #set defaults
    if rad is None: rad = np.min(nx,ny)/2
    if cenx is None: cenx = nx/2
    if ceny is None: ceny = ny/2
    radsq = rad**2
    if rady == -1.:
        radysq = radsq
    else:
        radysq = rady**2

    # define the circle
    x = np.outer(np.arange(0-nx/2,nx-nx/2,1),np.ones(ny))
    #print x.size, x.shape
    y = np.outer(np.ones(nx),np.arange(0-ny/2,ny-ny/2,1))
    #print y.size, y.shape

    a = np.zeros([nx,ny])
    a = np.exp(-(x**2)/radsq  - ( y**2)/radysq)
    a[ int(nx/2), int(ny/2) ] = 1.0

    a = array_shift(a,int(cenx-nx/2),int(ceny-ny/2))

    # normalise if required
    if norm == True: a *= 1./np.sum(a)

    return a
#end make_gaussian


# shift - a 2D version of numpy's roll
def array_shift(array,xshift=0,yshift=0):
    array = np.roll(array,xshift,0)
    array = np.roll(array,yshift,1)
    return array

def convolve_gaussian( image, rad=3, rady=1):

    c = make_gaussian( image.shape[0], image.shape[1], rad, rady, cenx=0, ceny=0, norm=True )
    fc = np.fft.fft2( c )
    fimage = np.fft.fft2( image )
    output = np.real(np.fft.ifft2( np.conjugate(fc)*fimage ))
    return output


fname = "cifs\\1zux-sf_less20noCom.cif"
#load array of h,k,l,Int.
data = np.loadtxt( fname, usecols=(3,4,5,7))

# #init image array
# array = np.zeros( (93,93) )
#
# #for every reflection
# for i in np.arange(data.shape[0]):
#
#     m = int(data[i,0] + 46)
#     n = int(data[i,1] + 46)
#
#     print(m,n)
#     array[m,n] = 1
#
# array += array[:,::-1]
#array += np.rot90(array)

#
# apply rotational symmetries
#


# init a list of lists
alist = []

#for every reflection h,k,l, inten
for i in np.arange(data.shape[0]):


    #init a single point + 4 sym
    a = []

    #get inital hkl, Inten
    b = [data[i,0],data[i,1],data[i,2],data[i,3]]
    #add to a
    a.append(b)


    # calc four fold rotation and that point to a


    a.append( [ -b[1], b[0], b[2], b[3] ] )
    a.append( [ -b[0], -b[1], b[2], b[3] ] )
    a.append( [ b[1], -b[0], b[2], b[3] ] )

    # two fold
    c = []
    print(len(a))
    for b in a:
        a.append( [-b[0],b[1],-b[2],b[3]] )
    a += c
    print('x')

    # second two fold
    c = []
    for b in a:
        c.append( [b[0],-b[1],-b[2],b[3]] )
    a += c

    alist += a


#reflections and intensity
alist = np.array(alist)

array = np.zeros( (93,93) )

for i in np.arange(alist.shape[0]):

    m =int( alist[i,0] + 46)
    n =int( alist[i,1] + 46)
    array[m,n] = alist[i,3]

#
# crystal parameters
#
a, b, c = 77.050,   77.050,   37.210
ast, bst, cst = 1/a, 1/b, 1/c

na, nb, nc = 5, 5, 5
qmax = np.max([na*ast, nb*bst, nc*cst])
nq = 16*np.max([na,nb,nc])
qstep = qmax / float(nq)

#
#make resolution shell lists
#
rlist = []
for i in np.arange(nq):
    rlist.append([])

for i in np.arange(len(alist)):
    ia, ib, ic = alist[i][0], alist[i][1], alist[i][2]
    if (ia>na) or (ib>nb) or (ic>nc):
        continue

    q = np.sqrt( (ia*ast)**2 + (ib*bst)**2 + (ic*cst)**2)
    iq = int(q / qstep)
    if (iq>=0) and (iq < nq):
        rlist[iq].append( alist[i] )

#
# correlate the q-vectors
#
nth = 360
icorr = np.zeros( (nq, nth) )
print(len(rlist))

for iq in np.arange(nq):
    print("Correlating:", iq+1, "/", nq, len(rlist[iq]))
    for i in np.arange(len(rlist[iq])):
        ia, ib, ic = rlist[iq][i][0], rlist[iq][i][1], rlist[iq][i][2]
        q1v = np.array([ia*ast, ib*bst, ic*cst])
        q1 = np.sqrt( (ia*ast)**2 + (ib*bst)**2 + (ic*cst)**2)
        for j in np.arange(len(rlist[iq])):
            if i==j:
                continue

            ia, ib, ic = rlist[iq][j][0], rlist[iq][j][1], rlist[iq][j][2]
            q2v = np.array([ia*ast, ib*bst, ic*cst])
            q2 = np.sqrt( (ia*ast)**2 + (ib*bst)**2 + (ic*cst)**2)

            if (q1>0) and (q2>0):
                arg = np.sum(q1v*q2v)/(q1*q2)
                if arg > 1.0:
                    arg = 1.0
                if arg < -1.0:
                    arg = -1.0
                theta = np.arccos( arg )
            else:
                theta = 0.0

            ith = int( (theta/np.pi) * nth )
            if (ith>=0) and (ith<nth):
                icorr[iq,ith] += 1

icorr += icorr[:,::-1]


#icorr[:,nth/2-3:nth/2+3] = 0.0

r = np.arange(nq)
th = np.arange(nth)*np.pi/float(nth)
sth = np.sin(th)
#filter = np.outer( np.ones(nq), sth**2 + 0.01 )
#filter = np.outer( r*r, np.ones(nth) )
filter = np.outer( r*r, np.abs(sth) + 0.001 )
ifilt = np.where( filter > 0.0)
icorr[ifilt] *= 1.0/filter[ifilt]

icorr = convolve_gaussian(icorr, 8, 8 )

print(np.max(icorr[:4,-4]))
plt.imshow(icorr[10:,4:-4], origin='lower')
plt.clim([0.0,0.015])

plt.figure()
plt.plot(icorr[75,:])
plt.draw()
plt.show()



