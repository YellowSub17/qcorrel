import plot_and_process_utils as ppu
import numpy as np
import matplotlib.pyplot as plt


class displayPADF():
    
    def __init__(self, gbf=False,abf=False,rbf=False,tf=False,
                     rf=False, fname='None',res=0.1, wid=15,widr=20,w=4,maxrad=20,
                     rmax=100,rmin=25,thmax=180,cmin=None,cmax=None,
                     xp=[],yp=[], arc_dists =[],title='', fig_size=None, dpi=100, save_fname=None):
        self.gaussBlurFlag = gbf #False
        self.angularBlurFlag = abf
        self.radialBlurFlag = rbf
        self.tflag = tf
        self.rflag = rf
        self.title = title  
        self.arc_dists = arc_dists
        self.res = res
        self.fname = fname
        self.dpi = dpi
        self.wid = wid
        self.widr = widr
        self.w = w
        self.save_fname=save_fname
        self.maxrad = maxrad
        self.rmax = rmax
        self.rmin = rmin
        self.thmax = thmax
        
        self.cmin = cmin
        self.cmax = cmax
        self.xpnt=xp
        self.ypnt=yp
        
        self.fig_size = fig_size

    def loadAndPlot(self):
        
        data = np.loadtxt( self.fname ) #.transpose()
        
        if self.tflag == True:
            data = data.transpose()
            
        s = data.shape
        #self.wid = 15
        #self.widr = 20
        
        th = np.outer(np.ones(s[0]), np.arange(s[1]))*np.pi/s[1]
        r = np.outer(np.arange(s[0]), np.ones(s[1]))
        
        sn = np.sin(th)
        
        
        #data *= 1.0/(sn + 1e-3)
        #data *= sn
        
        if self.rflag==True:
            data[10:,:] *= 1.0/(r[10:,:]**2)
        #plt.imshow(sn)
        
        #plt.figure()
        #for i in np.arange(s[0]):
        #    data[i,:] += -np.average(data[i,:])
        
        
        datafilt = data*0.0
        
        
        #for i in np.arange(s[0]):
        #    f = np.fft.ifft(data[i,:])
        #    f[wid:-wid] = 0.0
        #    datafilt[i,:] = np.real(np.fft.ifft(f))
        
        
        x = np.outer(np.arange(s[0]), np.ones(s[1])) - s[0]/2
        y = np.outer(np.ones(s[0]), np.arange(s[1])) - s[1]/2
        rsq = x*x + y*y
        #self.w = 4
        g = np.exp( - rsq/(self.w*self.w))
        g = np.roll(np.roll(g, -s[0]//2, 0), -s[1]//2, 1)
        fg = np.fft.fft2(g)
        
        f = np.fft.fft2(data)
        
        if self.angularBlurFlag == True:
            f[:,wid:-wid] = 0.0
        
        if self.radialBlurFlag == True:    
            f[widr:-widr,:] = 0.0
        
        if self.gaussBlurFlag == True:
            f *= fg
        datafilt = np.real(np.fft.ifft2(f))
        #maxrad = 60.0
        #rmin = 25
        r0 = self.maxrad*(self.rmin/s[0])
        #rmax = 100
        r1 = self.maxrad*(self.rmax/s[0])
        
        #thmax = s[1]//2
        ppu.plot_map(datafilt[self.rmin:self.rmax,:self.thmax],title=self.title,fig_size=self.fig_size,dpi=self.dpi, cb=False, cmap='viridis', xlabel='$\\theta$ / $ ^{\circ}$', ylabel='$r_1=r_2$ / $\AA$', origin="lower", extent=[0,180,r0,r1])
        #plt.imshow(datafilt[self.rmin:self.rmax,:self.thmax], origin="lower", extent=[0,180,r0,r1], aspect=9)

        
        plt.clim([self.cmin,self.cmax])
        
        clr = 'r'
        mkr = 'x'

        cols = plt.cm.gist_rainbow(np.linspace(0,1,len(self.xpnt)))
        for x,y, col in zip(self.xpnt,self.ypnt, cols):
#            plt.plot(x, y,'o',  ms=6, markeredgewidth=1.5, markerfacecolor=col, markeredgecolor=(0,0,0,1), label=f'{x} {y}')
            plt.plot(x, y,'o',  ms=6, markeredgewidth=1.5, markerfacecolor=(0.8,0.8,0.8,1), markeredgecolor=(0,0,0,1))
        
        plt.plot(65, 10, 'r.', ms=7)

        cols = plt.cm.autumn(np.linspace(0,0.9, len(self.arc_dists)))
        for arc_dist,col in zip(self.arc_dists, cols):
            t_space = np.linspace(1, 90, 200)
            plt.plot( t_space, (arc_dist/2)/np.sin(np.radians(t_space)/2),linestyle='dashed', color=col, label=f'${arc_dist}\AA$', dashes=(5,10))
        
        plt.xlim(0, 180)
        plt.ylim(self.rmin*self.res, 20)
        plt.legend() 
        plt.draw()

        if self.save_fname != None:
            plt.savefig(self.save_fname)
            




#arc_dists = [1.2, 1.3,1.5,2.5,6, 11]
#arc_dists = [ 1.3,2.2,6.2,8.5,11, 14]
arc_dists = [ 1.3,2.2,6.2,11]
# Jack's 1al1 model
fname = "1al1_ex_rcor.dat"

#xpnt = [65,65, 72, 65, 90, 90, 27, 45, 90, 23, 50, 90, 90, 75, 55, 80]
#ypnt = [8,10, 11.8, 5.8, 4.5, 6.8, 5, 9, 9, 7, 7, 15, 13, 13.5, 12, 10.75]


xpnt = [65,65, 72, 65, 90, 90, 27, 45, 90, 19, 53, 90, 90, 78, 55, 80, 64.2]
ypnt = [8,10, 11.8, 5.8, 4.1, 7.1, 5, 9, 9, 7, 7, 15, 13, 14.2, 12, 10.75,13.3]
fig_size = (6,4)
dpi=150

al_ex = displayPADF(gbf=True,abf=False,rbf=False,tf=True,rf=True,
                     fname=fname,wid=15,widr=20,w=5,maxrad=20,
                     rmax=201,rmin=30,thmax=89,cmin=-1e-3,cmax=1e-2,
                     xp=xpnt,yp=ypnt, arc_dists=arc_dists,dpi=dpi, fig_size=fig_size, save_fname='1al1_j_blur.png')
al_ex.loadAndPlot()

# Patrick's 1al1 models
fname = "1al1_padf_r1r2.dat"

al = displayPADF(gbf=False,res=0.2, abf=False,rbf=False,tf=False,rf=False,
                     fname=fname,wid=15,widr=20,w=4,maxrad=60,
                     rmax=100,rmin=15,thmax=180,cmin=-1e72,cmax=1e73,
                     xp=xpnt,yp=ypnt, arc_dists=arc_dists, dpi=dpi, fig_size=fig_size, save_fname='1al1_p_padf.png')
al.loadAndPlot()

# Jack's 1al1 model
fname = "1al1_ex_rcor.dat"


al_ex2 = displayPADF(gbf=True,abf=False,rbf=False,tf=True,rf=True,
                     fname=fname,wid=15,widr=20,w=0.25,maxrad=20,
                     rmax=201,rmin=30,thmax=89,cmin=-1e-7,cmax=1e-6,
                     xp=xpnt,yp=ypnt, arc_dists=arc_dists,dpi=dpi, fig_size=fig_size,save_fname='1al1_j_norm.png')

al_ex2.loadAndPlot()

plt.show()
