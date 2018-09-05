from numpy import *
from pylab import *
from IPython.html.widgets import interact, interactive
from IPython.html import widgets
from IPython.display import display
import os,sys,glob, copy
matplotlib.rcParams.update({'font.size': 12})

class langley_prog:
    
    def __init__(self):
        dfiles = glob.glob('DATA'+os.path.sep+'*.SPUV')
        dfiles.sort()

        self.slopes = zeros(6)
        self.dayfraction = ''
        # calibration constants
        self.spa = [0.0, 2.412917, 1.0417198/1.4, 2.60924/2.5, 0.92835/0.85, 1.32079/0.85, 9.28147]
        self.spb = [0.0, 0.001200, -0.0009700, -0.00040, 0.00050, -0.00197, 0.00030]
        self.spc =  1000.0
        self.spd =     0.0 
        self.w1 = interact(self.file_select,filename=dfiles)


        
    def file_select(self,filename=''):
        f = open(filename,'r')
        self.title = filename
        lines = f.readlines()
        spuv = []
        for line in lines[18:]:
            spuv.append([float(val) for val in line.split()[0:]]) 
        f.close()
        self.spuv=array(spuv)
        self.cspuv = copy.deepcopy(self.spuv)
        for i in range(0,6):
            self.cspuv[:,i+3] = (self.spuv[:,i+3]-(self.spb[i+1]*self.spc-self.spd))/(self.spa[i+1]*self.spc)
        self.cspuv = where(self.cspuv<0,1e-10,self.cspuv)
        self.costheta = self.spuv[:,2]
        self.costheta = where(self.costheta<0,1e-10,self.costheta)
        self.m = 1./self.costheta
        self.xtime = self.spuv[:,1]
        self.sza = arccos(self.costheta)*180/pi
        # remove possible existing widget before starting a new one
        try:
            self.w2.widget.close()
        except:
            None
        # setup widget:
        roptions = ['slopes','Langley slopes']
        self.w2 = interact(self.langley,representation=roptions, wl_368nm=True,
             wl_501nm=True,wl_670nm=True,wl_780nm=True,wl_870nm=True,wl_940nm=True,
             dayfraction = ['morning','afternoon'], mmin = (1,10,0.1), mmax=(1,10,0.1))

    
    def langley(self, representation = 'slopes',wl_368nm=True,
                wl_501nm=True,wl_670nm=True,wl_780nm=True,wl_870nm=True,wl_940nm=True,
               dayfraction = 'morning', mmin = 4, mmax=6):
        flags = [wl_368nm,wl_501nm,wl_670nm,wl_780nm,wl_870nm,wl_940nm]
        labels = ['368nm','501nm','670nm','780nm','870nm','940nm']
        waves = ['368 nm','501 nm','670 nm','780 nm','870 nm','940 nm']
        colors = ['r','b','g','k','m','c']
        if representation.endswith('slopes'):
            f, ax = subplots()
            f.set_figheight(7)
            f.set_figwidth(10)
            self.dayfraction = dayfraction
            # select points depending on settings:
            if mmax < mmax:
                mt = mmin
                mmin = mmax
                mmax = mt
            logical_mask = zeros((len(self.m),2), bool)
            logical_mask[:,0] = ma.masked_inside(self.m,mmin,mmax).mask
            if dayfraction == 'morning':
                logical_mask[:,1] = ma.masked_less(self.xtime,0.5).mask
            else:
                logical_mask[:,1] = ma.masked_greater(self.xtime,0.5).mask
            selected = where(all(logical_mask,axis=1))[0]
            for i,j in enumerate(flags):
                if j:
                    if representation == 'slopes':
                        ax.plot(self.m[selected],log(self.cspuv[selected,i+3]),label=labels[i])
                    else:
                        xp = self.m[selected]
                        yp = log(self.cspuv[selected,i+3])
                        res = polyfit(xp,yp,1)
                        ax.plot(xp,yp,'o'+colors[i],label=labels[i])
                        if (res[1] > 0):
                            txt = 'y = %6.3fx + %6.3f'%(res[0],res[1])
                        else:    
                            txt = 'y = %6.3fx - %6.3f'%(res[0],abs(res[1]))
                        ax.plot(xp,xp*res[0]+res[1],'-'+colors[i],label=txt)
                        self.slopes[i]=res[0]
            ax.set_title("Langley plot "+self.title)
            ax.set_xlabel('Air mass (-)')
            ax.set_ylabel("$ln{L_{\lambda}}$")
            ax.legend(loc='best')
