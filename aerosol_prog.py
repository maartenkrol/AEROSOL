from numpy import *
from pylab import *
from IPython.html.widgets import interact, interactive
from IPython.html import widgets
from IPython.display import display
import os,sys,glob, copy
matplotlib.rcParams.update({'font.size': 12})

class aerosol_prog:
    
    def __init__(self):
        dfiles = glob.glob('DATA'+os.path.sep+'*.SPUV')
        dfiles.sort()

        self.slopes = zeros(6)
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
        roptions = ['raw counts','calibrated radiances','airmass']
        self.w2 = interact(self.aerosol,representation=roptions, wl_368nm=True,
             wl_501nm=False,wl_670nm=False,wl_780nm=False,wl_870nm=False,wl_940nm=False)

    
    def aerosol(self, representation = 'raw counts',wl_368nm=True,
            wl_501nm=False,wl_670nm=False,wl_780nm=False,wl_870nm=False,wl_940nm=False):

        flags = [wl_368nm,wl_501nm,wl_670nm,wl_780nm,wl_870nm,wl_940nm]
        labels = ['368nm','501nm','670nm','780nm','870nm','940nm']
        waves = ['368 nm','501 nm','670 nm','780 nm','870 nm','940 nm']
        colors = ['r','b','g','k','m','c']
        if representation == 'raw counts':
            f, ax = subplots()
            f.set_figheight(7)
            f.set_figwidth(10)
            for i,j in enumerate(flags):
                if j: 
                    ax.plot(self.spuv[:,1],self.spuv[:,i+3],label=labels[i])
            ax.set_title("Raw counts "+self.title)
            ax.set_xlabel('Day fraction (0-1)')
            ax.set_ylabel("Raw counts")
            ax.legend(loc='best')
        elif representation == 'calibrated radiances':
            f, ax = subplots()
            f.set_figheight(7)
            f.set_figwidth(10)
            for i,j in enumerate(flags):
                if j: 
                    ax.plot(self.cspuv[:,1],self.cspuv[:,i+3],label=labels[i])
            ax.set_title("Calibrated Radiances "+self.title)
            ax.set_xlabel('Day fraction (0-1)')
            ax.set_ylabel("Radiances (W m$^{-2}$ nm$^{-1}$)")
            ax.legend(loc='best')
        elif representation == 'airmass':
            f, ax = subplots(2,sharex=True)
            f.set_figheight(7)
            f.set_figwidth(10)
            ax[0].plot(self.xtime,self.m,'ro',label = "morning")
            ax[1].plot(self.xtime,self.sza,'ro',label = "morning")
            afternoon = where(self.xtime>0.5)
            ax[0].plot(self.xtime[afternoon],self.m[afternoon],'bo',label = "afternoon")
            ax[1].plot(self.xtime[afternoon],self.sza[afternoon],'bo',label = "afternoon")
            ax[0].set_title("Air mass $m$ and solar zenith angle during the day "+self.title)
            ax[1].set_xlabel('Day fraction (0-1)')
            ax[0].set_ylabel("$m$ (-)")
            ax[1].set_ylabel("solar zenith angle (degrees)")
            ax[0].set_ylim(0,20)
            ax[0].legend(loc='best')
