import numpy as np
import matplotlib.pyplot as plt

class FootFall:
    def __init__(self, fs):
        self.set_fs(fs)
        self.cal_force_para()

    def set_fs(self, fs):
        self.fs = fs                                               # walking frequency
        self.vel = 1.67 * fs ** 2-4.83 * fs+4.5                    # walking velocity
        self.steplength = self.vel/self.fs                         # single step (two footprints) length
        self.timestep = 5.8531-6.9773*fs+3.2242*fs**2-0.515*fs**3  # single footprint duration
        self.timeoverlap = self.timestep - 1.0/self.fs             # overlap time of two footprint

    def cal_step_n(self, startoffset, distance):
        self.nstep = floor((distance-startoffset)/self.steplength-0.5)
        self.ttime = self.nstep/self.fs + self.timeoverlap

    def cal_force(self,T_NOW):

        if T_NOW <=1.0/self.fs+self.timeoverlap:
            self.force = 746*((self.k1*T_NOW**1)+(self.k2*T_NOW**2) +
                          (self.k3*T_NOW**3)+(self.k4*T_NOW**4)+(self.k5*T_NOW**5) +
                          (self.k6*T_NOW**6)+(self.k7*T_NOW**7)+(self.k8*T_NOW**8))
        else:
            self.force = 0

    def cal_force_para(self):
        if self.fs <= 1.75:
            self.k1 = -8.0 * self.fs + 38
            self.k2 = 376 * self.fs - 844
            self.k3 = -2804 * self.fs + 6025
            self.k4 = 6308 * self.fs - 16573
            self.k5 = 1732 * self.fs + 13619
            self.k6 = -24648 * self.fs + 16045
            self.k7 = 31836 * self.fs - 33614
            self.k8 = -12948 * self.fs + 15532

        elif self.fs <= 2:
            self.k1 = 24 * self.fs - 18
            self.k2 = -404 * self.fs + 521
            self.k3 = 4224 * self.fs - 6274
            self.k4 = -29144 * self.fs + 45468
            self.k5 = 109976 * self.fs - 175808
            self.k6 = -217424 * self.fs + 353403
            self.k7 = 212776 * self.fs - 350259
            self.k8 = -81572 * self.fs + 135624
        else:
            self.k1 = 75 * self.fs - 120
            self.k2 = -1720 * self.fs + 3153
            self.k3 = 17055 * self.fs - 31936
            self.k4 = -94265 * self.fs + 175710
            self.k5 = 298940 * self.fs - 553736
            self.k6 = -529390 * self.fs + 977335
            self.k7 = 481665 * self.fs - 888037
            self.k8 = -174265 * self.fs + 321008

class WalkDoubleStep:
    """this is the class to define the single step"""
    def __init__(self,initialpos, fs,xrange, yrange, mode='heel'):
        '''always assume left foot go before the right'''
        self.lf = FootFall(fs)
        self.rf = FootFall(fs)
        self.fs = float(fs)
        self.dT = float((self.rf.timestep-self.rf.timeoverlap)*2.0)     # time distance between two left steps
        self.ipos = initialpos
        self.xrange = xrange
        self.yrange = yrange
        self.mode = mode
        self.stepvel = xrange/self.lf.timestep

    def update(self,ct):
        self.get_left_status(ct)
        self.get_right_status(ct)

    def get_left_status(self,ct):
        '''Get the current position of the left step, '''
        nstep = np.floor(ct/self.dT)   # the step
        tstep = ct - nstep*self.dT  # the progress time within the step

        self.cleftposbase = self.ipos + self.lf.vel * nstep*1.0/self.fs
        self.cleftpos = self.cleftposbase + self.stepvel*tstep

        self.lf.cal_force(tstep)
        self.cleftforce = self.lf.force

    def get_right_status(self,ct):
        '''Get the current position of the left step, '''

        nstep = np.floor((ct-1.0/self.fs)/self.dT)   # the step
        tstep = (ct-1.0/self.fs) - nstep*self.dT  # the progress time within the step

        self.crightposbase = self.ipos + 1.0/self.fs*self.stepvel+ self.rf.vel * (nstep)*1.0/self.fs
        self.crightpos = self.crightposbase + self.stepvel * tstep

        self.rf.cal_force(tstep)
        if ct > 1.0 / self.fs:
            self.crightforce = self.rf.force
        else:
            self.crightforce = 0




if __name__ == '__main__':
    freq = 2
    initialpos = 0
    walk = WalkDoubleStep(initialpos, freq,0.2, 20, mode='heel')
    n = 1000
    ttime = 10
    tspace = np.linspace(0,ttime,n)
    lf = np.zeros((2,n))
    rf = np.zeros((2,n))

    ind = 0
    for ti in tspace:
        walk.update(ti)
        lf[0,ind] = walk.cleftpos
        lf[1,ind] = walk.cleftforce
        rf[0,ind] = walk.crightpos
        rf[1,ind] = walk.crightforce
        ind = ind + 1

    plt.plot(lf[0, :], lf[1, :], 'r*-')
    plt.plot(rf[0, :], rf[1, :], 'b*-')
    plt.show()



    print 1

"""
    freq = 2
    f1 = FootFall(freq)

    tspace = np.linspace(0,f1.timestep,1000)  # get the time split for single step
    flist = np.zeros(np.size(tspace))
    ind = 0
    for ti in tspace:
        flist[ind] = f1.cal_force(ti)
        ind = ind + 1
    plt.plot(tspace, flist, 'r-')
    plt.show()

"""