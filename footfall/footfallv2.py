import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

class FootFallPatch:
    """this defines footfallpatch"""
    def __init__(self, shape='rec',mode='uniform'):
        self.mode = mode
        self.shape = shape
        self.patcharea = None

    def setrec(self,rangeX,rangeY):
        self.recrx = rangeX
        self.recry = rangeY
        self.patcharea =[rangeX, rangeY]

    def getlengthdirect(self):
        """get the length in the direction"""

        if self.shape == 'rec':
            return self.recrx

    def getarea(self):
        return self.patcharea


class FootFall:
    def __init__(self, fs, patch):
        self.fs = fs      # left-to-right (single) frequency
        self.patch = patch
        self.duration = 5.8531 - 6.9773 * fs + 3.2242 * fs ** 2 - 0.515 * fs ** 3  # single footprint duration
        self.overlap = self.duration-1.0/self.fs
        self.lengthdirect = float(self.patch.getlengthdirect())
        self.dT = 2*self.duration - 2*self.overlap
        #
        self.vel_step =self.lengthdirect/self.duration                # moving velocity during the footfall duration
        self.vel_avg = 1.67 * fs ** 2 - 4.83 * fs + 4.5               # walking velocity
        self.vel_move = (self.vel_avg*1.0/self.fs -self.lengthdirect)/(self.duration-self.overlap*2)                                # foot move velocity

        self.cal_force_para()  # get force parameter based on fs

    def getpatcharea(self,basecoordx,basecoordy):
        if self.patch.shape == 'rec':
            width,height = self.patch.getarea()
            bottomleftcoordX = basecoordx
            bottomleftcoordY = basecoordy - height/2.0

        return {'base':(bottomleftcoordX,bottomleftcoordY),'width':width,'height':height}


    def cal_force(self, T_NOW):
        if T_NOW >0 and T_NOW<self.duration:

            self.force = 746 * ((self.k1 * T_NOW ** 1) + (self.k2 * T_NOW ** 2) +
                                (self.k3 * T_NOW ** 3) + (self.k4 * T_NOW ** 4) + (self.k5 * T_NOW ** 5) +
                                (self.k6 * T_NOW ** 6) + (self.k7 * T_NOW ** 7) + (self.k8 * T_NOW ** 8))
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


class Step:
    """this is the class for single step, either left or right"""
    def __init__(self,footfall,offsettime=0,offX=0,offY=0):
        self.footfall = footfall

        self.length = self.footfall.dT*self.footfall.vel_avg
        self.ot = offsettime
        self.odx = offX
        self.ody = offY
        self.nstep = 0
        self.tstep = 0
        self.cposbase = 0
        self.cpos = 0
        self.cforce = 0

    def update(self,ct):
        self.getstate(ct)
        self.getloc()
        self.getforce()

    def getstate(self,ct):
        if ct >= self.ot:
            self.nstep = np.floor((ct - self.ot) / self.footfall.dT)  # the step
            self.tstep = (ct - self.ot) - float(self.nstep) * self.footfall.dT  # the progress time within the step

    def getloc(self):
        """location of current time"""
        self.cposbase = self.odx + self.nstep * self.length
        self.patcharea = self.footfall.getpatcharea(self.cposbase,self.ody)


        if self.tstep <=self.footfall.duration:
            self.cpos = self.cposbase + self.footfall.vel_step * self.tstep
        else:
            self.cpos = self.cposbase + self.footfall.lengthdirect +\
                        self.footfall.vel_move * (self.tstep-self.footfall.duration)

    def getforce(self):
        """get current force"""
        self.footfall.cal_force(self.tstep)
        self.cforce =self.footfall.force

    def getvalues(self):
        return [self.nstep,self.tstep,self.cposbase,self.patcharea,self.cforce]

    def plt_patch(self):
        fig, ax = plt.subplots(1)
        rect = Rectangle(self.patcharea['base'],self.patcharea['width'],self.patcharea['height'])
        pc = PatchCollection([rect], facecolor='r', alpha=0.5, edgecolor=None)
        ax.add_collection(pc)

        plt.show()

class SingleWalk():
    """walk from single person"""
    def __init__(self,footfall):
        self.ls = Step(footfall,offsettime=0, offX=0,offY=0)
        self.rs = Step(footfall, offsettime=1.0/footfall.fs, offX=0.5,offY=0)
        self.lshist = []
        self.rshist = []
        self.cthist = []

    def update(self,ct):
        self.ls.update(ct)
        self.rs.update(ct)
        self.commit(ct)


        self. plt_time_patch()

    def commit(self,ct):
        self.lshist.append(self.ls.getvalues())
        self.rshist.append(self.rs.getvalues())
        self.cthist.append(ct)

    def postprocess(self):
        self.cthist = np.array(self.cthist)
        self.lshist = np.array(self.lshist)
        self.rshist = np.array(self.rshist)

    def plt_ls_time_force(self):
        fig = plt.figure()
        plt.plot(np.transpose(self.cthist), self.lshist[:,4], 'r*-')
        plt.show()

    def plt_rs_time_force(self):
        fig = plt.figure()
        plt.plot(np.transpose(self.cthist), self.rshist[:,4], 'r*-')
        plt.show()

    def plt_time_force(self):
        fig = plt.figure()
        plt.plot(np.transpose(self.cthist), self.lshist[:,4], 'r*-')
        plt.plot(np.transpose(self.cthist), self.rshist[:, 4], 'b*-')
        plt.show()

    def plt_time_loc(self):
        fig = plt.figure()
        plt.plot(np.transpose(self.cthist), self.lshist[:, 2], 'r*-')
        plt.plot(np.transpose(self.cthist), self.rshist[:, 2], 'b*-')
        plt.show()

    def plt_time_patch(self):
        self.ls.plt_patch()
        self.rs.plt_patch()

if __name__ == '__main__':

    # define pressure patch
    patch1 = FootFallPatch(shape='rec',mode='uniform')
    patch1.setrec(0.2,0.1)  # foot area 0.2m by 0.1 m

    # define footfall with specific frequency
    fs = 2
    FootFall1 = FootFall(fs,patch1)


    person1 = SingleWalk(FootFall1)
    n = 1000
    ttime = 10
    tspace = np.linspace(0,ttime,n)

    for ti in tspace:
        person1.update(ti)

    person1.postprocess()
    #person1.plt_time_force()
    #person1.plt_time_loc()
    print 1