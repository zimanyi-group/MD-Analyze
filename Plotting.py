import numpy as np
import sys
import time
import pdb
import math
import Atom_Manip as am

import scipy.optimize as opt

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure



import tkinter as tk


plt.style.use('seaborn-deep')
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['mediumblue', 'crimson','darkgreen', 'darkorange','crimson', 'darkorchid'])
plt.rcParams['figure.figsize'] = [12,8]
plt.rcParams['axes.linewidth'] = 1.7
plt.rcParams['lines.linewidth'] = 6.0
plt.rcParams['axes.grid'] = True
plt.rcParams['font.size'] = 22
plt.rcParams['font.family'] =  'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
#plt.rcParams['axes.facecolor'] = 'grey'
plt.ion()

m_colormap1 = np.array([[0.8,0.0,0.1], [0.1,0.0,0.8] ])


def clamp(val, min, max):
    if(val >= min and val<= max):
        return val
    return min if val < min else max


def Plot_Hist(dat):
    fig = plt.figure(figsize=[10,10])
    ax = fig.subplots(2,2, subplot_kw=dict(projection = '3d'))
    labels = ["Total", "Si","O", "x"]
    count = 0
    for i in range(2):
        ax[i,0].bar3d(dat[:,0], dat[:,1],0.0, 5.3,5.3, dat[:,2 + count], color = 'darkorchid')
        ax[i,0].set_title(labels[count])
        count+=1
        ax[i,1].bar3d(dat[:,0], dat[:,1],0.0, 5.3,5.3, dat[:,2 + count], color = 'crimson')
        ax[i,1].set_title(labels[count])
        count+=1
    plt.show(block = True)









class Dump_Plot():
   
    def __init__(self, sim):
        self.root = tk.Tk()
        self.fig = Figure(figsize=[10,8])
        self.ax = self.fig.subplots(1,2, subplot_kw=dict(projection = '3d'))
        self.labels = [ "Si","O", "x", "total"]
        self.can = FigureCanvasTkAgg(self.fig, master = self.root)
        self.quit = tk.Button(self.root, text="Quit", padx=50, pady=5, command= self.root.destroy)
        self.next = tk.Button(self.root, text="next",padx=50, pady=5, command=lambda: self.Update_Plot('next'))
        self.prev = tk.Button(self.root, text="prev",padx=50, pady=5, command=lambda: self.Update_Plot('prev'))
        self.TS_Label = tk.Label(self.root, bd=3, font = 'Times 18')
        self._sim = sim
        self.dat = am.Compute_Ratio(self._sim)
        self.p1 = None
        self.p2 = None
        self.mod = [0.2, 0.8, 1.0,1.0,1.0,0.5]
        self.clrmap1 = []
        self.clrmap2 = []

    def save(self, name):
        self.fig.savefig("fit_{0}.png".format(name))
        return

    def _clear(self):
        self.ax[0].clear()
        self.ax[1].clear()
        return

    def Set_Color1(self):
        zar = self.dat[:,2]
        self.clrmap1 = []
        delta = 18.0
        for z in zar:
            zp = clamp((z-np.min(zar))/delta, 0.0, 1.0)
            _r = 1.0 - zp
            _b = zp
            _color = np.array([_r, 0.0,_b])
            self.clrmap1.append(_color)


    def Set_Color2(self):
        zar = self.dat[:,3]
        self.clrmap2 = []
        delta = np.max(zar) - 8.0
        for z in zar:
            zp = clamp((z-8.0)/delta, 0.0, 1.0)
            _r = 1.0 - zp
            _g = zp
            _color = np.array([_r, 1.0,0.0])
            self.clrmap2.append(_color)
        return 

    def Update_Plot(self, dir):
        self._sim.Update(dir)
        self.dat = am.Compute_Ratio(self._sim)
        self.Set_Color1()
        self.Set_Color2()
        #pdb.set_trace()
        self._clear()
        self.ax[0].set_title(self.labels[0])
        self.ax[1].set_title(self.labels[1])
        self.ax[0].set_zlim(0,25)
        self.ax[1].set_zlim(0,25)
        self.ax[0].set_ylabel("Y")
        self.ax[1].set_ylabel("Y")
        self.ax[0].bar3d(self.dat[:,0], self.dat[:,1],0.0, 4.3,4.3, self.dat[:,2], color = self.clrmap1)
        self.ax[1].bar3d(self.dat[:,0], self.dat[:,1],0.0, 4.3,4.3, self.dat[:,3], color = self.clrmap2)
        self.TS_Label['text'] = "Timestep:\n{0:.0f}".format(self._sim.timestep)
        self.can.draw()
        return

    def plot_hist(self):
        self.TS_Label['text'] = "Timestep:\n{0:.0f}".format(self._sim.timestep)
        self.Set_Color1()
        self.Set_Color2()
        self.ax[0].set_zlim(0,25)
        self.ax[1].set_zlim(0,25)
        self.p1 = self.ax[0].bar3d(self.dat[:,0], self.dat[:,1],0.0, 4.3,4.3, self.dat[:,2], color = self.clrmap1)
        self.p2 = self.ax[1].bar3d(self.dat[:,0], self.dat[:,1],0.0, 4.3,4.3, self.dat[:,3], color = self.clrmap2)
        self.can.draw()
        self.prev.grid(row = 0, column = 1)
        self.next.grid(row = 0, column = 2)
        self.quit.grid(row = 1, column = 2)
        self.TS_Label.grid(row = 0, column = 0)
        self.can.get_tk_widget().grid(row = 1, column = 0)
        self.root.mainloop()
        return






class Dump_Plot_Map():
   
    def __init__(self, sim):
        self.root = tk.Tk()
        self.fig = Figure(figsize=[10,8])
        self.ax = self.fig.subplots(1,1)
        self.labels = [ "Si","O", "x", "total"]
        self.can = FigureCanvasTkAgg(self.fig, master = self.root)
        self.quit = tk.Button(self.root, text="Quit", padx=50, pady=5, command= self.root.destroy)
        self.next = tk.Button(self.root, text="next",padx=50, pady=5, command=lambda: self.Update_Plot('next'))
        self.prev = tk.Button(self.root, text="prev",padx=50, pady=5, command=lambda: self.Update_Plot('prev'))
        self.TS_Label = tk.Label(self.root, bd=3, font = 'Times 18')
        self._sim = sim
        self.dat = am.Compute_Ratio2(self._sim)
        self.p1 = None
        self.cb = None

    def save(self, name):
        self.fig.savefig("fit_{0}.png".format(name))
        return

    def _clear(self):
        self.ax.clear()
        #self.ax[1].clear()
        return



    def Update_Plot(self, dir):
        self._sim.Update(dir)
        self.dat = am.Compute_Ratio2(self._sim)
        si, ox, rat = self.dat
        #pdb.set_trace()
        self._clear()
        self.ax.set_title(self.labels[2])
        #self.ax[1].set_title(self.labels[2])
        self.p1 = self.ax.pcolormesh(rat,vmin = 0, vmax =2, cmap = 'viridis_r')
        self.TS_Label['text'] = "Timestep:\n{0:.0f}".format(self._sim.timestep)
        self.can.draw()
        return

    def Plot_Map(self):
        self.TS_Label['text'] = "Timestep:\n{0:.0f}".format(self._sim.timestep)
        si, ox, rat = self.dat
        self.p1 = self.ax.pcolormesh(rat, vmin = 0, vmax = 2, cmap = 'viridis_r')
        self.cb = self.fig.colorbar(self.p1, ax=self.ax)
        #self.p2 = self.ax[1].pcolormesh(rat, cmap='RdBu_r')
        #self.fig.colorbar(self.p2, ax=self.ax[1])
        self.can.draw()
        self.prev.grid(row = 0, column = 1)
        self.next.grid(row = 0, column = 2)
        self.quit.grid(row = 1, column = 2)
        self.TS_Label.grid(row = 0, column = 0)
        self.can.get_tk_widget().grid(row = 1, column = 0)
        self.root.mainloop()
        return




#top = tk.TopLevel()
#top.minsize(500,600)
#top.title('Scale Window')
#drop1 = tk.OptionMenu( top, self.opt , *self.labels )
#drop1.pack()
#drop1.wait_variable(opt)
#but = tk.Button(top, text="Compute Axis Scale", command = lambda: Compute_Scale(pts))
#but.pack()