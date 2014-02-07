#!/usr/bin/env python

import sys
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from scipy.io import netcdf

import numpy as np

from Tkinter import *
from tkFileDialog import askdirectory
import ttk

def read_ncdf(file,vars=[]):
    try:
        f=netcdf.netcdf_file(file,'r')
    except IOError:
        print('Error: Cannot open file'+file)
        return 0

    if vars == []: vars=f.variables.keys()
    d=dict((v,f.variables[v].data) for v in vars)
    f.close()
    return d

class spectra:

    def __init__(self,dir):
        runid=os.path.basename(os.path.normpath(dir))
        self.has_spectra=os.path.isfile(dir+runid+'_spectra.cdf')

        if self.has_spectra:
            spec = read_ncdf(dir+runid+'_spectra.cdf')
            self.lam = spec['lambda']
            self.brems = spec['brems']
            self.fida = spec['fida']
            self.full = spec['full']
            self.half = spec['half']
            self.third = spec['third']
            self.halo = spec['halo']

            self.channels=dict(('Channel '+str(i+1),i) for i in range(0,self.brems.shape[0]))
            self.wl_min = np.min(self.lam)
            self.wl_max = np.max(self.lam)
            self.dlam=np.abs(self.lam[1]-self.lam[0])

            self.chan = StringVar(value='Channel 1')
            self.nbi_on = BooleanVar(value=False)
            self.fida_on = BooleanVar(value=False)
            self.brems_on = BooleanVar(value=True)
            self.legend_on = BooleanVar(value=True)

            if sum(self.full[0,:]) != 0: 
                self.nbi_on.set(True)

            if sum(self.fida[0,:]) != 0:
                self.fida_on.set(True)

    def plot_spectra(self,ax,canvas):
        if self.has_spectra:
            ch=self.channels[self.chan.get()]
            lam=self.lam
            if self.brems_on.get():
                brems=self.brems[ch,:]
            else:
                brems=0.0

            ax.cla()
            if self.nbi_on.get():
               full=self.full[ch,:]+brems
               half=self.half[ch,:]+brems
               third=self.third[ch,:]+brems
               halo=self.halo[ch,:]+brems

               ax.plot(lam,full,label='Full')
               ax.plot(lam,half,label='Half')
               ax.plot(lam,third,label='Third')
               ax.plot(lam,halo,label='Halo')

            if self.fida_on.get():
               fida=self.fida[ch,:]+brems
               ax.plot(lam,fida,label='Fida')

            if self.legend_on.get(): ax.legend()

            ax.set_yscale('log')
            ax.set_xlabel('Wavelength [nm]')
            ax.set_ylabel('Ph/(s*nm*sr*m^2)')
            ax.set_title(self.chan.get())
            canvas.show()
        else: print('SPECTRA: No file to plot')

    def plot_intensity(self,ax,canvas):
        if self.has_spectra:
            w1 = self.lam >= self.wl_min
            w2 = self.lam <= self.wl_max
            w=np.logical_and(w1,w2)
            intens=np.sum(self.fida[:,w],axis=1)*self.dlam
            ch=range(1,len(intens)+1)
            ax.cla()
            ax.plot(ch,intens)
            ax.set_title('Intensity vs. Channel')
            ax.set_ylabel('Ph/(s*sr*m^2)')
            ax.set_xlabel('Channel Number')
            ax.set_yscale('log')
            canvas.show()

class npa:

    def __init__(self,dir):
        runid=os.path.basename(os.path.normpath(dir))
        self.has_npa=os.path.isfile(dir+runid+'_npa.cdf')
        self.has_wght=os.path.isfile(dir+runid+'_npa_weights.cdf')

        if self.has_npa:
            npa = read_ncdf(dir+runid+'_npa.cdf')
            self.npa_energy=npa['energy']
            self.npa_flux=npa['flux']
        if self.has_wght:
            wght = read_ncdf(dir+runid+'_npa_weights.cdf')
            self.w_energy=wght['energy']
            self.w_flux=wght['flux']

class weights:

    def __init__(self,dir):
        self.x=0

class neutrals:

    def __init__(self,dir):
        self.x=0
        
class viewer:

    def __init__(self,parent):

        self.load_dir()
        parent.title('FIDAviewer')

        #Make MenuBar
        self.MenuBar=Menu(parent)
        parent.config(menu=self.MenuBar)
        self.file=Menu(self.MenuBar,tearoff=False)
        self.file.add_command(label='Load Run',command=(lambda: self.load_dir()))
        self.file.add_command(label='Quit',command=(lambda: sys.exit()))
        self.MenuBar.add_cascade(label='File',menu=self.file,underline=0)

        #Make Notebook
        self.nb=ttk.Notebook(parent)
        self.spectra_frame=ttk.Frame(self.nb)         
        self.npa_frame=ttk.Frame(self.nb)         
        self.neutrals_frame=ttk.Frame(self.nb)         
        self.weights_frame=ttk.Frame(self.nb)
        self.nb.add(self.spectra_frame,text='Spectra')         
        self.nb.add(self.npa_frame,text='NPA')         
        self.nb.add(self.neutrals_frame,text='Neutrals')         
        self.nb.add(self.weights_frame,text='FIDA Weights')
        self.nb.pack(side=LEFT,expand=Y,fill=BOTH)

        self.fig=plt.Figure(figsize=(6,5),dpi=100)
        self.ax=self.fig.add_subplot(111)
        self.canvas=FigureCanvasTkAgg(self.fig,master=parent)
        self.canvas.get_tk_widget().pack(side=RIGHT)
        self.toolbar=NavigationToolbar2TkAgg(self.canvas,parent)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=TOP,expand=Y,fill=BOTH)      

        #Spectra Frame
        spec_chan=ttk.Combobox(self.spectra_frame,textvariable=self.spec.chan,\
            values=tuple(self.spec.channels.keys())).pack()
        ttk.Button(self.spectra_frame,text='Plot Spectra',\
            command=(lambda: self.spec.plot_spectra(self.ax,self.canvas))).pack(side=TOP,expand=Y,fill=BOTH)
        ttk.Button(self.spectra_frame,text='Plot Intensity',\
            command=(lambda: self.spec.plot_intensity(self.ax,self.canvas))).pack(side=TOP,expand=Y,fill=BOTH)

        #NPA Frame
        ttk.Button(self.npa_frame,text='Plot NPA').pack(expand=Y,fill=BOTH)
        #Neutrals Frame
        ttk.Button(self.neutrals_frame,text='Plot Neutrals').pack(expand=Y,fill=BOTH)
        #Weights Frame
        ttk.Button(self.weights_frame,text='Plot Weights').pack(expand=Y,fill=BOTH)
  

    def load_dir(self):
        self.dir = askdirectory()
        if self.dir == '': self.dir=os.getcwd()
        if self.dir[-1] != '/': self.dir+='/'

        self.runid = os.path.basename(os.path.normpath(self.dir))

        self.spec = spectra(self.dir)
        self.npa = npa(self.dir)

if __name__=='__main__':
    root=Tk()
    viewer(root)
    root.mainloop()

