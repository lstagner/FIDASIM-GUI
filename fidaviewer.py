#!/usr/bin/env python

import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

import numpy as np

from Tkinter import *
from tkFileDialog import askdirectory
import ttk

class viewer:

    def __init__(self,parent):

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

        #Spectra Frame
        ttk.Button(self.spectra_frame,text='Plot Spectra',command=(lambda: self.update_plot1())).pack(expand=Y,fill=BOTH)
        #NPA Frame
        ttk.Button(self.npa_frame,text='Plot NPA',command=(lambda: self.update_plot2())).pack(expand=Y,fill=BOTH)
        #Neutrals Frame
        ttk.Button(self.neutrals_frame,text='Plot Neutrals').pack(expand=Y,fill=BOTH)
        #Weights Frame
        ttk.Button(self.weights_frame,text='Plot Weights').pack(expand=Y,fill=BOTH)
  
        self.fig=plt.Figure(figsize=(5,4),dpi=100)
        self.ax=self.fig.add_subplot(111)
        self.ax.plot(np.arange(10),np.arange(10))
        self.canvas=FigureCanvasTkAgg(self.fig,master=parent)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=RIGHT,expand=Y,fill=BOTH)
        self.toolbar=NavigationToolbar2TkAgg(self.canvas,parent)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=TOP,fill=BOTH,expand=Y)      
        self.load_dir()

    def load_dir(self):
        self.dir=askdirectory()
 
    def update_plot1(self):
        self.ax.cla()
        self.ax.plot(np.arange(10),np.arange(10)**2)
        self.canvas.show()
 
    def update_plot2(self):
        self.ax.cla()
        x=np.linspace(0,2*np.pi)
        self.ax.plot(x,np.sin(x))
        self.canvas.show()

if __name__=='__main__':
    root=Tk()
    viewer(root)
    root.mainloop()

