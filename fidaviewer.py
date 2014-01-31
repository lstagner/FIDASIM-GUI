#!/usr/bin/env python

import sys
import fidaTools as ft
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
        self.nb.pack(expand=Y,fill=BOTH)

        #Spectra Frame
        ttk.Button(self.spectra_frame,text='Plot Spectra',command=(lambda: ft.plot_spectra(self.dir))).pack(expand=Y,fill=BOTH)
        #NPA Frame
        ttk.Button(self.npa_frame,text='Plot NPA',command=(lambda: ft.plot_npa(self.dir))).pack(expand=Y,fill=BOTH)
        #Neutrals Frame
        ttk.Button(self.neutrals_frame,text='Plot Neutrals',command=(lambda: ft.plot_neutrals(self.dir))).pack(expand=Y,fill=BOTH)
        #Weights Frame
        ttk.Button(self.weights_frame,text='Plot Weights',command=(lambda: ft.plot_fida_weights(self.dir))).pack(expand=Y,fill=BOTH)

        self.load_dir()

    def load_dir(self):
        self.dir=askdirectory()

if __name__=='__main__':
    root=Tk()
    viewer(root)
    root.mainloop()

