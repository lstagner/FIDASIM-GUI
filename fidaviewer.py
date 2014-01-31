#!/usr/bin/env python

import sys
import fidaTools as ft
from Tkinter import *
import ttk

root=Tk()
root.title('FIDAviewer')

mainframe = ttk.Frame(root,padding="3 3 12 12")
mainframe.grid(column=0,row=0,sticky=(N,W,E,S))
mainframe.columnconfigure(0,weight=1)
mainframe.rowconfigure(0,weight=1)

ttk.Button(mainframe,text='Plot Spectra',command=(lambda: ft.plot_spectra('153071H02'))).grid(column=0,row=0,sticky=W)
ttk.Button(mainframe,text='Plot NPA',command=(lambda: ft.plot_npa('153071H02'))).grid(column=0,row=1,sticky=W)
ttk.Button(mainframe,text='Plot Neutrals',command=(lambda: ft.plot_neutrals('153071H02'))).grid(column=0,row=2,sticky=W)
ttk.Button(mainframe,text='Plot Weights',command=(lambda: ft.plot_fida_weights('153071H02'))).grid(column=0,row=3,sticky=W)
ttk.Button(mainframe,text='Exit',command=(lambda: sys.exit())).grid(column=0,row=4,sticky=W)

for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

root.mainloop()


