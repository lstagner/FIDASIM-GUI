#!/usr/bin/env python

import sys
import os
import glob
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#from scipy.io import netcdf
import numpy as np
import tkinter as tk
from tkinter.filedialog import askdirectory
#from tkinter import *
from tkinter import ttk
import h5py

"""
Todo
----
* use f90nml python package to parse the fortran namelist file and find what was and wasn't calculated
"""

#def read_ncdf(file,vars=[]):
#    """ Reads a netCDF 3 file and returns a dict with its variables
#
#    Parameters
#    ----------
#    file : string
#        The netCDF file to be read
#    vars : string list, optional
#        List of variables to be read from file
#
#    Returns
#    -------
#    d : A python dict containing the files variables
#
#    Examples
#    --------
#    >>> d=read_ncdf('test.cdf')
#    >>> d1=read_ncdf('test.cdf',vars=['var1'])
#    >>> d2=read_ncdf('test.cdf',vars=['var1','var2'])
#
#    """
#    try:
#        f=netcdf.netcdf_file(file,'r',mmap=False)
#    except IOError:
#        print('Error: Cannot open file'+file)
#        return 0
#
#    if vars == []: vars=f.variables.keys()
#    d=dict((v,f.variables[v].data) for v in vars)
#    f.close()
#    return d

def load_dict_from_hdf5(filename):
    """
    ....
    """
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')


def recursively_load_dict_contents_from_group(h5file, path):
    """
    ....
    """
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item.value
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans


class Spectra:
    """ Spectra object that contains plot methods and parameters"""
    def __init__(self,dir):
        self._has_spectra = True if glob.glob(dir+'*_spectra.h5') else False
        print(self._has_spectra)
        if self._has_spectra:
            spec = load_dict_from_hdf5(glob.glob(dir+'*_spectra.h5')[0])
            self.lam = spec['lambda']
            if 'brems' in spec:
                self.brems = spec['brems']
            if 'fida' in spec:
                self.fida = spec['fida']
            if 'full' in spec:
                self.full = spec['full']
            if
            self.half = spec['half']
            self.third = spec['third']
            self.halo = spec['halo']

            self.channels=dict(('Channel '+str(i+1),i) for i in range(0,self.brems.shape[0]))
            self.wl_min = np.min(self.lam)
            self.wl_max = np.max(self.lam)
            self.dlam=np.abs(self.lam[1]-self.lam[0])

            self.chan = tk.StringVar(value='Channel 1')
            self.nbi_on = tk.BooleanVar(value=False)
            self.fida_on = tk.BooleanVar(value=False)
            self.brems_on = tk.BooleanVar(value=True)
            self.legend_on = tk.BooleanVar(value=True)

            if sum(self.full[0,:]) != 0:
                self.nbi_on.set(True)

            if sum(self.fida[0,:]) != 0:
                self.fida_on.set(True)

    def plot_spectra(self,fig,canvas):
        if self._has_spectra:
            ch=self.channels[self.chan.get()]
            lam=self.lam
            if self.brems_on.get():
                brems=self.brems[ch,:]
            else:
                brems=0.0

            fig.clf()
            ax=fig.add_subplot(111)
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


            if self.fida_on.get() or self.nbi_on.get():
                if self.legend_on.get(): ax.legend()
                ax.set_yscale('log')
                ax.set_xlabel('Wavelength [nm]')
                ax.set_ylabel('Ph/(s*nm*sr*m^2)')
                ax.set_title(self.chan.get())
                canvas.show()
            else: print('SPECTRA: No Spectra Selected')

        else: print('SPECTRA: No file')

    def plot_intensity(self,fig,canvas):
        if self._has_spectra:
            w1 = self.lam >= self.wl_min
            w2 = self.lam <= self.wl_max
            w=np.logical_and(w1,w2)
            intens=np.sum(self.fida[:,w],axis=1)*self.dlam
            ch=range(1,len(intens)+1)
            fig.clf()
            ax=fig.add_subplot(111)
            ax.plot(ch,intens)
            ax.set_title('FIDA Intensity vs. Channel')
            ax.set_ylabel('Ph/(s*sr*m^2)')
            ax.set_xlabel('Channel Number')
            ax.set_yscale('log')
            canvas.show()
        else: print('SPECTRA: No file')

class NPA:
    """ NPA object that contains plot methods and parameters"""
    def __init__(self,dir):
        self._has_npa= True if glob.glob(dir+'*_npa.cdf') else False
        self._has_wght= True if glob.glob(dir+'*_npa_weights.cdf') else False
        self._has_neut= True if glob.glob(dir+'*_neutrals.cdf') else False
        self._has_geo= True if glob.glob(dir+'*_inputs.cdf') else False

        if self._has_npa:
            npa = load_dict_from_hdf5(glob.glob(dir+'*_npa.cdf')[0])
            self.npa_energy=npa['energy']
            self.npa_flux=npa['flux']
            self.ipos=npa['ipos']
            self.counts=npa['counts']
        if self._has_wght:
            wght = load_dict_from_hdf5(glob.glob(dir+'*_npa_weights.cdf')[0])
            self.w_energy=wght['energy']
            self.w_flux=wght['flux']
        if self._has_neut:
            neut = load_dict_from_hdf5(glob.glob(dir+'*_neutrals.cdf')[0])
            self.dens=neut['fdens'].sum(0).sum(0)+neut['hdens'].sum(0).sum(0)+\
                      neut['tdens'].sum(0).sum(0)+neut['halodens'].sum(0).sum(0)
        if self._has_geo:
            geo = load_dict_from_hdf5(glob.glob(dir+'*_inputs.cdf')[0],vars=['x_grid','y_grid','xlos','ylos','xlens','ylens','chan_id'])
            self.x_grid=geo['x_grid']
            self.y_grid=geo['y_grid']
            chan_id=geo['chan_id']
            w = chan_id == 1
            self.xlos=geo['xlos'][w]
            self.ylos=geo['ylos'][w]
            self.xlens=geo['xlens'][w]
            self.ylens=geo['ylens'][w]

        if (self._has_npa or self._has_wght):
            self.channels=dict(('Channel '+str(i+1),i) for i in range(0,3))

        self.chan=tk.StringVar(value='Channel 1')

    def plot_neutral_birth(self,fig,canvas):
        if self._has_npa:
            fig.clf()
            ax=fig.add_subplot(111)
            ch=self.channels[self.chan.get()]
            if self._has_neut and self._has_geo:
                ax.plot(self.x_grid[0,:,:],self.y_grid[0,:,:],'k,')
                ax.contour(self.x_grid[0,:,:],self.y_grid[0,:,:],self.dens,20)
                ax.plot([self.xlos[ch],self.xlens[ch]],[self.ylos[ch],self.ylens[ch]],'k')
            ax.plot(self.ipos[ch,0,0:self.counts[ch]],self.ipos[ch,1,0:self.counts[ch]],'k,',alpha=.3)
            ax.set_title('Neutral Birth Position')
            ax.set_xlim(min(self.x_grid[0,0,:]),max(self.x_grid[0,0,:]))
            ax.set_ylim(min(self.y_grid[0,:,0]),max(self.y_grid[0,:,0]))
            ax.set_xlabel('x [cm]')
            ax.set_ylabel('y [cm]')
            canvas.show()
        else: print('NPA: No file')

    def plot_flux(self,fig,canvas):
        if self._has_npa or self._has_wght:
            fig.clf()
            ax=fig.add_subplot(111)
            ch=self.channels[self.chan.get()]
            if self._has_npa:
                ax.step(self.npa_energy,self.npa_flux[ch,:],label='MC Flux')
            if self._has_wght:
                ax.plot(self.w_energy,self.w_flux[ch,:],label='WF Flux')

            ax.legend()
            ax.set_title('Neutral Flux: '+self.chan.get())
            ax.set_ylabel('Flux')
            ax.set_xlabel('Energy [keV]')
            canvas.show()
        else: print('NPA: No file')

class Weights:
    """ Weights object that contains plot methods and parameters"""
    def __init__(self,dir):
        self._has_npa_wght= True if glob.glob(dir+'*_npa_weights.cdf') else False
        self._has_fida_wght= True if glob.glob(dir+'*_fida_weights.cdf') else False

        if self._has_fida_wght:
            fida=load_dict_from_hdf5(glob.glob(dir+'*_fida_weights.cdf')[0])
            self.f_energy=fida['energy']
            self.f_pitch=fida['pitch']
            self.lam=fida['lambda']
            self.dlam=np.abs(self.lam[1]-self.lam[0])
            self.wl_max=np.max(self.lam)
            self.wl_min=np.min(self.lam)
            self.f_rad=fida['radius']
            self.f_wght=fida['wfunct']
            self.f_chan=len(self.f_rad)
            self.fida_chans=dict(('Channel '+str(i+1),i) for i in range(0,self.f_chan))

        if self._has_npa_wght:
            npa=load_dict_from_hdf5(glob.glob(dir+'*_npa_weights.cdf')[0])
            self.n_energy=npa['energy']
            self.n_pitch=npa['pitch']
            self.n_wght=npa['wfunct']
            self.n_rad=npa['radius']
            self.n_chan=len(self.n_rad)
            self.npa_chans=dict(('Channel '+str(i+1),i) for i in range(0,self.n_chan))

        self.lam_val=tk.DoubleVar(value=655.0)
        self.fida_chan=tk.StringVar(value='Channel 1')
        self.npa_chan=tk.StringVar(value='Channel 1')

    def plot_npa_weights(self,fig,canvas):
        if self._has_npa_wght:
            ch=self.npa_chans[self.npa_chan.get()]
            fig.clf()
            ax=fig.add_subplot(111)
            c=ax.contourf(self.n_energy,self.n_pitch,self.n_wght[ch,:,:],50)
            fig.colorbar(c)
            ax.set_title('NPA Weight')
            ax.set_ylabel('Pitch')
            ax.set_xlabel('Energy [keV]')
            canvas.show()

    def plot_fida_weights(self,fig,canvas):
        if self._has_fida_wght:
            ch=self.fida_chans[self.fida_chan.get()]
            wl=float(self.lam_val.get())
            ind=np.argmin(np.abs(self.lam-wl))
            fig.clf()
            ax=fig.add_subplot(111)
            c=ax.contourf(self.f_energy,self.f_pitch,self.f_wght[ch,:,:,ind],30)
            fig.colorbar(c)
            ax.set_xlabel('Energy [keV]')
            ax.set_ylabel('Pitch')
            ax.set_title('FIDA Weight')
            canvas.show()

class Neutrals:
    """ Neutrals object that contains plot methods and parameters"""
    def __init__(self,dir):
        self._has_neut= True if glob.glob(dir+'*_neutrals.h5') else False
        self._has_geo= True if glob.glob(dir+'*_inputs.h5') else False

        if self._has_neut and self._has_geo:
            neut=load_dict_from_hdf5(glob.glob(dir+'*_neutrals.h5')[0])
            geo=load_dict_from_hdf5(glob.glob(dir+'*_inputs.h5')[0],vars=['x_grid','y_grid','z_grid','u_grid','v_grid','w_grid'])
            self.fdens=neut['fdens'].sum(0)
            self.hdens=neut['hdens'].sum(0)
            self.tdens=neut['tdens'].sum(0)
            self.halodens=neut['halodens'].sum(0)
            self.x_grid=geo['x_grid']
            self.y_grid=geo['y_grid']
            self.z_grid=geo['z_grid']
            self.u_grid=geo['u_grid']
            self.v_grid=geo['v_grid']
            self.w_grid=geo['w_grid']

        ##Radio Buttons Variable
        self.plot_type=tk.StringVar(value='XY')

        ##Checkbox Variables
        self.use_uvw = tk.BooleanVar(value=False)
        self.full_on = tk.BooleanVar(value=True)
        self.half_on = tk.BooleanVar(value=True)
        self.third_on = tk.BooleanVar(value=True)
        self.halo_on = tk.BooleanVar(value=True)

    def plot_neutrals(self,fig,canvas):
        full_on = self.full_on.get()
        half_on = self.half_on.get()
        third_on = self.third_on.get()
        halo_on = self.halo_on.get()
        torf=lambda T: 1 if T else 0
        if (self._has_neut and self._has_geo) and (full_on or half_on or third_on or halo_on):
            fig.clf()
            ax=fig.add_subplot(111)
            pt=self.plot_type.get()
            if pt == 'X':
                if self.use_uvw.get():
                    x=self.u_grid[0,0,:]
                    ax.set_xlabel('U [cm]')
                else:
                    x=self.x_grid[0,0,:]
                    ax.set_xlabel('X [cm]')
                fdens=self.fdens.sum(0).sum(0)
                hdens=self.hdens.sum(0).sum(0)
                tdens=self.tdens.sum(0).sum(0)
                halodens=self.halodens.sum(0).sum(0)
                if full_on: ax.plot(x,fdens,label='Full')
                if half_on: ax.plot(x,hdens,label='Half')
                if third_on: ax.plot(x,tdens,label='Third')
                if halo_on: ax.plot(x,halodens,label='Halo')
                ax.legend()
                ax.set_title('Neutral Density')
                ax.set_ylabel('Density [cm^-3]')
                canvas.show()
            if pt == 'Y':
                if self.use_uvw.get():
                    x=self.v_grid[0,:,0]
                    ax.set_xlabel('V [cm]')
                else:
                    x=self.y_grid[0,:,0]
                    ax.set_xlabel('Y [cm]')
                fdens=self.fdens.sum(0).sum(1)
                hdens=self.hdens.sum(0).sum(1)
                tdens=self.tdens.sum(0).sum(1)
                halodens=self.halodens.sum(0).sum(1)
                if full_on: ax.plot(x,fdens,label='Full')
                if half_on: ax.plot(x,hdens,label='Half')
                if third_on: ax.plot(x,tdens,label='Third')
                if halo_on: ax.plot(x,halodens,label='Halo')
                ax.legend()
                ax.set_title('Neutral Density')
                ax.set_ylabel('Density [cm^-3]')
                canvas.show()
            if pt == 'Z':
                if self.use_uvw.get():
                    x=self.w_grid[:,0,0]
                    ax.set_xlabel('W [cm]')
                else:
                    x=self.z_grid[:,0,0]
                    ax.set_xlabel('Z [cm]')
                fdens=self.fdens.sum(1).sum(1)
                hdens=self.hdens.sum(1).sum(1)
                tdens=self.tdens.sum(1).sum(1)
                halodens=self.halodens.sum(1).sum(1)
                if full_on: ax.plot(x,fdens,label='Full')
                if half_on: ax.plot(x,hdens,label='Half')
                if third_on: ax.plot(x,tdens,label='Third')
                if halo_on: ax.plot(x,halodens,label='Halo')
                ax.legend()
                ax.set_title('Neutral Density')
                ax.set_ylabel('Density [cm^-3]')
                canvas.show()
            if pt == 'XY':
                if self.use_uvw.get():
                    x=self.u_grid[0,:,:]
                    y=self.v_grid[0,:,:]
                    ax.set_xlabel('U [cm]')
                    ax.set_ylabel('V [cm]')
                else:
                    x=self.x_grid[0,:,:]
                    y=self.y_grid[0,:,:]
                    ax.set_xlabel('X [cm]')
                    ax.set_ylabel('Y [cm]')

                fdens=self.fdens.sum(0)
                hdens=self.hdens.sum(0)
                tdens=self.tdens.sum(0)
                halodens=self.halodens.sum(0)

                dens=fdens*torf(full_on)+hdens*torf(half_on)+tdens*torf(third_on)+halodens*torf(halo_on)
                c=ax.contourf(x,y,dens,50)
                fig.colorbar(c)
                ax.set_title('Neutral Density')
                canvas.show()
            if pt == 'XZ':
                if self.use_uvw.get():
                    x=self.u_grid[:,0,:]
                    y=self.w_grid[:,0,:]
                    ax.set_xlabel('U [cm]')
                    ax.set_ylabel('W [cm]')
                else:
                    x=self.x_grid[:,0,:]
                    y=self.z_grid[:,0,:]
                    ax.set_xlabel('X [cm]')
                    ax.set_ylabel('Z [cm]')

                fdens=self.fdens.sum(1)
                hdens=self.hdens.sum(1)
                tdens=self.tdens.sum(1)
                halodens=self.halodens.sum(1)

                dens=fdens*torf(full_on)+hdens*torf(half_on)+tdens*torf(third_on)+halodens*torf(halo_on)
                c=ax.contourf(x,y,dens,50)
                fig.colorbar(c)
                ax.set_title('Neutral Density')
                canvas.show()
            if pt == 'YZ':
                if self.use_uvw.get():
                    x=self.v_grid[:,:,0]
                    y=self.w_grid[:,:,0]
                    ax.set_xlabel('V [cm]')
                    ax.set_ylabel('W [cm]')
                else:
                    x=self.y_grid[:,:,0]
                    y=self.z_grid[:,:,0]
                    ax.set_xlabel('Y [cm]')
                    ax.set_ylabel('Z [cm]')

                fdens=self.fdens.sum(2)
                hdens=self.hdens.sum(2)
                tdens=self.tdens.sum(2)
                halodens=self.halodens.sum(2)

                dens=fdens*torf(full_on)+hdens*torf(half_on)+tdens*torf(third_on)+halodens*torf(halo_on)
                c=ax.contourf(x,y,dens,50)
                fig.colorbar(c)
                ax.set_title('Neutral Density')
                canvas.show()

class Viewer:
    """Class that contains FIDAsim result viewer window"""
    def __init__(self,parent):

        self.load_dir()
        parent.title('FIDAviewer')

        #Make MenuBar
        self.MenuBar=tk.Menu(parent)
        parent.config(menu=self.MenuBar)
        self.file=tk.Menu(self.MenuBar,tearoff=False)
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
        self.nb.add(self.weights_frame,text='Weights')
        self.nb.pack(side=tk.LEFT,expand=tk.Y,fill=tk.BOTH)

        self.fig=plt.Figure(figsize=(6,5),dpi=100)
        self.ax=self.fig.add_subplot(111)
        self.canvas=FigureCanvasTkAgg(self.fig,master=parent)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT)
        self.toolbar=NavigationToolbar2TkAgg(self.canvas,parent)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP,expand=tk.Y,fill=tk.BOTH)

        #Spectra Frame
        if self.spec._has_spectra:
	        ttk.Combobox(self.spectra_frame,textvariable=self.spec.chan,\
    	        values=tuple(self.spec.channels.keys())).pack()
        	ttk.Checkbutton(self.spectra_frame,text='Hide NBI', variable=self.spec.nbi_on,\
            	onvalue=False,offvalue=True).pack()
        	ttk.Checkbutton(self.spectra_frame,text='Hide FIDA', variable=self.spec.fida_on,\
            	onvalue=False,offvalue=True).pack()
        	ttk.Checkbutton(self.spectra_frame,text='No Bremsstrahlung', variable=self.spec.brems_on,\
            	onvalue=False,offvalue=True).pack()
        	ttk.Checkbutton(self.spectra_frame,text='Hide Legend', variable=self.spec.legend_on,\
            	onvalue=False,offvalue=True).pack()

        	ttk.Button(self.spectra_frame,text='Plot Spectra',\
            	command=(lambda: self.spec.plot_spectra(self.fig,self.canvas))).pack(side=tk.TOP,expand=tk.Y,fill=tk.BOTH)
        	ttk.Button(self.spectra_frame,text='Plot Intensity',\
            	command=(lambda: self.spec.plot_intensity(self.fig,self.canvas))).pack(side=tk.TOP,expand=tk.Y,fill=tk.BOTH)

        #NPA Frame
        if self.npa._has_npa:
        	ttk.Combobox(self.npa_frame,textvariable=self.npa.chan,\
            	values=tuple(self.npa.channels.keys())).pack()
        	ttk.Button(self.npa_frame,text='Plot Neutral Birth',\
            	command=(lambda: self.npa.plot_neutral_birth(self.fig,self.canvas))).pack(side=tk.TOP,expand=tk.Y,fill=tk.BOTH)
        	ttk.Button(self.npa_frame,text='Plot Flux',\
            	command=(lambda: self.npa.plot_flux(self.fig,self.canvas))).pack(side=tk.TOP,expand=tk.Y,fill=tk.BOTH)

        #Neutrals Frame
        ttk.Radiobutton(self.neutrals_frame,text='Density vs X',variable=self.neut.plot_type,value='X').pack()
        ttk.Radiobutton(self.neutrals_frame,text='Density vs Y',variable=self.neut.plot_type,value='Y').pack()
        ttk.Radiobutton(self.neutrals_frame,text='Density vs Z',variable=self.neut.plot_type,value='Z').pack()
        ttk.Radiobutton(self.neutrals_frame,text='Contour XY',variable=self.neut.plot_type,value='XY').pack()
        ttk.Radiobutton(self.neutrals_frame,text='Contour XZ',variable=self.neut.plot_type,value='XZ').pack()
        ttk.Radiobutton(self.neutrals_frame,text='Contour YZ',variable=self.neut.plot_type,value='YZ').pack()


        ttk.Checkbutton(self.neutrals_frame,text='Use Machine Coordinates', variable=self.neut.use_uvw,\
            onvalue=True,offvalue=False).pack()
        ttk.Checkbutton(self.neutrals_frame,text='Hide Full', variable=self.neut.full_on,\
            onvalue=False,offvalue=True).pack()
        ttk.Checkbutton(self.neutrals_frame,text='Hide Half', variable=self.neut.half_on,\
            onvalue=False,offvalue=True).pack()
        ttk.Checkbutton(self.neutrals_frame,text='Hide Third', variable=self.neut.third_on,\
            onvalue=False,offvalue=True).pack()
        ttk.Checkbutton(self.neutrals_frame,text='Hide Halo', variable=self.neut.halo_on,\
            onvalue=False,offvalue=True).pack()

        ttk.Button(self.neutrals_frame,text='Plot',\
            command=(lambda: self.neut.plot_neutrals(self.fig,self.canvas))).pack(expand=tk.Y,fill=tk.BOTH)

        #Weights Frame
        if self.wght._has_fida_wght:
        	ttk.Combobox(self.weights_frame,textvariable=self.wght.fida_chan,\
            	values=tuple(self.wght.fida_chans.keys())).pack()

        	tk.Scale(self.weights_frame,orient=tk.HORIZONTAL,length=200,\
            	from_=self.wght.wl_min,to=self.wght.wl_max,resolution=self.wght.dlam,variable=self.wght.lam_val).pack()

        	ttk.Button(self.weights_frame,text='Plot FIDA Weights',\
            	command=(lambda: self.wght.plot_fida_weights(self.fig,self.canvas))).pack(side=tk.TOP,expand=tk.Y,fill=tk.BOTH)

        if self.wght._has_npa_wght:
        	ttk.Combobox(self.weights_frame,textvariable=self.wght.npa_chan,\
            	values=tuple(self.wght.npa_chans.keys())).pack()

        	ttk.Button(self.weights_frame,text='Plot NPA Weights',\
            	command=(lambda: self.wght.plot_npa_weights(self.fig,self.canvas))).pack(side=tk.TOP,expand=tk.Y,fill=tk.BOTH)


    def load_dir(self):
        self.dir = askdirectory()
        if self.dir == '': self.dir=os.getcwd()
        if self.dir[-1] != '/': self.dir+='/'

        self.spec = Spectra(self.dir)
        self.npa = NPA(self.dir)
        self.neut = Neutrals(self.dir)
        self.wght = Weights(self.dir)

if __name__=='__main__':
    root=tk.Tk()
    Viewer(root)
    root.mainloop()

