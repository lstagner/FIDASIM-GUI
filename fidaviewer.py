#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import numpy as np
import tkinter as tk
from tkinter.filedialog import askdirectory
from tkinter import ttk
import h5py
import f90nml
import collections


"""
Todo
----
* in taking mean of beam densities, should it only be for non-zero elements?
* optimize: can more stuff be loaded only when used? can more stuff be saved and not recalculated (ie set/get)?
* option to change volume element in neutral plotting for better fidelity in going from beam to mach coords
* find out if histogram2d give left edges or right
* rerun sample file sim setting all bools to true
* implement multiple matching filenames
* Make Spectra wl_min and wl_max changeable from gui
* get more intellegent h5 reader to just grab what's needed
* NPA needs work. I haven't used NPA data before - NGB
* currently seems to load neutrals (more?) twice. check this and fix
* add another tab to gui "Imaging" w/ "Lens" drop down. Choose spectra and wavelength range to integrate and make contour
* DONE - display msg like "No NPA file found" under each tab for clarity
* DONE - use f90nml python package to parse the fortran namelist file and find what was and wasn't calculated
* DONE - check for multiple matching filenames
* DONE - what is .get() business? This is tk setter/getter feature
* DONE - Make brems separate signal and stop adding to other spectra
* DONE - change xyz to uvw and vise versa in Neutrals
* DONE - neutral density legend units
* DONE - use histogram2d when ~np.array_equal(x, uniq(x_grid)), etc. ie beam coords != mach coords
* DONE - give window name of dir
* DONE - how to sort channels in drop down box?
* DONE - put NB name on plot (in geo file)

"""

def load_dict_from_hdf5(h5_filepath):
    """
    Load h5 file as a dict
    """
    def recursively_load_dict_contents_from_group(h5_obj, path):
        """
        Recursively load a dict from h5 file
        """
        ans = {}
        for key, item in h5_obj[path].items():
            if isinstance(item, h5py._hl.dataset.Dataset):
                ans[key] = item.value
            elif isinstance(item, h5py._hl.group.Group):
                ans[key] = recursively_load_dict_contents_from_group(h5_obj, path + key + '/')
        return ans

    with h5py.File(h5_filepath, 'r') as h5_obj:
        return recursively_load_dict_contents_from_group(h5_obj, '/')


class Spectra:
    """ Spectra object that contains plot methods and parameters"""
    def __init__(self, dir, nml):
        spec_files = glob.glob(dir + '*_spectra.h5')
        self._has_spectra = (len(spec_files) > 0)

        if self._has_spectra:
            print('Loading spectra')

            if len(spec_files) > 1:
                raise NotImplementedError('Multiple spectra files found')
            else:
                spec = load_dict_from_hdf5(spec_files[0])

            self.lam = spec['lambda']
            self.nchan = spec['nchan']
            self.channels = collections.OrderedDict(('Channel ' + str(i + 1), i) for i in range(self.nchan))
            self.wl_min = np.min(self.lam)
            self.wl_max = np.max(self.lam)
            self.dlam = np.abs(self.lam[1] - self.lam[0])
            self.chan = tk.StringVar(value = 'Channel 1')
            self.nbi_on = tk.BooleanVar(value = nml['calc_bes'] > 0)
            self.fida_on = tk.BooleanVar(value = nml['calc_fida'] > 0)
            self.brems_on = tk.BooleanVar(value = nml['calc_brems'] > 0)
            self.legend_on = tk.BooleanVar(value = True)
            if self.brems_on.get():
                self.brems = spec['brems']
            if self.fida_on.get():
                self.fida = spec['fida']
            if self.nbi_on.get():
                self.full = spec['full']
                self.half = spec['half']
                self.third = spec['third']
                self.halo = spec['halo']
        else:
            print('No spectra found')

    def plot_spectra(self, fig, canvas):
        if self._has_spectra:
            ch = self.channels[self.chan.get()]
            lam = self.lam

            fig.clf()
            ax = fig.add_subplot(111)

            if self.brems_on.get():
                brems = self.brems[ch, :]
                ax.plot(lam, brems, label = 'Brems')

            if self.nbi_on.get():
                full = self.full[ch, :]
                half = self.half[ch, :]
                third = self.third[ch, :]
                halo = self.halo[ch, :]

                ax.plot(lam, full, label = 'Full')
                ax.plot(lam, half, label = 'Half')
                ax.plot(lam, third, label = 'Third')
                ax.plot(lam, halo, label = 'Halo')

            if self.fida_on.get():
                fida = self.fida[ch, :]
                ax.plot(lam, fida, label = 'Fida')

            if self.brems_on.get() or self.fida_on.get() or self.nbi_on.get():
                if self.legend_on.get(): ax.legend()
                ax.set_yscale('log')
                ax.set_xlabel('Wavelength [nm]')
                ax.set_ylabel('$Ph\ /\ (s\ nm\ sr\ m^2)$')
                ax.set_title(self.chan.get())
                canvas.show()
            else:
                print('SPECTRA: No Spectra Selected')
        else:
            print('SPECTRA: No file')

    def plot_intensity(self, fig, canvas):
        if self._has_spectra:
            w1 = (self.lam >= self.wl_min)
            w2 = (self.lam <= self.wl_max)
            w = np.logical_and(w1, w2)
            intens = np.sum(self.fida[:, w], axis = 1) * self.dlam
            ch = range(1, len(intens) + 1)
            fig.clf()
            ax = fig.add_subplot(111)
            ax.plot(ch, intens)
            ax.set_title('FIDA Intensity vs. Channel')
            ax.set_ylabel('$Ph\ /\ (s\ sr\ m^2)$')
            ax.set_xlabel('Channel Number')
            ax.set_yscale('log')
            canvas.show()
        else: print('SPECTRA: No file')

class NPA:
    """ NPA object that contains plot methods and parameters"""
    def __init__(self,dir):
        npa_files = glob.glob(dir + '*_npa.h5')
        wght_files = glob.glob(dir + '*_npa_weights.h5')
        neut_files = glob.glob(dir + '*_neutrals.h5')

        self._has_npa = (len(npa_files) > 0)
        self._has_wght = (len(wght_files) > 0)
        self._has_neut = (len(neut_files) > 0)

        if self._has_npa:
            print('Loading NPA')

            if len(npa_files) > 1:
                raise NotImplementedError('Multiple NPA files found')
            else:
                npa = load_dict_from_hdf5(npa_files[0])

            self.npa_energy = npa['energy']
            self.npa_flux = npa['flux']
            self.nchan = npa['nchan']
        else:
            print('No NPA found')

        if self._has_wght:
            print('Loading NPA weights')

            if len(wght_files) > 1:
                raise NotImplementedError('Multiple NPA weight files found')
            else:
                wght = load_dict_from_hdf5(wght_files[0])

            self.w_energy = wght['energy']
            self.w_flux = wght['flux']
        else:
            print('No NPA weights found')

        if self._has_neut:

            if len(neut_files) > 1:
                raise NotImplementedError('Multiple neutrals files found')
            else:
                neut = load_dict_from_hdf5(neut_files[0])

            self.dens = neut['fdens'].sum(0).sum(0) + neut['hdens'].sum(0).sum(0) + \
                        neut['tdens'].sum(0).sum(0) + neut['halodens'].sum(0).sum(0)

#        if self._has_geo:
#            geo = load_dict_from_hdf5(glob.glob(dir + '*_geometry.h5')[0])  #,vars = ['x_grid','y_grid','xlos','ylos','xlens','ylens','chan_id'])
#            self.x_grid = geo['x_grid']
#            self.y_grid = geo['y_grid']
#            chan_id = geo['chan_id']
#            w = chan_id == 1
#            self.xlos = geo['xlos'][w]
#            self.ylos = geo['ylos'][w]
#            self.xlens = geo['xlens'][w]
#            self.ylens = geo['ylens'][w]

        if (self._has_npa or self._has_wght):
            self.channels = collections.OrderedDict(('Channel ' + str(i + 1), i) for i in range(0, self.nchan))  # should it be nchan not 3???

        self.chan = tk.StringVar(value = 'Channel 1')

    def plot_neutral_birth(self, fig, canvas):
        if self._has_npa:
            fig.clf()
            ax = fig.add_subplot(111)
            ch = self.channels[self.chan.get()]

            if self._has_neut:
                ax.plot(self.x_grid[0,:,:],self.y_grid[0,:,:],'k,')
                ax.contour(self.x_grid[0,:,:],self.y_grid[0,:,:],self.dens,20)
                ax.plot([self.xlos[ch],self.xlens[ch]],[self.ylos[ch],self.ylens[ch]],'k')

            ax.set_title('Neutral Birth Position')
            ax.set_xlim(min(self.x_grid[0,0,:]) ,max(self.x_grid[0,0,:]))
            ax.set_ylim(min(self.y_grid[0,:,0]),max(self.y_grid[0,:,0]))
            ax.set_xlabel('x [cm]')
            ax.set_ylabel('y [cm]')
            canvas.show()
        else:
            print('NPA: No file')

    def plot_flux(self,fig,canvas):
        if self._has_npa or self._has_wght:
            fig.clf()
            ax = fig.add_subplot(111)
            ch = self.channels[self.chan.get()]
            if self._has_npa:
                ax.step(self.npa_energy,self.npa_flux[ch,:],label = 'MC Flux')
            if self._has_wght:
                ax.plot(self.w_energy,self.w_flux[ch,:],label = 'WF Flux')

            ax.legend()
            ax.set_title('Neutral Flux: '+self.chan.get())
            ax.set_ylabel('Flux')
            ax.set_xlabel('Energy [keV]')
            canvas.show()
        else: print('NPA: No file')

class Weights:
    """ Weights object that contains plot methods and parameters"""
    def __init__(self,dir):
        npa_wght_files = glob.glob(dir + '*_npa_weights.h5')
        fida_wght_files = glob.glob(dir + '*_fida_weights.h5')

        self._has_npa_wght = (len(npa_wght_files) > 0)
        self._has_fida_wght = (len(fida_wght_files) > 0)

        if self._has_fida_wght:
            print('Loading FIDA weights')

            if len(npa_wght_files) > 1:
                raise NotImplementedError('Multiple FIDA weight files found')
            else:
                fida = load_dict_from_hdf5(fida_wght_files[0])
            self.f_energy = fida['energy']
            self.f_pitch = fida['pitch']
            self.lam = fida['lambda']
            self.dlam = np.abs(self.lam[1] - self.lam[0])
            self.wl_max = np.max(self.lam)
            self.wl_min = np.min(self.lam)
            self.f_rad = fida['radius']
            self.f_wght = fida['weight']
            self.f_chan = len(self.f_rad)
            self.fida_chans = dict(('Channel '+str(i+1),i) for i in range(0,self.f_chan))
        else:
            print('No FIDA weights found')

        if self._has_npa_wght:
            if len(npa_wght_files) > 1:
                raise NotImplementedError('Multiple NPA weight files found')
            else:
                npa = load_dict_from_hdf5(npa_wght_files[0])
            self.n_energy = npa['energy']
            self.n_pitch = npa['pitch']
            self.n_wght = npa['weight']
            self.n_rad = npa['radius']
            self.n_nchan = npa['nchan']  #len(self.n_rad)
            self.npa_chans = dict(('Channel ' + str(i + 1), i) for i in range(0, self.n_nchan))

        self.lam_val = tk.DoubleVar(value = 655.0)
        self.fida_chan = tk.StringVar(value = 'Channel 1')
        self.npa_chan = tk.StringVar(value = 'Channel 1')

    def plot_npa_weights(self,fig,canvas):
        if self._has_npa_wght:
            ch = self.npa_chans[self.npa_chan.get()]
            fig.clf()
            ax = fig.add_subplot(111)
            c = ax.contourf(self.n_energy, self.n_pitch, self.n_wght[ch,:,:], 50)
            fig.colorbar(c)
            ax.set_title('NPA Weight')
            ax.set_ylabel('Pitch')
            ax.set_xlabel('Energy [keV]')
            canvas.show()

    def plot_fida_weights(self,fig,canvas):
        if self._has_fida_wght:
            ch = self.fida_chans[self.fida_chan.get()]
            wl = float(self.lam_val.get())
            ind = np.argmin(np.abs(self.lam-wl))
            fig.clf()
            ax = fig.add_subplot(111)
            c = ax.contourf(self.f_energy,self.f_pitch,self.f_wght[ch,:,:,ind],30)
            fig.colorbar(c)
            ax.set_xlabel('Energy [keV]')
            ax.set_ylabel('Pitch')
            ax.set_title('FIDA Weight')
            canvas.show()

class Neutrals:
    """ Neutrals object that contains plot methods and parameters"""
    def __init__(self ,dir):
        neut_files = glob.glob(dir + '*_neutrals.h5')
        geo_files = glob.glob(dir + '*_geometry.h5')

        self._has_neut = (len(neut_files) > 0)
        self._has_geo = (len(geo_files) > 0)

        if self._has_geo:
            print('Loading geometry')

            if len(geo_files) > 1:
                raise NotImplementedError('Multiple geometry files found')
            else:
                geo = load_dict_from_hdf5(geo_files[0])

            self.beam_name = geo['nbi']['name'].decode('UTF-8')

        if self._has_neut:
            print('Loading neutrals')

            if len(neut_files) > 1:
                raise NotImplementedError('Multiple neutrals files found')
            else:
                neut = load_dict_from_hdf5(neut_files[0])

            # All grids and gridded data to --> (nx, ny, nz)
            self.fdens = neut['fdens'].sum(3).T   # sum over energy state
            self.hdens = neut['hdens'].sum(3).T
            self.tdens = neut['tdens'].sum(3).T
            self.halodens = neut['halodens'].sum(3).T
            self.x_grid = neut['grid']['x_grid'].T     # mach coords
            self.y_grid = neut['grid']['y_grid'].T     # mach coords
            self.z_grid = neut['grid']['z_grid'].T     # mach coords
            self.nx = neut['grid']['nx']
            self.ny = neut['grid']['ny']
            self.nz = neut['grid']['nz']

            # beam coords
            self.x_grid_beam, self.y_grid_beam, self.z_grid_beam = np.meshgrid(neut['grid']['x'], neut['grid']['y'], neut['grid']['z'], indexing='ij')

            # Are beam and machine coordinates the same?
            self.beam_mach_same = np.array_equal(self.x_grid, self.x_grid_beam) and np.array_equal(self.y_grid, self.y_grid_beam) and np.array_equal(self.z_grid, self.z_grid_beam)
        else:
            print('No neutrals found')

        ##Radio Buttons Variable
        self.plot_type = tk.StringVar(value = 'XY')

        ##Checkbox Variables
        self.use_mach_coords = tk.BooleanVar(value = False)
        self.full_on = tk.BooleanVar(value = True)
        self.half_on = tk.BooleanVar(value = True)
        self.third_on = tk.BooleanVar(value = True)
        self.halo_on = tk.BooleanVar(value = True)

    def plot_neutrals(self,fig,canvas):
        full_on = self.full_on.get()
        half_on = self.half_on.get()
        third_on = self.third_on.get()
        halo_on = self.halo_on.get()
        torf = lambda T: 1 if T else 0

        if self._has_neut and (full_on or half_on or third_on or halo_on):
            fig.clf()
            ax = fig.add_subplot(111)
            pt = self.plot_type.get()

            if pt == 'X':
                if self.use_mach_coords.get() and not self.beam_mach_same:
                    # Use machine coords and they're not the same as beam coords

                    ax.set_xlabel('X [cm]')

                    # Need to bin data onto mach regular grid before taking projections
                    fdens_hist = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (self.nx, self.ny), weights=self.fdens.flatten())
                    fdens = fdens_hist[0]

                    # Histogram returns edges of shape (nx+1). Convert to centers
                    xedges = fdens_hist[1]
                    yedges = fdens_hist[2]
                    dx = xedges[1] - xedges[0]
                    x = xedges[0:-1] + dx / 2.

                    hdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.hdens.flatten())[0]
                    tdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.tdens.flatten())[0]
                    halodens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.halodens.flatten())[0]

                    # histogram2d sums weights, need mean
                    fdens = fdens.mean(1) / self.nz
                    hdens = hdens.mean(1) / self.nz
                    tdens = tdens.mean(1) / self.nz
                    halodens = halodens.mean(1) / self.nz
                else:
                    # Use beam coords or beam and machine coords are the same
                    if self.use_mach_coords.get():
                        ax.set_xlabel('X [cm]')
                    elif self.beam_mach_same:
                        ax.set_xlabel('$X = X_{beam}$ [cm]')
                    else:
                        ax.set_xlabel('$X_{beam}$ [cm]')

                    # Use data as is for beam coords or when coord systems are the same
                    x = self.x_grid_beam[:, 0, 0]
                    fdens = self.fdens.mean(1).mean(1)
                    hdens = self.hdens.mean(1).mean(1)
                    tdens = self.tdens.mean(1).mean(1)
                    halodens = self.halodens.mean(1).mean(1)

                if full_on: ax.plot(x, fdens, label = 'Full')
                if half_on: ax.plot(x, hdens, label = 'Half')
                if third_on: ax.plot(x, tdens, label = 'Third')
                if halo_on: ax.plot(x, halodens, label = 'Halo')
                ax.legend()
                ax.set_title('Neutral Density. NB {}'.format(self.beam_name))
                ax.set_ylabel('Mean Density [$cm^{-3}$]')
                canvas.show()

            if pt == 'Y':
                if self.use_mach_coords.get() and not self.beam_mach_same:
                    # Use machine coords and they're not the same as beam coords

                    ax.set_xlabel('Y [cm]')

                    # Need to bin data onto mach regular grid before taking projections
                    fdens_hist = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (self.nx, self.ny), weights=self.fdens.flatten())
                    fdens = fdens_hist[0]

                    # Histogram returns edges of shape (nx+1). Convert to centers
                    xedges = fdens_hist[1]
                    yedges = fdens_hist[2]
                    dx = yedges[1] - yedges[0]
                    x = yedges[0:-1] + dx / 2.

                    hdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.hdens.flatten())[0]
                    tdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.tdens.flatten())[0]
                    halodens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.halodens.flatten())[0]

                    # histogram2d sums weights, need mean
                    fdens = fdens.mean(0) / self.nz
                    hdens = hdens.mean(0) / self.nz
                    tdens = tdens.mean(0) / self.nz
                    halodens = halodens.mean(0) / self.nz
                else:
                    # Use beam coords or beam and machine coords are the same
                    if self.use_mach_coords.get():
                        ax.set_xlabel('Y [cm]')
                    elif self.beam_mach_same:
                        ax.set_xlabel('$Y = Y_{beam}$ [cm]')
                    else:
                        ax.set_xlabel('$Y_{beam}$ [cm]')

                    # Use data as is for beam coords or when coord systems are the same
                    x = self.y_grid_beam[0, :, 0]
                    fdens = self.fdens.mean(0).mean(1)
                    hdens = self.hdens.mean(0).mean(1)
                    tdens = self.tdens.mean(0).mean(1)
                    halodens = self.halodens.mean(0).mean(1)

                if full_on: ax.plot(x, fdens, label = 'Full')
                if half_on: ax.plot(x, hdens, label = 'Half')
                if third_on: ax.plot(x, tdens, label = 'Third')
                if halo_on: ax.plot(x, halodens, label = 'Halo')
                ax.legend()
                ax.set_title('Neutral Density. NB {}'.format(self.beam_name))
                ax.set_ylabel('Mean Density [$cm^{-3}$]')
                canvas.show()

            if pt == 'Z':
                if self.use_mach_coords.get() and not self.beam_mach_same:
                    # Use machine coords and they're not the same as beam coords
                    ax.set_xlabel('Z [cm]')

                    # Need to bin data onto mach regular grid before taking projections
                    fdens_hist = np.histogram2d(self.x_grid.flatten(), self.z_grid.flatten(), bins = (self.nx, self.nz), weights=self.fdens.flatten())
                    fdens = fdens_hist[0]

                    # Histogram returns edges of shape (nx+1). Convert to centers
                    xedges = fdens_hist[1]
                    yedges = fdens_hist[2]
                    dx = yedges[1] - yedges[0]
                    x = yedges[0:-1] + dx / 2.

                    hdens = np.histogram2d(self.x_grid.flatten(), self.z_grid.flatten(), bins = (xedges, yedges), weights=self.hdens.flatten())[0]
                    tdens = np.histogram2d(self.x_grid.flatten(), self.z_grid.flatten(), bins = (xedges, yedges), weights=self.tdens.flatten())[0]
                    halodens = np.histogram2d(self.x_grid.flatten(), self.z_grid.flatten(), bins = (xedges, yedges), weights=self.halodens.flatten())[0]

                    # histogram2d sums weights, need mean
                    fdens = fdens.mean(0) / self.ny
                    hdens = hdens.mean(0) / self.ny
                    tdens = tdens.mean(0) / self.ny
                    halodens = halodens.mean(0) / self.ny
                else:
                    # Use beam coords or beam and machine coords are the same
                    if self.use_mach_coords.get():
                        ax.set_xlabel('Z [cm]')
                    elif self.beam_mach_same:
                        ax.set_xlabel('$Z = Z_{beam}$ [cm]')
                    else:
                        ax.set_xlabel('$Z_{beam}$ [cm]')

                    # Use data as is for beam coords or when coord systems are the same
                    x = self.z_grid_beam[0, 0, :]
                    fdens = self.fdens.mean(0).mean(0)
                    hdens = self.hdens.mean(0).mean(0)
                    tdens = self.tdens.mean(0).mean(0)
                    halodens = self.halodens.mean(0).mean(0)

                if full_on: ax.plot(x, fdens, label = 'Full')
                if half_on: ax.plot(x, hdens, label = 'Half')
                if third_on: ax.plot(x, tdens, label = 'Third')
                if halo_on: ax.plot(x, halodens, label = 'Halo')
                ax.legend()
                ax.set_title('Neutral Density. NB {}'.format(self.beam_name))
                ax.set_ylabel('Mean Density [$cm^{-3}$]')
                canvas.show()

            if pt == 'XY':
                if self.use_mach_coords.get() and not self.beam_mach_same:
                    # Use machine coords and they're not the same as beam coords
                    ax.set_xlabel('X [cm]')
                    ax.set_ylabel('Y [cm]')

                    # Need to bin data onto mach regular grid before taking projections
                    fdens_hist = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (self.nx, self.ny), weights=self.fdens.flatten())
                    fdens = fdens_hist[0]

                    # Histogram returns edges of shape (nx+1). Convert to centers
                    xedges = fdens_hist[1]
                    yedges = fdens_hist[2]
                    dx = xedges[1] - xedges[0]
                    dy = yedges[1] - yedges[0]
                    x = xedges[0:-1] + dx / 2.
                    y = yedges[0:-1] + dy / 2.

                    x, y = np.meshgrid(x, y, indexing='ij')

                    hdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.hdens.flatten())[0]
                    tdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.tdens.flatten())[0]
                    halodens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.halodens.flatten())[0]

                    # histogram2d sums weights, need mean
                    fdens = fdens / self.nz
                    hdens = hdens / self.nz
                    tdens = tdens / self.nz
                    halodens = halodens / self.nz
                else:
                    # Use beam coords or beam and machine coords are the same
                    if self.use_mach_coords.get():
                        ax.set_xlabel('X [cm]')
                        ax.set_ylabel('Y [cm]')
                    elif self.beam_mach_same:
                        ax.set_xlabel('$X = X_{beam}$ [cm]')
                        ax.set_ylabel('$Y = Y_{beam}$ [cm]')
                    else:
                        ax.set_xlabel('$X_{beam}$ [cm]')
                        ax.set_ylabel('$Y_{beam}$ [cm]')

                    # Use data as is for beam coords or when coord systems are the same
                    x = self.x_grid_beam[:, :, 0]
                    y = self.y_grid_beam[:, :, 0]
                    fdens = self.fdens.mean(2)
                    hdens = self.hdens.mean(2)
                    tdens = self.tdens.mean(2)
                    halodens = self.halodens.mean(2)

                dens = fdens * torf(full_on) + hdens * torf(half_on) + tdens * torf(third_on) + halodens * torf(halo_on)

                c = ax.contourf(x, y, dens, 50)
                cb = fig.colorbar(c)
                cb.ax.set_ylabel('[$cm^{-3}$]')
                ax.set_title('Mean Neutral Density. NB {}'.format(self.beam_name))
                canvas.show()

            if pt == 'XZ':
                if self.use_mach_coords.get() and not self.beam_mach_same:
                    # Use machine coords and they're not the same as beam coords
                    ax.set_xlabel('X [cm]')
                    ax.set_ylabel('Z [cm]')

                    # Need to bin data onto mach regular grid before taking projections
                    fdens_hist = np.histogram2d(self.x_grid.flatten(), self.z_grid.flatten(), bins = (self.nx, self.nz), weights=self.fdens.flatten())
                    fdens = fdens_hist[0]

                    # Histogram returns edges of shape (nx+1). Convert to centers
                    xedges = fdens_hist[1]
                    yedges = fdens_hist[2]
                    dx = xedges[1] - xedges[0]
                    dy = yedges[1] - yedges[0]
                    x = xedges[0:-1] + dx / 2.
                    y = yedges[0:-1] + dy / 2.

                    x, y = np.meshgrid(x, y, indexing='ij')

                    hdens = np.histogram2d(self.x_grid.flatten(), self.z_grid.flatten(), bins = (xedges, yedges), weights=self.hdens.flatten())[0]
                    tdens = np.histogram2d(self.x_grid.flatten(), self.z_grid.flatten(), bins = (xedges, yedges), weights=self.tdens.flatten())[0]
                    halodens = np.histogram2d(self.x_grid.flatten(), self.z_grid.flatten(), bins = (xedges, yedges), weights=self.halodens.flatten())[0]

                    # histogram2d sums weights, need mean
                    fdens = fdens / self.ny
                    hdens = hdens / self.ny
                    tdens = tdens / self.ny
                    halodens = halodens / self.ny
                else:
                    # Use beam coords or beam and machine coords are the same
                    if self.use_mach_coords.get():
                        ax.set_xlabel('X [cm]')
                        ax.set_ylabel('Z [cm]')
                    elif self.beam_mach_same:
                        ax.set_xlabel('$X = X_{beam}$ [cm]')
                        ax.set_ylabel('$Z = Z_{beam}$ [cm]')
                    else:
                        ax.set_xlabel('$X_{beam}$ [cm]')
                        ax.set_ylabel('$Z_{beam}$ [cm]')

                    # Use data as is for beam coords or when coord systems are the same
                    x = self.x_grid_beam[:, 0, :]
                    y = self.z_grid_beam[:, 0, :]
                    fdens = self.fdens.mean(1)
                    hdens = self.hdens.mean(1)
                    tdens = self.tdens.mean(1)
                    halodens = self.halodens.mean(1)

                dens = fdens * torf(full_on) + hdens * torf(half_on) + tdens * torf(third_on) + halodens * torf(halo_on)

                c = ax.contourf(x,y,dens,50)
                cb = fig.colorbar(c)
                cb.ax.set_ylabel('[$cm^{-3}$]')
                ax.set_title('Mean Neutral Density. NB {}'.format(self.beam_name))
                canvas.show()

            if pt == 'YZ':
                if self.use_mach_coords.get() and not self.beam_mach_same:
                    # Use machine coords and they're not the same as beam coords
                    ax.set_xlabel('Y [cm]')
                    ax.set_ylabel('Z [cm]')

                    # Need to bin data onto mach regular grid before taking projections
                    fdens_hist = np.histogram2d(self.y_grid.flatten(), self.z_grid.flatten(), bins = (self.ny, self.nz), weights=self.fdens.flatten())
                    fdens = fdens_hist[0]

                    # Histogram returns edges of shape (nx+1). Convert to centers
                    xedges = fdens_hist[1]
                    yedges = fdens_hist[2]
                    dx = xedges[1] - xedges[0]
                    dy = yedges[1] - yedges[0]
                    x = xedges[0:-1] + dx / 2.
                    y = yedges[0:-1] + dy / 2.

                    x, y = np.meshgrid(x, y, indexing='ij')

                    hdens = np.histogram2d(self.y_grid.flatten(), self.z_grid.flatten(), bins = (xedges, yedges), weights=self.hdens.flatten())[0]
                    tdens = np.histogram2d(self.y_grid.flatten(), self.z_grid.flatten(), bins = (xedges, yedges), weights=self.tdens.flatten())[0]
                    halodens = np.histogram2d(self.y_grid.flatten(), self.z_grid.flatten(), bins = (xedges, yedges), weights=self.halodens.flatten())[0]

                    # histogram2d sums weights, need mean
                    fdens = fdens / self.nx
                    hdens = hdens / self.nx
                    tdens = tdens / self.nx
                    halodens = halodens / self.nx
                else:
                    # Use beam coords or beam and machine coords are the same
                    if self.use_mach_coords.get():
                        ax.set_xlabel('Y [cm]')
                        ax.set_ylabel('Z [cm]')
                    elif self.beam_mach_same:
                        ax.set_xlabel('$Y = Y_{beam}$ [cm]')
                        ax.set_ylabel('$Z = Z_{beam}$ [cm]')
                    else:
                        ax.set_xlabel('$Y_{beam}$ [cm]')
                        ax.set_ylabel('$Z_{beam}$ [cm]')

                    # Use data as is for beam coords or when coord systems are the same
                    x = self.y_grid_beam[0, :, :]
                    y = self.z_grid_beam[0, :, :]
                    fdens = self.fdens.mean(0)
                    hdens = self.hdens.mean(0)
                    tdens = self.tdens.mean(0)
                    halodens = self.halodens.mean(0)

                dens = fdens * torf(full_on) + hdens * torf(half_on) + tdens * torf(third_on) + halodens * torf(halo_on)

                c = ax.contourf(x, y, dens, 50)
                cb = fig.colorbar(c)
                cb.ax.set_ylabel('[$cm^{-3}$]')
                ax.set_title('Mean Neutral Density. NB {}'.format(self.beam_name))
                canvas.show()

class Viewer:
    """Class that contains FIDAsim result viewer window"""
    def __init__(self, parent):

        self.load_dir()
        parent.title('FIDAviewer. {}'.format(self.dir))

        #Make MenuBar
        self.MenuBar = tk.Menu(parent)
        parent.config(menu = self.MenuBar)
        self.file = tk.Menu(self.MenuBar, tearoff = False)
        self.file.add_command(label = 'Load Run', command = (lambda: self.load_dir()))
        self.file.add_command(label = 'Quit', command = (lambda: sys.exit()))
        self.MenuBar.add_cascade(label = 'File', menu = self.file, underline = 0)

        #Make Notebook
        self.nb = ttk.Notebook(parent)
        self.spectra_frame = ttk.Frame(self.nb)
        self.npa_frame = ttk.Frame(self.nb)
        self.neutrals_frame = ttk.Frame(self.nb)
        self.weights_frame = ttk.Frame(self.nb)
#        self.imaging_frame = ttk.Frame(self.nb)
        self.nb.add(self.spectra_frame, text = 'Spectra')
        self.nb.add(self.npa_frame ,text = 'NPA')
        self.nb.add(self.neutrals_frame, text = 'Neutrals')
        self.nb.add(self.weights_frame, text = 'Weights')
#        self.nb.add(self.imaging_frame, text = 'Imaging')
        self.nb.pack(side = tk.LEFT , expand = tk.Y, fill = tk.BOTH)
        self.fig = plt.Figure(figsize = (6, 5), dpi = 100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master = parent)
        self.canvas.get_tk_widget().pack(side = tk.RIGHT)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, parent)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)

        #Spectra Frame
        if self.spec._has_spectra:
            ttk.Combobox(self.spectra_frame, textvariable = self.spec.chan,
                         values = list(self.spec.channels.keys())).pack()

            ttk.Checkbutton(self.spectra_frame, text = 'Hide NBI', variable = self.spec.nbi_on,
                            onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.spectra_frame,text = 'Hide FIDA', variable = self.spec.fida_on,
                            onvalue = False,offvalue = True).pack()

            ttk.Checkbutton(self.spectra_frame,text = 'No Bremsstrahlung', variable = self.spec.brems_on,\
            	             onvalue = False,offvalue = True).pack()

            ttk.Checkbutton(self.spectra_frame,text = 'Hide Legend', variable = self.spec.legend_on,\
            	             onvalue = False,offvalue = True).pack()

            ttk.Button(self.spectra_frame,text = 'Plot Spectra',\
            	        command = (lambda: self.spec.plot_spectra(self.fig,self.canvas))).pack(side = tk.TOP,expand = tk.Y,fill = tk.BOTH)

            ttk.Button(self.spectra_frame,text = 'Plot Intensity',\
            	        command = (lambda: self.spec.plot_intensity(self.fig,self.canvas))).pack(side = tk.TOP,expand = tk.Y,fill = tk.BOTH)
        else:
            ttk.Label(self.spectra_frame, text = '\n\nNo spectral data found').pack()

        #NPA Frame
        if self.npa._has_npa:
            ttk.Combobox(self.npa_frame, textvariable = self.npa.chan, values = tuple(self.npa.channels.keys())).pack()

            ttk.Button(self.npa_frame, text = 'Plot Neutral Birth',\
                       command = (lambda: self.npa.plot_neutral_birth(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y,fill = tk.BOTH)

            ttk.Button(self.npa_frame, text = 'Plot Flux',\
                       command = (lambda: self.npa.plot_flux(self.fig, self.canvas))).pack(side = tk.TOP,expand = tk.Y, fill = tk.BOTH)
        else:
            ttk.Label(self.npa_frame, text = '\n\nNo NPA data found').pack()

        #Neutrals Frame
        ttk.Radiobutton(self.neutrals_frame,text = 'Density vs X',variable = self.neut.plot_type,value = 'X').pack()
        ttk.Radiobutton(self.neutrals_frame,text = 'Density vs Y',variable = self.neut.plot_type,value = 'Y').pack()
        ttk.Radiobutton(self.neutrals_frame,text = 'Density vs Z',variable = self.neut.plot_type,value = 'Z').pack()
        ttk.Radiobutton(self.neutrals_frame,text = 'Contour XY',variable = self.neut.plot_type,value = 'XY').pack()
        ttk.Radiobutton(self.neutrals_frame,text = 'Contour XZ',variable = self.neut.plot_type,value = 'XZ').pack()
        ttk.Radiobutton(self.neutrals_frame,text = 'Contour YZ',variable = self.neut.plot_type,value = 'YZ').pack()


        ttk.Checkbutton(self.neutrals_frame,text = 'Use Machine Coordinates', variable = self.neut.use_mach_coords,\
                        onvalue = True,offvalue = False).pack()
        ttk.Checkbutton(self.neutrals_frame,text = 'Hide Full', variable = self.neut.full_on,\
                        onvalue = False,offvalue = True).pack()
        ttk.Checkbutton(self.neutrals_frame,text = 'Hide Half', variable = self.neut.half_on,\
                        onvalue = False,offvalue = True).pack()
        ttk.Checkbutton(self.neutrals_frame,text = 'Hide Third', variable = self.neut.third_on,\
                        onvalue = False,offvalue = True).pack()
        ttk.Checkbutton(self.neutrals_frame,text = 'Hide Halo', variable = self.neut.halo_on,\
                        onvalue = False,offvalue = True).pack()

        ttk.Button(self.neutrals_frame,text = 'Plot',\
                   command = (lambda: self.neut.plot_neutrals(self.fig,self.canvas))).pack(expand = tk.Y,fill = tk.BOTH)

        #Weights Frame
        if self.wght._has_fida_wght:
            ttk.Combobox(self.weights_frame,textvariable = self.wght.fida_chan,\
                         values = tuple(self.wght.fida_chans.keys())).pack()

            tk.Scale(self.weights_frame,orient = tk.HORIZONTAL,length = 200,\
                     from_ = self.wght.wl_min,to = self.wght.wl_max,resolution = self.wght.dlam,variable = self.wght.lam_val).pack()

            ttk.Button(self.weights_frame,text = 'Plot FIDA Weights',\
                       command = (lambda: self.wght.plot_fida_weights(self.fig,self.canvas))).pack(side = tk.TOP,expand = tk.Y,fill = tk.BOTH)
        else:
            ttk.Label(self.weights_frame, text = '\n\nNo FIDA weight data found').pack()

        if self.wght._has_npa_wght:
            ttk.Combobox(self.weights_frame,textvariable = self.wght.npa_chan,\
                         values = tuple(self.wght.npa_chans.keys())).pack()

            ttk.Button(self.weights_frame,text = 'Plot NPA Weights',\
                       command = (lambda: self.wght.plot_npa_weights(self.fig,self.canvas))).pack(side = tk.TOP,expand = tk.Y,fill = tk.BOTH)
        else:
            ttk.Label(self.weights_frame, text = '\n\nNo NPA weight data found').pack()

    def load_nml(self, dir):
        nml = f90nml.read(glob.glob(dir + '*_inputs.dat')[0])['fidasim_inputs']

        return nml

    def load_dir(self):
        self.dir = askdirectory()
        if self.dir == '': self.dir = os.getcwd()
        if self.dir[-1] != '/': self.dir += '/'
#        if self.dir[-1] != os.path.sep: self.dir += os.path.sep

        self.nml = self.load_nml(self.dir)
        self.spec = Spectra(self.dir, self.nml)
        self.npa = NPA(self.dir)
        self.neut = Neutrals(self.dir)
        self.wght = Weights(self.dir)

if __name__ == '__main__':
    root = tk.Tk()
    Viewer(root)
    root.mainloop()

