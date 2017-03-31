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
#from scipy.integrate import simps
import scipy.integrate as integrate


"""
Todo
----
* with smart h5 reader, could load only info needed to make gui first, then only get data when called, and then save for later use
* can cache pickle files to make much faster
* add validation of wavelength min and max to not be beyond data
* take units from files, don't hardcode. Low priority future-proofing
* in taking mean of beam densities, should it only be for non-zero elements? As grid vol --> inf, density --> 0 otherwise
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


def vec2vec_rotate(vec1, vec2):
    """Returns the rotation matrix to rotate from vec1 to parallel to vec2.

    From the cited paper:
    "We have a unit vector, f, that we wish to rotate to another
    unit vector, t, by rotation in a plane containing both; in other words, we
    seek a rotation matrix R(f,t) such that R(f, t)f = t."

    Notes
    -----
    See paper: "Efficiently Building a Matrix to Rotate One Vector to Another"
    http://www.cs.brown.edu/~jfh/papers/Moller-EBA-1999/main.htm

    Error arises when vec1 has two equal length components (cannot define p) (exception created for vec1 on an axis)
    Error arises when vec2 is on an axis and p is the same axis (divide by zero)

    Parameters
    ----------
    vec1 : float array (3)
        Original vector

    vec2 : float array (3)
        Vector to rotate vec1 to

    History
    -------
    2013-11-03 Nathan Bolte copied from Mike Van Zeeland's get_beam_rays.pro
    2013-12-01 NGB changed name from vec2vec_rotation to vec2vec_rotate to avoid letting IDL use copy from
               get_beam_rays.pro, which it was doing.
    2014-01-27 NGB added reform of inputs. Changed from () to [] for arrays.
    2015-01-20 NGB added error output
    2016-12-21 NGB converted to python. Changed variable names to those in paper.

    To Do
    -----
    * Generalize so that either vec in can be any axis
    * DONE - Are the if statements mutally exclusive? If so, make them if, elif, else
    """
    # Make np arrays
    f = np.array(vec1, dtype=float, ndmin=1)
    t = np.array(vec2, dtype=float, ndmin=1)

    # Make sure inputs are only 3 long
    if (f.ndim != t.ndim != 1) or (f.size != t.size != 3):
        raise ValueError('Inputs must be 1D arrays or lists of length 3')

    # Normalize input vectors
    f = f / np.sqrt(np.sum(f ** 2))
    t = t / np.sqrt(np.sum(t ** 2))

    f_abs = np.abs(f)

    # Define p vector
    if (f_abs[0] < f_abs[1]) and (f_abs[0] < f_abs[2]):
        p = [1., 0., 0.]

    elif (f_abs[1] < f_abs[0]) and (f_abs[1] < f_abs[2]):
        p = [0., 1., 0.]

    elif (f_abs[2] < f_abs[0]) and (f_abs[2] < f_abs[1]):
        p = [0., 0., 1.]

    # This allows vec1 to be an axis and p to still be defined
    # will choose y axis if f=xhat
    elif np.sum(f * [1., 0., 0.]) == 1.:
        p = [0., 1., 0.]

    # will choose x axis if f=yhat
    elif np.sum(f * [0., 1., 0.]) == 1.:
        p = [1., 0., 0.]

    # will choose z axis if f=zhat
    elif np.sum(f * [0., 0., 1.]) == 1.:
        p = [0., 0., 1.]

    else:
        raise ValueError('p undefined')

    u = p - f
    v = p - t

    # Dot products
    uu  = np.sum(u * u)
    vv  = np.sum(v * v)
    uv  = np.sum(u * v)

    if uu * vv == 0.:
        raise ValueError('Zero in denominator')

    # Rotation matrix output container
    rot_mat = np.zeros((3,3))

    # Assemble rotation matrix
    for i in range(3):
        for j in range(3):
            if i == j:
                dij = 1.
            else:
                dij = 0.

            rot_mat[i, j] = dij - (2. / uu) * u[i] * u[j] - (2. / vv) * v[i] * v[j] + (4. * uv) / (uu * vv) * v[i] * u[j]

    return rot_mat


def intersect_line_plane(plane_pt1, plane_pt2, plane_pt3, line_pt, line_axis):
        '''Calculate the intersection location between line and plane

        Parameters
        ----------
        Plane object
            Plane to find intersection with this line

        Returns
        -------
        list or None
            Two element list: point and axis of line itself (ie line is in plane)
            Three element list: Coordinates of intersection point
            None: Line does not intersect plane

        Notes
        -----
        Not implemented for multiple lines or planes

        * For testing for cases in line-plane intersection, see [1]_
        * For the cases where the line-plane intersection is a point, see [2]_

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection#Algebraic_form
        .. [2] http://mathworld.wolfram.com/Line-PlaneIntersection.html
        '''
        # other = plane, self = line

        X1 = plane_pt1
        X2 = plane_pt2
        X3 = plane_pt3
        X4 = line_pt
        line_axis = line_axis

        # Vector normal to plane
        plane_norm_vec = np.cross(X1 - X2, X3 - X2)
        plane_norm_vec /= np.linalg.norm(plane_norm_vec)

        # Avoid using same point on line and plane. Just move further along line (arbitrarily let t = 1.)
        if np.array_equal(X1, X4):
            X4 = X4 + line_axis * 1.

        # Test for different cases.
        # Since vec1, plane_norm_vec, and line_axis are all normalized, the following dot products are [0, 1]. So can
        # use a tolerance instead of comparing to zero.
        tol = 1e-15
        vec1 = (X4 - X1)
        vec1 /= np.linalg.norm(vec1)
#        print()
        if np.abs(np.dot(line_axis, plane_norm_vec)) < tol:
            # Line and plane are parallel
#            print('Line and plane are parallel', np.dot(line_axis, plane_norm_vec))
            if np.abs(np.dot(vec1, plane_norm_vec)) < tol:
                # Line is in the plane. Intersection is the line itself
#                print('Line is in the plane.', np.abs(np.dot(vec1, plane_norm_vec)))
                return [X4, line_axis]
            else:
                # Line does not intersect plane
#                print('Line does not intersect plane.' , np.abs(np.dot(vec1, plane_norm_vec)))
                return None
        else:
            # Intersection is a point
#            print('Intersection is a point')
            mat1 = np.ones((4, 4), dtype=float)
            mat1[1:4, 0] = X1
            mat1[1:4, 1] = X2
            mat1[1:4, 2] = X3
            mat1[1:4, 3] = X4

            mat2 = np.copy(mat1)
            mat2[0, 3] = 0.
            mat2[1:4, 3] = line_axis

            t = -np.linalg.det(mat1) / np.linalg.det(mat2)

            x = X4[0] + line_axis[0] * t
            y = X4[1] + line_axis[1] * t
            z = X4[2] + line_axis[2] * t

            return [x, y, z]


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


def find_lenses(nchan, lens_loc):
    """Find locations for unique lenses in fidasim run

    Parameters
    ----------
    nchan : int
        Number of spectral channels (lines of sight)

    lens_loc : 2D array
        Cartesian coords of lenses in machine coords, (nchan, 3)

    Returns
    -------
    uniq_lens_indeces : list
        Indeces to locate spectra for each unique lens location

    nlenses : int
        Number of unique len locations

    Todo
    ----
    * Can't seem to do w/ np.isclose. Make this work
    """
    uniq_lens_indeces = list()
    master_ind = np.arange(nchan)
    nlos = 0
    ic = 0
    iter_count = -1
    while True:
        iter_count += 1
        this_lens_loc = lens_loc[ic, :]
        w = (lens_loc[:, 0] == this_lens_loc[0]) & (lens_loc[:, 1] == this_lens_loc[1]) & (lens_loc[:, 2] == this_lens_loc[2])
        uniq_lens_indeces.append(master_ind[w])
        nlos += uniq_lens_indeces[-1].size
        if (nlos >= nchan) or (iter_count >= nchan):
            break
        else:
            # next index not in w that hasn't been covered yet (ie, still need to examine)
            ic = np.min(np.setdiff1d(master_ind, np.array(uniq_lens_indeces).flatten()))
    nlenses = len(uniq_lens_indeces)

    return uniq_lens_indeces, nlenses


class Spectra:
    """ Spectra object that contains plot methods and parameters"""
    def __init__(self, dir, nml):
        spec_files = glob.glob(dir + '*_spectra.h5')
        geo_files = glob.glob(dir + '*_geometry.h5')
        self._has_spectra = (len(spec_files) > 0)
        self._has_geo = (len(geo_files) > 0)

        if self._has_spectra:
            print('Loading spectra')

            if len(spec_files) > 1:
                raise NotImplementedError('Multiple spectra files found')
            else:
                spec = load_dict_from_hdf5(spec_files[0])

            self.lam = spec['lambda']
            self.nchan = spec['nchan']
            self.channels = collections.OrderedDict(('Channel ' + str(i + 1), i) for i in range(self.nchan))

            self.dlam = np.abs(self.lam[1] - self.lam[0])

            # Spectra frame variables
            self.wl_min = tk.StringVar(value = str(np.min(self.lam)))
            self.wl_max = tk.StringVar(value = str(np.max(self.lam)))
            self.chan = tk.StringVar(value = 'Channel 1')
            self.nbi_on = tk.BooleanVar(value = nml['calc_bes'] > 0)
            self.fida_on = tk.BooleanVar(value = nml['calc_fida'] > 0)
            self.brems_on = tk.BooleanVar(value = nml['calc_brems'] > 0)
            self.legend_on = tk.BooleanVar(value = True)

            # Imaging frame checkbox variables
            self.wl_min_imaging = tk.StringVar(value = str(np.min(self.lam)))
            self.wl_max_imaging = tk.StringVar(value = str(np.max(self.lam)))
            self.full_on_imaging = tk.BooleanVar(value = nml['calc_bes'] > 0)
            self.half_on_imaging = tk.BooleanVar(value = nml['calc_bes'] > 0)
            self.third_on_imaging = tk.BooleanVar(value = nml['calc_bes'] > 0)
            self.halo_on_imaging = tk.BooleanVar(value = nml['calc_bes'] > 0)
            self.fida_on_imaging = tk.BooleanVar(value = nml['calc_fida'] > 0)

            if self.brems_on.get():
                self.brems = spec['brems']
            else:
                self.brems = None

            if self.fida_on.get():
                self.fida = spec['fida']
            else:
                self.fida = None

            if self.nbi_on.get():
                self.full = spec['full']
                self.half = spec['half']
                self.third = spec['third']
                self.halo = spec['halo']
            else:
                self.full = None

            if self._has_geo:
                print('Loading geometry')

                if len(geo_files) > 1:
                    raise NotImplementedError('Multiple geometry files found')
                else:
                    geo = load_dict_from_hdf5(geo_files[0])

                self.lens_loc = geo['spec']['lens']    # (nchan, 3)
                self.lens_axis = geo['spec']['axis']   # (nchan, 3)

                self.uniq_lens_indeces, nlenses = find_lenses(self.nchan, self.lens_loc)

                self.lenses = collections.OrderedDict(('Lens ' + str(i + 1), i) for i in range(nlenses))
                self.lens = tk.StringVar(value = 'Lens 1')
            else:
                print('No geometry file found')
        else:
            print('No spectra found')

    def plot_spectra(self, fig, canvas):
        if self._has_spectra:
            ch = self.channels[self.chan.get()]
            lam = self.lam

            fig.clf()
            ax = fig.add_subplot(111)

            if self.brems_on.get():
                if self.brems is None:
                    print('No brems spectra found')
                else:
                    brems = self.brems[ch, :]
                    ax.plot(lam, brems, label = 'Brems')

            if self.nbi_on.get():
                if self.full is None:
                    print('No beam spectra found')
                else:
                    full = self.full[ch, :]
                    half = self.half[ch, :]
                    third = self.third[ch, :]
                    halo = self.halo[ch, :]

                    ax.plot(lam, full, label = 'Full')
                    ax.plot(lam, half, label = 'Half')
                    ax.plot(lam, third, label = 'Third')
                    ax.plot(lam, halo, label = 'Halo')

            if self.fida_on.get():
                if self.fida is None:
                    print('No FIDA spectra found')
                else:
                    fida = self.fida[ch, :]
                    ax.plot(lam, fida, label = 'Fida')

            if self.brems_on.get() or self.fida_on.get() or self.nbi_on.get():
                if self.legend_on.get(): ax.legend()
                ax.set_yscale('log')
                ax.set_xlabel('Wavelength [nm]')
                ax.set_ylabel('$Ph\ /\ (s\ nm\ sr\ m^2)$')
                ax.set_title(self.chan.get())
                ax.set_xlim([float(self.wl_min.get()), float(self.wl_max.get())])
                canvas.show()
            else:
                print('SPECTRA: No Spectra Selected')
        else:
            print('SPECTRA: No file')

    def plot_intensity(self, fig, canvas):
        if self._has_spectra:
            w1 = (self.lam >= float(self.wl_min.get()))
            w2 = (self.lam <= float(self.wl_max.get()))
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

    def plot_spec_image(self, fig, canvas):
        """Plot 2D contour of line-integrated spectra
        """
        torf = lambda T: 1. if T else 0.

        lens = self.lenses[self.lens.get()]     # this lens index (0 to nlenses-1)
        ch = self.uniq_lens_indeces[lens]       # (this_nchan), indeces for this lens
        this_nchan = ch.size                    # number of channels for this lens

        fig.clf()
        ax = fig.add_subplot(111)

        if self.full_on_imaging.get():
            if self.full is None:
                print('No beam spectra found')
                full = 0.
                half = 0.
                third = 0.
                halo = 0.
            else:
                full = self.full[ch, :]
                half = self.half[ch, :]
                third = self.third[ch, :]
                halo = self.halo[ch, :]

        if self.fida_on_imaging.get():
            if self.fida is None:
                print('No FIDA spectra found')
                fida = 0.
            else:
                fida = self.fida[ch, :]

        if (self.fida is not None) or (self.full is not None):
            spec = full * torf(self.full_on_imaging.get()) + half * torf(self.half_on_imaging.get()) + \
                   third * torf(self.third_on_imaging.get()) + halo * torf(self.halo_on_imaging.get()) + \
                   fida * torf(self.fida_on_imaging.get())

            # Integrate over wavelengths
            w = (self.lam >= float(self.wl_min_imaging.get())) & (self.lam <= float(self.wl_max_imaging.get()))
            spec = integrate.simps(spec[:, w], x = self.lam[w], axis = 1)  # (this_nchan)

            # General method: Project views onto plane projection_dist away from lens along avg lens axis
            projection_dist = 100.                      # arbitary for now, make tk variable
            lens_axis = self.lens_axis[ch, :]           # (this_nchan, 3)
            lens_axis_avg = lens_axis.mean(0)           # (3)
            lens_loc = self.lens_loc[ch[0], :]          # (3), same for all in ch

#            print(lens_axis.shape, lens_axis_avg.shape, lens_loc.shape, spec.shape)

            # Find point projection_dist along lens axis (ie point on plane pierced by lens axis line)
            t = np.sqrt(projection_dist ** 2 / np.sum(lens_axis_avg ** 2))
            plane_pt1 = lens_loc + lens_axis_avg * t

            # Find any vector perp to lens_axis_avg (by crossing w/ any non-colinear vector) to define the plane
            any_vec = np.array([lens_axis_avg[0] + 5., lens_axis_avg[1], lens_axis_avg[2]])   # 5. is arbitrary
            plane_vec1 = np.cross(lens_axis_avg, any_vec)

            # Find second plane vector
            plane_vec2 = np.cross(lens_axis_avg, plane_vec1)

            # Find two more points to define plane
            plane_pt2 = plane_pt1 + plane_vec1 * 5.         # 5. is arbitrary
            plane_pt3 = plane_pt1 + plane_vec2 * 5.         # 5. is arbitrary

            # Step thru each channel for this lens and find intersection with plane (call 'target' points)
            target = list()
            valid_ic = list()
            for ic in range(this_nchan):
                res = intersect_line_plane(plane_pt1, plane_pt2, plane_pt3, lens_loc, lens_axis[ic, :])

                if res is None:
                    msg = 'Warning: {}, channel {}, overall channel {}, does not intersect projection plane. Ignoring.'
                    print(msg.format(self.lens.get(), ic, ch[ic] + 1))
                elif len(res) == 2:
                    msg = 'Warning: {}, channel {}, overall channel {}, lies in projection plane. Ignoring.'
                    print(msg.format(self.lens.get(), ic, ch[ic] + 1))
                elif len(res) == 3:
                    # Intersection is a point
                    target.append(res)
                    valid_ic.append(ic)

            nvalid = len(valid_ic)
            target = np.array(target)       # (nvalid, 3), all on plane perp to lens_axis_avg

            # are target points all actually on one plane --> YES!
#            for i in range(nvalid - 1):
#                print(np.dot(target[i, :] - target[i + 1, :], lens_axis_avg))

            # Remove channels that wont be imaged
            spec = spec[valid_ic]   # (nvalid)
            ch = ch[valid_ic]


            # Now need to make uniform grid in plane and triangulate spec onto that grid
            rot_mat = vec2vec_rotate([1., 0., 0.], lens_axis_avg)

            target_rotated = np.copy(target)

            # Shift origin
#            target_rotated = target_rotated + np.tile(plane_pt1, (nvalid, 1))

            # Apply rotation matrix
            target_rotated = np.dot(rot_mat, target_rotated.T).T

            # Shift origin
            target_rotated = target_rotated + np.tile(plane_pt1, (nvalid, 1))

#            print(target_rotated)
#
#            print(np.dot(target_rotated[1, :] - target_rotated[0, :], lens_axis_avg))

#            for i in range(nvalid - 1):
#                print(np.dot(target_rotated[i, :] - target_rotated[i + 1, :], [1., 0., 0.]))


            # Plot contour
#            c = ax.contourf(x, y, spec, 50)
#            cb = fig.colorbar(c)
#            cb.ax.set_ylabel('[$Ph\ /\ (s\ sr\ m^2)$]')
#            ax.set_title('Intensity')
#            canvas.show()


#            if self.legend_on.get(): ax.legend()
#            ax.set_yscale('log')
#            ax.set_xlabel('Wavelength [nm]')
#            ax.set_ylabel('$Ph\ /\ (s\ nm\ sr\ m^2)$')
#            ax.set_title(self.chan.get())
#            ax.set_xlim([float(self.wl_min.get()), float(self.wl_max.get())])
#            canvas.show()

    def plot_brems_image(self, fig, canvas):
        print('Not implemented yet')
#        if self.brems_on.get():
#            if self.brems is None:
#                print('No brems spectra found')
#            else:
#                brems = self.brems[ch, :]
#                ax.plot(lam, brems, label = 'Brems')


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
            self.fida_chans = collections.OrderedDict(('Channel '+str(i+1),i) for i in range(0,self.f_chan))
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
            self.npa_chans = collections.OrderedDict(('Channel ' + str(i + 1), i) for i in range(0, self.n_nchan))

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
    def __init__(self, dir):
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
        else:
            print('No geometry file found')

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

        ## Radio Buttons Variable
        self.plot_type = tk.StringVar(value = 'XY')

        ## Checkbox Variables
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
        self.imaging_frame = ttk.Frame(self.nb)
        self.nb.add(self.spectra_frame, text = 'Spectra')
        self.nb.add(self.npa_frame ,text = 'NPA')
        self.nb.add(self.neutrals_frame, text = 'Neutrals')
        self.nb.add(self.weights_frame, text = 'Weights')
        self.nb.add(self.imaging_frame, text = 'Imaging')
        self.nb.pack(side = tk.LEFT , expand = tk.Y, fill = tk.BOTH)
        self.fig = plt.Figure(figsize = (6, 5), dpi = 100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master = parent)
        self.canvas.get_tk_widget().pack(side = tk.RIGHT)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, parent)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)

        # Spectra Frame
        if self.spec._has_spectra:
            ttk.Combobox(self.spectra_frame, textvariable = self.spec.chan,
                         values = list(self.spec.channels.keys())).pack()

            ttk.Checkbutton(self.spectra_frame, text = 'Hide NBI', variable = self.spec.nbi_on,
                            onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.spectra_frame,text = 'Hide FIDA', variable = self.spec.fida_on,
                            onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.spectra_frame,text = 'Hide Bremsstrahlung', variable = self.spec.brems_on,\
            	             onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.spectra_frame, text = 'Hide Legend', variable = self.spec.legend_on,\
            	             onvalue = False, offvalue = True).pack()

            ttk.Label(self.spectra_frame, text = 'Wavelength Min.').pack()
            ttk.Entry(self.spectra_frame, textvariable = self.spec.wl_min, state = tk.NORMAL, width = 10).pack()

            ttk.Label(self.spectra_frame, text = 'Wavelength Max.').pack()
            ttk.Entry(self.spectra_frame, textvariable = self.spec.wl_max, state = tk.NORMAL, width = 10).pack()

            ttk.Button(self.spectra_frame, text = 'Plot Spectra',\
            	        command = (lambda: self.spec.plot_spectra(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)

            ttk.Button(self.spectra_frame,text = 'Plot Intensity',\
            	        command = (lambda: self.spec.plot_intensity(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)
        else:
            ttk.Label(self.spectra_frame, text = '\n\nNo spectral data found').pack()

        # NPA Frame
        if self.npa._has_npa:
            ttk.Combobox(self.npa_frame, textvariable = self.npa.chan, values = tuple(self.npa.channels.keys())).pack()

            ttk.Button(self.npa_frame, text = 'Plot Neutral Birth',\
                       command = (lambda: self.npa.plot_neutral_birth(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y,fill = tk.BOTH)

            ttk.Button(self.npa_frame, text = 'Plot Flux',\
                       command = (lambda: self.npa.plot_flux(self.fig, self.canvas))).pack(side = tk.TOP,expand = tk.Y, fill = tk.BOTH)
        else:
            ttk.Label(self.npa_frame, text = '\n\nNo NPA data found').pack()

        # Neutrals Frame
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

        # Weights Frame
        if self.wght._has_fida_wght:
            ttk.Combobox(self.weights_frame,textvariable = self.wght.fida_chan,\
                         values = tuple(self.wght.fida_chans.keys())).pack()

            tk.Scale(self.weights_frame,orient = tk.HORIZONTAL, length = 200,\
                     from_ = self.wght.wl_min, to = self.wght.wl_max, resolution = self.wght.dlam, variable = self.wght.lam_val).pack()

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

        # Imaging frame
        if self.spec._has_spectra and self.spec._has_geo:
            ttk.Combobox(self.imaging_frame, textvariable = self.spec.lens,
                         values = list(self.spec.lenses.keys())).pack()

            ttk.Checkbutton(self.imaging_frame,text = 'Exclude FIDA', variable = self.spec.fida_on_imaging,
                            onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.imaging_frame,text = 'Exclude Full', variable = self.spec.full_on_imaging,
                            onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.imaging_frame,text = 'Exclude Half', variable = self.spec.half_on_imaging,
                            onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.imaging_frame,text = 'Exclude Third', variable = self.spec.third_on_imaging,
                            onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.imaging_frame,text = 'Exclude Halo', variable = self.spec.halo_on_imaging,
                            onvalue = False, offvalue = True).pack()

            ttk.Label(self.imaging_frame, text = 'Wavelength Min.').pack()
            ttk.Entry(self.imaging_frame, textvariable = self.spec.wl_min_imaging, state = tk.NORMAL, width = 10).pack()

            ttk.Label(self.imaging_frame, text = 'Wavelength Max.').pack()
            ttk.Entry(self.imaging_frame, textvariable = self.spec.wl_max_imaging, state = tk.NORMAL, width = 10).pack()

            ttk.Button(self.imaging_frame, text = 'Plot Image',\
            	        command = (lambda: self.spec.plot_spec_image(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)

            ttk.Button(self.imaging_frame, text = 'Plot Brems',\
            	        command = (lambda: self.spec.plot_brems_image(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)

        else:
            ttk.Label(self.imaging_frame, text = '\n\nNo imaging data found').pack()

    def load_nml(self, dir):
        nml = f90nml.read(glob.glob(dir + '*_inputs.dat')[0])['fidasim_inputs']

        return nml

    def load_dir(self):
        self.dir = askdirectory()
        if self.dir == '': self.dir = os.getcwd()
        if self.dir[-1] != '/': self.dir += '/'

        self.nml = self.load_nml(self.dir)
        self.spec = Spectra(self.dir, self.nml)
        self.npa = NPA(self.dir)
        self.neut = Neutrals(self.dir)
        self.wght = Weights(self.dir)

if __name__ == '__main__':
    root = tk.Tk()
    Viewer(root)
    root.mainloop()

