#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import numpy as np
import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import ttk
import h5py
import f90nml
import collections
#from scipy.integrate import simps
import scipy.integrate as integrate
import fidasim as fs
import scipy.interpolate as interpolate


"""
Todo
----
* consider showing images vs angles instead of distance (ie independent of projection_dist)
* add beam centerline to imaging contour plots
* cannot edit wavelengths until after changing channel and replotting. Why? Fix this.
* fix bad plots of neutrals in mach coords
* clean plots when all data turned off
* with smart h5 reader, could load only info needed to make gui first, then only get data when called, and then save for later use
* take units from files, don't hardcode. Low priority future-proofing
* in taking mean of beam densities, should it only be for non-zero elements? As grid vol --> inf, density --> 0 otherwise
* optimize: can more stuff be loaded only when used? can more stuff be saved and not recalculated (ie set/get)?
* option to change volume element in neutral plotting for better fidelity in going from beam to mach coords
* rerun sample file sim setting all bools to true
* get more intellegent h5 reader to just grab what's needed
* NPA needs work. I haven't used NPA data before - NGB
* currently seems to load neutrals twice. check this and fix

* DONE - add validation of wavelength min and max to not be beyond data. Don't do this
* DONE - can cache pickle files to make much faster. Don't do this.
* DONE - make 'reset wavelength range' button for spectra and imaging frames
* DONE - find out if histogram2d gives left edges or right. A: right side. doing it right
* DONE - separate concept of nbi_on and full_on from has_nb_spec and has_full_spec to separate spectra from imaging frames
* DONE - add another tab to gui "Imaging" w/ "Lens" drop down. Choose spectra and wavelength range to integrate and make contour
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
* DONE - implement multiple matching filenames
* DONE - Make Spectra wl_min and wl_max changeable from gui
"""

def project_image(projection_dist, axes, aperture, data, beam_src, beam_axis):
    """Given several lines of sight and an intensity per LOS, project an image on a plane perpendicular
    to the average LOS axis.

    Parameters
    ----------
    projection_dist : float
        Distance along average LOS axis for projection plane

    axes : array (nchan, 3)
        Normalized axes vectors defining LOS

    aperture : array (3)
        Common location for all LOS (aperature or lens)

    data : array (nchan)
        Data to be projected

    Returns
    -------
    x1_grid : float array (100, 101)
        Relative coordinates for grid_data (perpendicular to average LOS axis)

    x2_grid : float array (100, 101)
        Second set of relative coordinates for grid_data (perpendicular to average LOS axis)

    grid_data : float array (100, 101)
        data interpolated onto a uniform grid on a plane perpendicular to the average LOS axis and
        projection_dist from the aperature

    valid_ic : int array (<=nchan)
        Indeces (relative to an (nchan) array) indicating which LOS made a valid projection

    Todo
    ----
    * Generalize with points (nchan, 3) distinct for each LOS
    """
    avg_los_axis = axes.mean(0)           # (3)
    nchan = axes.shape[0]

    # Find point projection_dist along lens axis (ie point on plane pierced by lens axis line)
    t = np.sqrt(projection_dist ** 2 / np.sum(avg_los_axis ** 2))
    plane_pt1 = aperture + avg_los_axis * t

    # Find any vector perp to avg_los_axis (by crossing w/ any non-colinear vector) to define the plane
    any_vec = np.array([avg_los_axis[0] + 5., avg_los_axis[1], avg_los_axis[2]])   # 5. is arbitrary
    plane_vec1 = np.cross(avg_los_axis, any_vec)

    # Find second plane vector
    plane_vec2 = np.cross(avg_los_axis, plane_vec1)

    # Find two more points to define plane
    plane_pt2 = plane_pt1 + plane_vec1 * 5.         # 5. is arbitrary
    plane_pt3 = plane_pt1 + plane_vec2 * 5.         # 5. is arbitrary

    # Step thru each LOS and find intersection with plane (call them 'target' points)
    target = list()         # locations where LOS hit projection plane
    valid_ic = list()       # indeces (relative to (nchan) array) where LOS makes valid projection
    for ic in range(nchan):
        res = intersect_line_plane(plane_pt1, plane_pt2, plane_pt3, aperture, axes[ic, :])

        if res is None:
            print('Warning: LOS {}, does not intersect projection plane. Ignoring'.format(ic))
        elif len(res) == 2:
            print('Warning: LOS {} lies in projection plane. Ignoring'.format(ic))
        elif len(res) == 3:
            # Intersection is a point
            target.append(res)
            valid_ic.append(ic)

    target = np.array(target)       # (nvalid, 3), all on plane perp to avg_los_axis

    # Remove channels that wont be imaged
    data = data[valid_ic]   # (nvalid)

    # Rotate target locations into coord sys aligned with avg_los_axis
    dis = np.sqrt(np.sum((aperture - plane_pt1) ** 2.0))
    beta = np.arcsin((aperture[2] - plane_pt1[2]) / dis)
    alpha = np.arctan2((plane_pt1[1] - aperture[1]), (plane_pt1[0] - aperture[0]))
    gamma = 0.
    target_rotated = fs.preprocessing.uvw_to_xyz(alpha, beta, gamma, target.T, origin=plane_pt1).T  # (nvalid, 3)

    # Interpolate data onto uniform grid along target plane
    n1d = 100    # no. of grid points in each direction
    x1 = np.linspace(target_rotated[:, 1].min(), target_rotated[:, 1].max(), num = n1d)
    x2 = np.linspace(target_rotated[:, 2].min(), target_rotated[:, 2].max(), num = n1d + 1)
    x1_grid, x2_grid = np.meshgrid(x1, x2, indexing='ij')
    grid_data = interpolate.griddata(target_rotated[:, 1:3], data, (x1_grid, x2_grid), fill_value = 0.)

    # Find two pts on beam centerline in mach coords
    beam_pt1 = beam_src
    beam_pt2 = beam_pt1 + beam_axis * 10.

    # Rotate beam to relative coords
    beam_rotated = fs.preprocessing.uvw_to_xyz(alpha, beta, gamma, np.array([beam_pt1, beam_pt2]).T, origin=plane_pt1)  # (3, 2)
    beam_rotated_pt1 = beam_rotated[1:3, 0]
    beam_rotated_pt2 = beam_rotated[1:3, 1]

#    return x1_grid, x2_grid, grid_data, valid_ic
    return x1, x2, grid_data, beam_rotated_pt1, beam_rotated_pt2 - beam_rotated_pt1

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
        if np.abs(np.dot(line_axis, plane_norm_vec)) < tol:
            # Line and plane are parallel
            if np.abs(np.dot(vec1, plane_norm_vec)) < tol:
                # Line is in the plane. Intersection is the line itself
                return [X4, line_axis]
            else:
                # Line does not intersect plane
                return None
        else:
            # Intersection is a point
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
    uniq_lens_indices : list
        Indeces to locate spectra for each unique lens location

    nlenses : int
        Number of unique len locations

    """
    nchan = lens_loc.shape[0]
    lens_list = [tuple(lens_loc[i,:]) for i in range(nchan)]
    lens_set = set(lens_list)

    uniq_lens_indices = []
    for i,lens in enumerate(lens_list):
        if lens in lens_set:
            uniq_lens_indices.append(i)
            lens_set.remove(lens)

    return uniq_lens_indices, len(uniq_lens_indices)


class Spectra:
    """ Spectra object that contains plot methods and parameters"""
    def __init__(self, nml):
        dir = nml["result_dir"]
        runid = nml["runid"]
        spec_file = os.path.join(dir, runid+'_spectra.h5')
        geo_file = nml["geometry_file"]
        self._has_spectra = os.path.isfile(spec_file)
        self._has_geo = os.path.isfile(geo_file)

        if self._has_spectra:
            print('Loading spectra')

            spec = load_dict_from_hdf5(spec_file)

            self.lam = spec['lambda']
            self.nchan = spec['nchan']
            self.channels_spectra = collections.OrderedDict(('Channel ' + str(i + 1), i) for i in range(self.nchan))

            self.dlam = np.abs(self.lam[1] - self.lam[0])

            # Availability booleans
            self.has_bes = ('full' in spec)
            self.has_fida = ('fida' in spec)
            self.has_brems = ('brems' in spec)

            # Spectra frame variables (with initial values)
            self.wl_min_spectra = tk.StringVar(value = str(np.min(self.lam)))
            self.wl_max_spectra = tk.StringVar(value = str(np.max(self.lam)))
            self.chan_spectra = tk.StringVar(value = 'Channel 1')
            self.bes_on_spectra = tk.BooleanVar(value = self.has_bes)
            self.fida_on_spectra = tk.BooleanVar(value = self.has_fida)
            self.brems_on_spectra = tk.BooleanVar(value = self.has_brems)
            self.legend_on = tk.BooleanVar(value = True)

            # Imaging frame variables (with initial values)
            self.wl_min_imaging = tk.StringVar(value = str(np.min(self.lam)))
            self.wl_max_imaging = tk.StringVar(value = str(np.max(self.lam)))
            self.full_on_imaging = tk.BooleanVar(value = self.has_bes)
            self.half_on_imaging = tk.BooleanVar(value = self.has_bes)
            self.third_on_imaging = tk.BooleanVar(value = self.has_bes)
            self.halo_on_imaging = tk.BooleanVar(value = self.has_bes)
            self.fida_on_imaging = tk.BooleanVar(value = self.has_fida)
            self.brems_on_imaging = tk.BooleanVar(value = self.has_brems)
            self.projection_dist = tk.StringVar(value = 100.)

            if self.has_brems:
                self.brems = spec['brems']

            if self.has_fida:
                self.fida = spec['fida']

            if self.has_bes:
                self.full = spec['full']
                self.half = spec['half']
                self.third = spec['third']
                self.halo = spec['halo']

            if self._has_geo:
                print('Loading geometry')
                geo = load_dict_from_hdf5(geo_file)

                self.lens_loc = geo['spec']['lens']    # (nchan, 3)
                self.lens_axis = geo['spec']['axis']   # (nchan, 3)
                self.beam_src = geo['nbi']['src']
                self.beam_axis = geo['nbi']['axis']

                self.uniq_lens_indeces, nlenses = find_lenses(self.nchan, self.lens_loc)

                self.lenses = collections.OrderedDict(('Lens ' + str(i + 1), i) for i in range(nlenses))
                self.lens = tk.StringVar(value = 'Lens 1')

#                for i in range(nlenses):
#                    print('Lens {}: {}'.format(i + 1, self.lens_loc[self.uniq_lens_indeces[i][0], :]))
            else:
                print('No geometry file found')
        else:
            print('No Spectra File Found')

    def plot_spectra(self, fig, canvas):
        if self._has_spectra:
            ch = self.channels_spectra[self.chan_spectra.get()]
            lam = self.lam

            fig.clf()
            ax = fig.add_subplot(111)

            if self.brems_on_spectra.get():
                if self.has_brems:
                    ax.plot(lam, self.brems[ch, :], label = 'Brems')
                else:
                    print('No brems spectra available')

            if self.bes_on_spectra.get():
                if self.has_bes:
                    ax.plot(lam, self.full[ch, :], label = 'Full')
                    ax.plot(lam, self.half[ch, :], label = 'Half')
                    ax.plot(lam, self.third[ch, :], label = 'Third')
                    ax.plot(lam, self.halo[ch, :], label = 'Halo')
                else:
                    print('No beam spectra available')

            if self.fida_on_spectra.get():
                if self.has_fida:
                    ax.plot(lam, self.fida[ch, :], label = 'Fida')
                else:
                    print('No FIDA spectra available')

            if self.brems_on_spectra.get() or self.fida_on_spectra.get() or self.bes_on_spectra.get():
                if self.legend_on.get():
                    ax.legend()
                ax.set_yscale('log')
                ax.set_xlabel('Wavelength [nm]')
                ax.set_ylabel('$Ph\ /\ (s\ nm\ sr\ m^2)$')
                ax.set_title(self.chan_spectra.get())
                ax.set_xlim([float(self.wl_min_spectra.get()), float(self.wl_max_spectra.get())])
                canvas.show()
            else:
                print('SPECTRA: No Spectra Selected')
        else:
            print('No Spectra File Found')

    def plot_intensity(self, fig, canvas):
        if self._has_spectra:
            w1 = (self.lam >= float(self.wl_min_spectra.get()))
            w2 = (self.lam <= float(self.wl_max_spectra.get()))
            w = np.logical_and(w1, w2)
            intens = integrate.simps(self.fida[:, w], x = self.lam[w], axis = 1)
            ch = range(1, len(intens) + 1)
            fig.clf()
            ax = fig.add_subplot(111)
            ax.plot(ch, intens)
            ax.set_title('FIDA Intensity vs. Channel')
            ax.set_ylabel('$Ph\ /\ (s\ sr\ m^2)$')
            ax.set_xlabel('Channel Number')
            ax.set_yscale('log')
            canvas.show()
        else: print('No Spectra File Found')

    def plot_spec_image(self, fig, canvas):
        """Plot 2D contour of line-integrated spectra excluding brems
        """
        torf = lambda T: 1. if T else 0.

        lens = self.lenses[self.lens.get()]     # this lens index (0 to nlenses-1)
        ch = self.uniq_lens_indeces[lens]       # (this_nchan), indeces for this lens

        full_on = self.full_on_imaging.get()
        half_on = self.half_on_imaging.get()
        third_on = self.third_on_imaging.get()
        halo_on = self.halo_on_imaging.get()
        fida_on = self.fida_on_imaging.get()

        fig.clf()
        ax = fig.add_subplot(111)
        ax.axis('equal')

        if self.has_bes:
            full = self.full[ch, :]
            half = self.half[ch, :]
            third = self.third[ch, :]
            halo = self.halo[ch, :]
        else:
            full = 0.
            half = 0.
            third = 0.
            halo = 0.
            if full_on or half_on or third_on or halo_on:
                print('No beam spectra available')

        if self.has_fida:
            fida = self.fida[ch, :]
        else:
            fida = 0.
            if fida_on:
                print('No FIDA spectra available')

        if (fida_on) or (full_on) or (half_on) or (third_on) or (halo_on):
            spec = full * torf(full_on) + half * torf(half_on) + third * torf(third_on) + \
                   halo * torf(halo_on) + fida * torf(fida_on)

            # Integrate over wavelengths
            w = (self.lam >= float(self.wl_min_imaging.get())) & (self.lam <= float(self.wl_max_imaging.get()))
            spec = integrate.simps(spec[:, w], x = self.lam[w], axis = 1)  # (this_nchan)

            lens_axis = self.lens_axis[ch, :]           # (this_nchan, 3), all LOS axes for this lens
            lens_loc = self.lens_loc[ch[0], :]          # (3), same for all in ch

#            yp_grid, zp_grid, grid_spec, valid_ic = project_image(float(self.projection_dist.get()), lens_axis, lens_loc, spec)
            x1, x2, grid_spec, beam_pt, beam_axis = project_image(float(self.projection_dist.get()), lens_axis,
                                                                  lens_loc, spec, self.beam_src, self.beam_axis)

            # Find where beam hits edges of target plane (Assumes crosses left and right, not general solution)
            t1 = (np.min(x1) - beam_pt[0]) / beam_axis[0]
            t2 = (np.max(x1) - beam_pt[0]) / beam_axis[0]
            beam_pt1 = beam_pt + beam_axis * t1
            beam_pt2 = beam_pt + beam_axis * t2

            # Plot contour
#            c = ax.contourf(yp_grid, zp_grid, grid_spec, 50)
            c = ax.contourf(x1, x2, grid_spec.T, 50)
            cb = fig.colorbar(c)
            cb.ax.set_ylabel('[$Ph\ /\ (s\ sr\ m^2)$]')
            ax.set_title('Intensity\nLens at [{:4.0f},{:4.0f},{:4.0f}]'.format(lens_loc[0], lens_loc[1], lens_loc[2]))
            ax.set_xlabel('X1 [cm]')
            ax.set_ylabel('X2 [cm]')

            # Overplot beam centerline
            ax.plot([beam_pt1[0], beam_pt2[0]], [beam_pt1[1], beam_pt2[1]], color = 'magenta')
#            arr = plt.Arrow(beam_pt1[0], beam_pt1[1], 0.5, 0.5)
#            ax.add_patch(arr)
#            ax.arrow(beam_pt1[0], beam_pt1[1], 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')
            canvas.show()
        else:
            print('No spectra selected to plot')

    def plot_brems_image(self, fig, canvas):
        """Plot 2D contour of line-integrated brems
        """
        lens = self.lenses[self.lens.get()]     # this lens index (0 to nlenses-1)
        ch = self.uniq_lens_indeces[lens]       # (this_nchan), indeces for this lens

        fig.clf()
        ax = fig.add_subplot(111)
        ax.axis('equal')

        if self.has_brems:
            brems = self.brems[ch, :]

            # Integrate over wavelengths
            w = (self.lam >= float(self.wl_min_imaging.get())) & (self.lam <= float(self.wl_max_imaging.get()))
            spec = integrate.simps(brems[:, w], x = self.lam[w], axis = 1)  # (this_nchan)

            lens_axis = self.lens_axis[ch, :]           # (this_nchan, 3), all LOS axes for this lens
            lens_loc = self.lens_loc[ch[0], :]          # (3), same for all in ch (for TAE data)

#            yp_grid, zp_grid, grid_spec, valid_ic = project_image(float(self.projection_dist.get()), lens_axis, lens_loc, spec)
            x1, x2, grid_spec, beam_pt, beam_axis = project_image(float(self.projection_dist.get()), lens_axis,
                                                                  lens_loc, spec, self.beam_src, self.beam_axis)

            # Find where beam hits edges of target plane (Assumes crosses left and right, not general solution)
            t1 = (np.min(x1) - beam_pt[0]) / beam_axis[0]
            t2 = (np.max(x1) - beam_pt[0]) / beam_axis[0]
            beam_pt1 = beam_pt + beam_axis * t1
            beam_pt2 = beam_pt + beam_axis * t2

             # Plot contour
#            c = ax.contourf(yp_grid, zp_grid, grid_spec, 50)
            c = ax.contourf(x1, x2, grid_spec.T, 50)
            cb = fig.colorbar(c)
            cb.ax.set_ylabel('[$Ph\ /\ (s\ sr\ m^2)$]')
            ax.set_title('Intensity\nLens at [{:4.0f},{:4.0f},{:4.0f}]'.format(lens_loc[0], lens_loc[1], lens_loc[2]))
            ax.set_xlabel('X1 [cm]')
            ax.set_ylabel('X2 [cm]')

            # Overplot beam centerline
            ax.plot([beam_pt1[0], beam_pt2[0]], [beam_pt1[1], beam_pt2[1]], color = 'magenta')
            canvas.show()
        else:
            print('No brems spectra available')

    def reset_wave_spectra(self):
        self.wl_min_spectra.set(np.min(self.lam))
        self.wl_max_spectra.set(np.max(self.lam))

    def reset_wave_imaging(self):
        self.wl_min_imaging.set(np.min(self.lam))
        self.wl_max_imaging.set(np.max(self.lam))


class NPA:
    """ NPA object that contains plot methods and parameters"""
    def __init__(self, nml):
        dir = nml["result_dir"]
        runid = nml["runid"]
        npa_file = os.path.join(dir,runid+'_npa.h5')
        wght_file = os.path.join(dir,runid+'_npa_weights.h5')
        neut_file = os.path.join(dir,runid+'_neutrals.h5')
        geo_file = nml["geometry_file"]

        self._has_npa = os.path.isfile(npa_file)
        self._has_wght = os.path.isfile(wght_file)
        self._has_neut = os.path.isfile(neut_file)
        self._has_geo = os.path.isfile(geo_file)

        if self._has_npa:
            print('Loading NPA')
            npa = load_dict_from_hdf5(npa_file)

            self.npa_energy = npa['energy']
            self.npa_flux = npa['flux']
            self.nchan = npa['nchan']
        else:
            print('No NPA file found')

        if self._has_wght:
            print('Loading NPA weights')
            wght = load_dict_from_hdf5(wght_file)

            self.w_energy = wght['energy']
            self.w_flux = wght['flux']
        else:
            print('No NPA weights file found')

        if self._has_neut:
            neut = load_dict_from_hdf5(neut_file)

            self.dens = neut['fdens'].sum(0).sum(0) + neut['hdens'].sum(0).sum(0) + \
                        neut['tdens'].sum(0).sum(0) + neut['halodens'].sum(0).sum(0)
        else:
            print('No neutrals file found')

        if (self._has_npa or self._has_wght):
            self.channels_npa = collections.OrderedDict(('Channel ' + str(i + 1), i) for i in range(0, self.nchan))  # should it be nchan not 3???

        self.chan_npa = tk.StringVar(value = 'Channel 1')

    def plot_neutral_birth(self, fig, canvas):
        if self._has_npa:
            fig.clf()
            ax = fig.add_subplot(111)
            ch = self.channels_npa[self.chan_npa.get()]

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

    def plot_flux(self, fig, canvas):
        if self._has_npa or self._has_wght:
            fig.clf()
            ax = fig.add_subplot(111)
            ch = self.channels_npa[self.chan_npa.get()]
            if self._has_npa:
                ax.step(self.npa_energy,self.npa_flux[ch,:],label = 'MC Flux')
            if self._has_wght:
                ax.plot(self.w_energy,self.w_flux[ch,:],label = 'WF Flux')

            ax.legend()
            ax.set_title('Neutral Flux: '+self.chan_npa.get())
            ax.set_ylabel('Flux')
            ax.set_xlabel('Energy [keV]')
            canvas.show()
        else: print('No NPA file found')


class Weights:
    """ Weights object that contains plot methods and parameters"""
    def __init__(self,nml):
        dir = nml["result_dir"]
        runid = nml["runid"]
        npa_wght_file = os.path.join(dir,runid+'_npa_weights.h5')
        fida_wght_file = os.path.join(dir,runid+'_fida_weights.h5')

        self._has_npa_wght = os.path.isfile(npa_wght_file)
        self._has_fida_wght = os.path.isfile(fida_wght_file)

        if self._has_fida_wght:
            print('Loading FIDA weights')
            fida = load_dict_from_hdf5(fida_wght_file)

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
            npa = load_dict_from_hdf5(npa_wght_file)
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
    def __init__(self, nml):
        dir = nml["result_dir"]
        runid = nml["runid"]
        neut_file = os.path.join(dir,runid+'_neutrals.h5')
        geo_file = nml["geometry_file"]

        self._has_neut = os.path.isfile(neut_file)
        self._has_geo = os.path.isfile(geo_file)

        if self._has_geo:
            print('Loading geometry')
            geo = load_dict_from_hdf5(geo_file)

            self.beam_name = geo['nbi']['name'].decode('UTF-8')
        else:
            print('No geometry file found')

        if self._has_neut:
            print('Loading neutrals')
            neut = load_dict_from_hdf5(neut_file)

            # All grids and gridded data to --> (nx, ny, nz)
            self.fdens = neut['fdens'].sum(3).T        # sum over energy state
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
            print('No neutrals file found')

        ## Radio Buttons Variable
        self.plot_type = tk.StringVar(value = 'XY')

        ## Checkbox Variables
        self.use_mach_coords = tk.BooleanVar(value = False)
        self.full_on_neutrals = tk.BooleanVar(value = True)
        self.half_on_neutrals = tk.BooleanVar(value = True)
        self.third_on_neutrals = tk.BooleanVar(value = True)
        self.halo_on_neutrals = tk.BooleanVar(value = True)
        self.transpose = tk.BooleanVar(value = False)

    def plot_neutrals(self, fig, canvas):
        full_on = self.full_on_neutrals.get()
        half_on = self.half_on_neutrals.get()
        third_on = self.third_on_neutrals.get()
        halo_on = self.halo_on_neutrals.get()
        torf = lambda T: 1 if T else 0

        if self._has_neut:
            if not (full_on or half_on or third_on or halo_on):
                print('No neutrals selected to plot')
            else:
                fig.clf()
                ax = fig.add_subplot(111)

                pt = self.plot_type.get()

                if pt == 'X':
                    if self.use_mach_coords.get() and not self.beam_mach_same:
                        # Use machine coords and they're not the same as beam coords (so must rebin)
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
                        # Use machine coords and they're not the same as beam coords (so must rebin)

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
                        # Use machine coords and they're not the same as beam coords (so must rebin)
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
                        # Use machine coords and they're not the same as beam coords (so must rebin)
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

#                        x, y = np.meshgrid(x, y, indexing='ij')

                        hdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.hdens.flatten())[0]
                        tdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.tdens.flatten())[0]
                        halodens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(), bins = (xedges, yedges), weights=self.halodens.flatten())[0]

                        # histogram2d sums weights, need mean
                        fdens = fdens / self.nz
                        hdens = hdens / self.nz
                        tdens = tdens / self.nz
                        halodens = halodens / self.nz
                    else:
                        # Use data as is for beam coords or when coord systems are the same
                        x = self.x_grid_beam[:, 0, 0]
                        y = self.y_grid_beam[0, :, 0]
                        fdens = self.fdens.mean(2)
                        hdens = self.hdens.mean(2)
                        tdens = self.tdens.mean(2)
                        halodens = self.halodens.mean(2)

                    if self.transpose.get():
                        if self.use_mach_coords.get():
                            ax.set_xlabel('Y [cm]')
                            ax.set_ylabel('X [cm]')
                        elif self.beam_mach_same:
                            ax.set_xlabel('$Y = Y_{beam}$ [cm]')
                            ax.set_ylabel('$X = X_{beam}$ [cm]')
                        else:
                            ax.set_xlabel('$Y_{beam}$ [cm]')
                            ax.set_ylabel('$X_{beam}$ [cm]')
                    else:
                        if self.use_mach_coords.get():
                            ax.set_xlabel('X [cm]')
                            ax.set_ylabel('Y [cm]')
                        elif self.beam_mach_same:
                            ax.set_xlabel('$X = X_{beam}$ [cm]')
                            ax.set_ylabel('$Y = Y_{beam}$ [cm]')
                        else:
                            ax.set_xlabel('$X_{beam}$ [cm]')
                            ax.set_ylabel('$Y_{beam}$ [cm]')

                    dens = fdens * torf(full_on) + hdens * torf(half_on) + tdens * torf(third_on) + halodens * torf(halo_on)

                    ax.axis('equal')
                    if self.transpose.get():
                        c = ax.contourf(y, x, dens, 50)
                    else:
                        c = ax.contourf(x, y, dens.T, 50)
                    cb = fig.colorbar(c)
                    cb.ax.set_ylabel('[$cm^{-3}$]')
                    ax.set_title('Mean Neutral Density. {}'.format(self.beam_name))
                    canvas.show()

                if pt == 'XZ':
                    if self.use_mach_coords.get() and not self.beam_mach_same:
                        # Use machine coords and they're not the same as beam coords (so must rebin)
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

#                        x, y = np.meshgrid(x, y, indexing='ij')

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
                        x = self.x_grid_beam[:, 0, 0]
                        y = self.z_grid_beam[0, 0, :]
                        fdens = self.fdens.mean(1)
                        hdens = self.hdens.mean(1)
                        tdens = self.tdens.mean(1)
                        halodens = self.halodens.mean(1)

                    if self.transpose.get():
                        if self.use_mach_coords.get():
                            ax.set_xlabel('Z [cm]')
                            ax.set_ylabel('X [cm]')
                        elif self.beam_mach_same:
                            ax.set_xlabel('$Z = Z_{beam}$ [cm]')
                            ax.set_ylabel('$X = X_{beam}$ [cm]')
                        else:
                            ax.set_xlabel('$Z_{beam}$ [cm]')
                            ax.set_ylabel('$X_{beam}$ [cm]')
                    else:
                        if self.use_mach_coords.get():
                            ax.set_xlabel('X [cm]')
                            ax.set_ylabel('Z [cm]')
                        elif self.beam_mach_same:
                            ax.set_xlabel('$X = X_{beam}$ [cm]')
                            ax.set_ylabel('$Z = Z_{beam}$ [cm]')
                        else:
                            ax.set_xlabel('$X_{beam}$ [cm]')
                            ax.set_ylabel('$Z_{beam}$ [cm]')

                    dens = fdens * torf(full_on) + hdens * torf(half_on) + tdens * torf(third_on) + halodens * torf(halo_on)

                    ax.axis('equal')
                    if self.transpose.get():
                        c = ax.contourf(y, x, dens, 50)
                    else:
                        c = ax.contourf(x, y, dens.T, 50)
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

#                        x, y = np.meshgrid(x, y, indexing='ij')

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
                        x = self.y_grid_beam[0, :, 0]
                        y = self.z_grid_beam[0, 0, :]
                        fdens = self.fdens.mean(0)
                        hdens = self.hdens.mean(0)
                        tdens = self.tdens.mean(0)
                        halodens = self.halodens.mean(0)

                    if self.transpose.get():
                        if self.use_mach_coords.get():
                            ax.set_xlabel('Z [cm]')
                            ax.set_ylabel('Y [cm]')
                        elif self.beam_mach_same:
                            ax.set_xlabel('$Z = Z_{beam}$ [cm]')
                            ax.set_ylabel('$Y = Y_{beam}$ [cm]')
                        else:
                            ax.set_xlabel('$Z_{beam}$ [cm]')
                            ax.set_ylabel('$Y_{beam}$ [cm]')
                    else:
                        if self.use_mach_coords.get():
                            ax.set_xlabel('Y [cm]')
                            ax.set_ylabel('Z [cm]')
                        elif self.beam_mach_same:
                            ax.set_xlabel('$Y = Y_{beam}$ [cm]')
                            ax.set_ylabel('$Z = Z_{beam}$ [cm]')
                        else:
                            ax.set_xlabel('$Y_{beam}$ [cm]')
                            ax.set_ylabel('$Z_{beam}$ [cm]')

                    dens = fdens * torf(full_on) + hdens * torf(half_on) + tdens * torf(third_on) + halodens * torf(halo_on)

                    ax.axis('equal')
                    if self.transpose.get():
                        c = ax.contourf(y, x, dens, 50)
                    else:
                        c = ax.contourf(x, y, dens.T, 50)
                    cb = fig.colorbar(c)
                    cb.ax.set_ylabel('[$cm^{-3}$]')
                    ax.set_title('Mean Neutral Density. NB {}'.format(self.beam_name))
                    canvas.show()


class Viewer:
    """Class that contains FIDAsim result viewer window"""
    def __init__(self, parent):

        self.load_namelist()
        parent.title('FIDAviewer. {}'.format(self.namelistfile))

        # Make MenuBar
        self.MenuBar = tk.Menu(parent)
        parent.config(menu = self.MenuBar)
        self.file = tk.Menu(self.MenuBar, tearoff = False)
        self.file.add_command(label = 'Load Run', command = (lambda: self.load_namelist()))
        self.file.add_command(label = 'Quit', command = (lambda: sys.exit()))
        self.MenuBar.add_cascade(label = 'File', menu = self.file, underline = 0)

        # Make Notebook
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
            ttk.Combobox(self.spectra_frame, textvariable = self.spec.chan_spectra,
                         values = list(self.spec.channels_spectra.keys())).pack()

            ttk.Checkbutton(self.spectra_frame, text = 'Hide BES', variable = self.spec.bes_on_spectra,
                            onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.spectra_frame,text = 'Hide FIDA', variable = self.spec.fida_on_spectra,
                            onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.spectra_frame,text = 'Hide Bremsstrahlung', variable = self.spec.brems_on_spectra,\
            	             onvalue = False, offvalue = True).pack()

            ttk.Checkbutton(self.spectra_frame, text = 'Hide Legend', variable = self.spec.legend_on,\
            	             onvalue = False, offvalue = True).pack()

            ttk.Label(self.spectra_frame, text = 'Wavelength Min (nm)').pack()
            ttk.Entry(self.spectra_frame, textvariable = self.spec.wl_min_spectra, state = tk.NORMAL, width = 10).pack()

            ttk.Label(self.spectra_frame, text = 'Wavelength Max (nm)').pack()
            ttk.Entry(self.spectra_frame, textvariable = self.spec.wl_max_spectra, state = tk.NORMAL, width = 10).pack()

            ttk.Button(self.spectra_frame, text = 'Reset Wavelength',\
            	        command = (lambda: self.spec.reset_wave_spectra())).pack(side = tk.TOP)

            ttk.Button(self.spectra_frame, text = 'Plot Spectra',\
            	        command = (lambda: self.spec.plot_spectra(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)

            ttk.Button(self.spectra_frame,text = 'Plot Intensity',\
            	        command = (lambda: self.spec.plot_intensity(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)
        else:
            ttk.Label(self.spectra_frame, text = '\n\nNo spectral data found').pack()

        # NPA Frame
        if self.npa._has_npa:
            ttk.Combobox(self.npa_frame, textvariable = self.npa.chan_npa, values = tuple(self.npa.channels_npa.keys())).pack()

            ttk.Button(self.npa_frame, text = 'Plot Neutral Birth',\
                       command = (lambda: self.npa.plot_neutral_birth(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y,fill = tk.BOTH)

            ttk.Button(self.npa_frame, text = 'Plot Flux',\
                       command = (lambda: self.npa.plot_flux(self.fig, self.canvas))).pack(side = tk.TOP,expand = tk.Y, fill = tk.BOTH)
        else:
            ttk.Label(self.npa_frame, text = '\n\nNo NPA data found').pack()

        # Neutrals Frame
        ttk.Radiobutton(self.neutrals_frame, text = 'Density vs X', variable = self.neut.plot_type, value = 'X').pack()
        ttk.Radiobutton(self.neutrals_frame, text = 'Density vs Y', variable = self.neut.plot_type, value = 'Y').pack()
        ttk.Radiobutton(self.neutrals_frame, text = 'Density vs Z', variable = self.neut.plot_type, value = 'Z').pack()
        ttk.Radiobutton(self.neutrals_frame, text = 'Contour XY', variable = self.neut.plot_type, value = 'XY').pack()
        ttk.Radiobutton(self.neutrals_frame, text = 'Contour XZ', variable = self.neut.plot_type, value = 'XZ').pack()
        ttk.Radiobutton(self.neutrals_frame, text = 'Contour YZ', variable = self.neut.plot_type, value = 'YZ').pack()


        ttk.Checkbutton(self.neutrals_frame, text = 'Use Machine Coordinates', variable = self.neut.use_mach_coords,\
                        onvalue = True, offvalue = False).pack()

        ttk.Checkbutton(self.neutrals_frame, text = 'Hide Full', variable = self.neut.full_on_neutrals,\
                        onvalue = False, offvalue = True).pack()

        ttk.Checkbutton(self.neutrals_frame, text = 'Hide Half', variable = self.neut.half_on_neutrals,\
                        onvalue = False, offvalue = True).pack()

        ttk.Checkbutton(self.neutrals_frame, text = 'Hide Third', variable = self.neut.third_on_neutrals,\
                        onvalue = False, offvalue = True).pack()

        ttk.Checkbutton(self.neutrals_frame, text = 'Hide Halo', variable = self.neut.halo_on_neutrals,\
                        onvalue = False, offvalue = True).pack()

        ttk.Checkbutton(self.neutrals_frame, text = 'Transpose', variable = self.neut.transpose,\
                        onvalue = True, offvalue = False).pack()

        ttk.Button(self.neutrals_frame, text = 'Plot',\
                   command = (lambda: self.neut.plot_neutrals(self.fig, self.canvas))).pack(expand = tk.Y, fill = tk.BOTH)

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

            ttk.Label(self.imaging_frame, text = 'Wavelength Min (nm)').pack()
            ttk.Entry(self.imaging_frame, textvariable = self.spec.wl_min_imaging, state = tk.NORMAL, width = 10).pack()

            ttk.Label(self.imaging_frame, text = 'Wavelength Max (nm)').pack()
            ttk.Entry(self.imaging_frame, textvariable = self.spec.wl_max_imaging, state = tk.NORMAL, width = 10).pack()

            ttk.Button(self.imaging_frame, text = 'Reset Wavelength',\
            	        command = (lambda: self.spec.reset_wave_imaging())).pack(side = tk.TOP)

            ttk.Button(self.imaging_frame, text = 'Plot Image',\
            	        command = (lambda: self.spec.plot_spec_image(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)

            ttk.Button(self.imaging_frame, text = 'Plot Brems',\
            	        command = (lambda: self.spec.plot_brems_image(self.fig, self.canvas))).pack(side = tk.TOP, expand = tk.Y, fill = tk.BOTH)

            ttk.Label(self.imaging_frame, text = 'Projection Distance (cm)').pack()
            ttk.Entry(self.imaging_frame, textvariable = self.spec.projection_dist, state = tk.NORMAL, width = 10).pack()
        else:
            ttk.Label(self.imaging_frame, text = '\n\nNo imaging data found').pack()

    def read_nml(self, filename):
        nml = f90nml.read(filename)['fidasim_inputs']

        # Use nml dir if results_dir is invalid
        if not os.path.isdir(nml['result_dir']):
            nml['result_dir'] = os.path.dirname(filename)

        if not os.path.isfile(nml['geometry_file']):
            nml['geometry_file'] = os.path.join(nml['result_dir'],os.path.basename(nml['geometry_file']))

        return nml

    def load_namelist(self):
        self.namelistfile = askopenfilename(initialdir=os.path.expanduser('~'),filetypes=[('Namelist Files','*.dat')])
        self.nml = self.read_nml(self.namelistfile)
        self.spec = Spectra(self.nml)
        self.npa = NPA(self.nml)
        self.neut = Neutrals(self.nml)
        self.wght = Weights(self.nml)

if __name__ == '__main__':
    root = tk.Tk()
    Viewer(root)
    root.mainloop()

