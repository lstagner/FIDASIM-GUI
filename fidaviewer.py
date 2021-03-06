#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
#import glob
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
#import fidasim as fs
import scipy.interpolate as interpolate
#import numpy as np
from scipy.spatial.distance import cdist


"""
Todo
----
* use http://nbviewer.jupyter.org/gist/tillahoffmann/f844bce2ec264c1c8cb5 for binning neutral density and spec (project_image)
* change images to use angles instead of distance (ie independent of projection_dist)
* add beam centerline to imaging contour plots
* cannot edit wavelengths until after changing channel and replotting. Why? Fix this.
* fix bad plots of neutrals in mach coords by using kde
* clean plots when all data turned off
* with smart h5 reader, could load only info needed to make gui first, then only get data when called, and then save for later use
* take units from files, don't hardcode. Low priority future-proofing
* in taking mean of beam densities, should it only be for non-zero elements? As grid vol --> inf, density --> 0 otherwise
* optimize: can more stuff be loaded only when used? can more stuff be saved and not recalculated (ie set/get)?
* option to change volume element in neutral plotting for better fidelity in going from beam to mach coords
* get more intellegent h5 reader to just grab what's needed
* NPA needs work
* currently seems to load neutrals twice. check this and fix
* separate beam from project_image(). get beam angle coords in plot function itself
* add tab for plotting beam grid, los, beam centerline
"""


"""Taken from http://nbviewer.jupyter.org/gist/tillahoffmann/f844bce2ec264c1c8cb5 or
https://gist.github.com/tillahoffmann/f844bce2ec264c1c8cb5
"""
class gaussian_kde(object):
    """Representation of a kernel-density estimate using Gaussian kernels.

    Kernel density estimation is a way to estimate the probability density
    function (PDF) of a random variable in a non-parametric way.
    `gaussian_kde` works for both uni-variate and multi-variate data.   It
    includes automatic bandwidth determination.  The estimation works best for
    a unimodal distribution; bimodal or multi-modal distributions tend to be
    oversmoothed.

    Parameters
    ----------
    dataset : array_like
        Datapoints to estimate from. In case of univariate data this is a 1-D
        array, otherwise a 2-D array with shape (# of dims, # of data).
    bw_method : str, scalar or callable, optional
        The method used to calculate the estimator bandwidth.  This can be
        'scott', 'silverman', a scalar constant or a callable.  If a scalar,
        this will be used directly as `kde.factor`.  If a callable, it should
        take a `gaussian_kde` instance as only parameter and return a scalar.
        If None (default), 'scott' is used.  See Notes for more details.
    weights : array_like, shape (n, ), optional, default: None
        An array of weights, of the same shape as `x`.  Each value in `x`
        only contributes its associated weight towards the bin count
        (instead of 1).

    Attributes
    ----------
    dataset : ndarray
        The dataset with which `gaussian_kde` was initialized.
    d : int
        Number of dimensions.
    n : int
        Number of datapoints.
    neff : float
        Effective sample size using Kish's approximation.
    factor : float
        The bandwidth factor, obtained from `kde.covariance_factor`, with which
        the covariance matrix is multiplied.
    covariance : ndarray
        The covariance matrix of `dataset`, scaled by the calculated bandwidth
        (`kde.factor`).
    inv_cov : ndarray
        The inverse of `covariance`.

    Methods
    -------
    kde.evaluate(points) : ndarray
        Evaluate the estimated pdf on a provided set of points.
    kde(points) : ndarray
        Same as kde.evaluate(points)
    kde.pdf(points) : ndarray
        Alias for ``kde.evaluate(points)``.
    kde.set_bandwidth(bw_method='scott') : None
        Computes the bandwidth, i.e. the coefficient that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        .. versionadded:: 0.11.0
    kde.covariance_factor : float
        Computes the coefficient (`kde.factor`) that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        The default is `scotts_factor`.  A subclass can overwrite this method
        to provide a different method, or set it through a call to
        `kde.set_bandwidth`.

    Notes
    -----
    Bandwidth selection strongly influences the estimate obtained from the KDE
    (much more so than the actual shape of the kernel).  Bandwidth selection
    can be done by a "rule of thumb", by cross-validation, by "plug-in
    methods" or by other means; see [3]_, [4]_ for reviews.  `gaussian_kde`
    uses a rule of thumb, the default is Scott's Rule.

    Scott's Rule [1]_, implemented as `scotts_factor`, is::

        n**(-1./(d+4)),

    with ``n`` the number of data points and ``d`` the number of dimensions.
    Silverman's Rule [2]_, implemented as `silverman_factor`, is::

        (n * (d + 2) / 4.)**(-1. / (d + 4)).

    Good general descriptions of kernel density estimation can be found in [1]_
    and [2]_, the mathematics for this multi-dimensional implementation can be
    found in [1]_.

    References
    ----------
    .. [1] D.W. Scott, "Multivariate Density Estimation: Theory, Practice, and
           Visualization", John Wiley & Sons, New York, Chicester, 1992.
    .. [2] B.W. Silverman, "Density Estimation for Statistics and Data
           Analysis", Vol. 26, Monographs on Statistics and Applied Probability,
           Chapman and Hall, London, 1986.
    .. [3] B.A. Turlach, "Bandwidth Selection in Kernel Density Estimation: A
           Review", CORE and Institut de Statistique, Vol. 19, pp. 1-33, 1993.
    .. [4] D.M. Bashtannyk and R.J. Hyndman, "Bandwidth selection for kernel
           conditional density estimation", Computational Statistics & Data
           Analysis, Vol. 36, pp. 279-298, 2001.

    Examples
    --------
    Generate some random two-dimensional data:

    >>> from scipy import stats
    >>> def measure(n):
    >>>     "Measurement model, return two coupled measurements."
    >>>     m1 = np.random.normal(size=n)
    >>>     m2 = np.random.normal(scale=0.5, size=n)
    >>>     return m1+m2, m1-m2

    >>> m1, m2 = measure(2000)
    >>> xmin = m1.min()
    >>> xmax = m1.max()
    >>> ymin = m2.min()
    >>> ymax = m2.max()

    Perform a kernel density estimate on the data:

    >>> X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    >>> positions = np.vstack([X.ravel(), Y.ravel()])
    >>> values = np.vstack([m1, m2])
    >>> kernel = stats.gaussian_kde(values)
    >>> Z = np.reshape(kernel(positions).T, X.shape)

    Plot the results:

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
    ...           extent=[xmin, xmax, ymin, ymax])
    >>> ax.plot(m1, m2, 'k.', markersize=2)
    >>> ax.set_xlim([xmin, xmax])
    >>> ax.set_ylim([ymin, ymax])
    >>> plt.show()

    """
    def __init__(self, dataset, bw_method=None, weights=None):
        self.dataset = np.atleast_2d(dataset)
        if not self.dataset.size > 1:
            raise ValueError("`dataset` input should have multiple elements.")
        self.d, self.n = self.dataset.shape

        if weights is not None:
            self.weights = weights / np.sum(weights)
        else:
            self.weights = np.ones(self.n) / self.n

        # Compute the effective sample size
        # http://surveyanalysis.org/wiki/Design_Effects_and_Effective_Sample_Size#Kish.27s_approximate_formula_for_computing_effective_sample_size
        self.neff = 1.0 / np.sum(self.weights ** 2)

        self.set_bandwidth(bw_method=bw_method)

    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.

        Returns
        -------
        values : (# of points,)-array
            The values at each point.

        Raises
        ------
        ValueError : if the dimensionality of the input points is different than
                     the dimensionality of the KDE.

        """
        points = np.atleast_2d(points)

        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = np.reshape(points, (self.d, 1))
                m = 1
            else:
                msg = "points have dimension %s, dataset has dimension %s" % (d,
                    self.d)
                raise ValueError(msg)

        # compute the normalised residuals
        chi2 = cdist(points.T, self.dataset.T, 'mahalanobis', VI=self.inv_cov) ** 2
        # compute the pdf
        result = np.sum(np.exp(-.5 * chi2) * self.weights, axis=1) / self._norm_factor

        return result

    __call__ = evaluate

    def scotts_factor(self):
        return np.power(self.neff, -1./(self.d+4))

    def silverman_factor(self):
        return np.power(self.neff*(self.d+2.0)/4.0, -1./(self.d+4))

    #  Default method to calculate bandwidth, can be overwritten by subclass
    covariance_factor = scotts_factor

    def set_bandwidth(self, bw_method=None):
        """Compute the estimator bandwidth with given method.

        The new bandwidth calculated after a call to `set_bandwidth` is used
        for subsequent evaluations of the estimated density.

        Parameters
        ----------
        bw_method : str, scalar or callable, optional
            The method used to calculate the estimator bandwidth.  This can be
            'scott', 'silverman', a scalar constant or a callable.  If a
            scalar, this will be used directly as `kde.factor`.  If a callable,
            it should take a `gaussian_kde` instance as only parameter and
            return a scalar.  If None (default), nothing happens; the current
            `kde.covariance_factor` method is kept.

        Notes
        -----
        .. versionadded:: 0.11

        Examples
        --------
        >>> x1 = np.array([-7, -5, 1, 4, 5.])
        >>> kde = stats.gaussian_kde(x1)
        >>> xs = np.linspace(-10, 10, num=50)
        >>> y1 = kde(xs)
        >>> kde.set_bandwidth(bw_method='silverman')
        >>> y2 = kde(xs)
        >>> kde.set_bandwidth(bw_method=kde.factor / 3.)
        >>> y3 = kde(xs)

        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111)
        >>> ax.plot(x1, np.ones(x1.shape) / (4. * x1.size), 'bo',
        ...         label='Data points (rescaled)')
        >>> ax.plot(xs, y1, label='Scott (default)')
        >>> ax.plot(xs, y2, label='Silverman')
        >>> ax.plot(xs, y3, label='Const (1/3 * Silverman)')
        >>> ax.legend()
        >>> plt.show()

        """
        if bw_method is None:
            pass
        elif bw_method == 'scott':
            self.covariance_factor = self.scotts_factor
        elif bw_method == 'silverman':
            self.covariance_factor = self.silverman_factor
        elif np.isscalar(bw_method) and not isinstance(bw_method, string_types):
            self._bw_method = 'use constant'
            self.covariance_factor = lambda: bw_method
        elif callable(bw_method):
            self._bw_method = bw_method
            self.covariance_factor = lambda: self._bw_method(self)
        else:
            msg = "`bw_method` should be 'scott', 'silverman', a scalar " \
                  "or a callable."
            raise ValueError(msg)

        self._compute_covariance()

    def _compute_covariance(self):
        """Computes the covariance matrix for each Gaussian kernel using
        covariance_factor().
        """
        self.factor = self.covariance_factor()
        # Cache covariance and inverse covariance of the data
        if not hasattr(self, '_data_inv_cov'):
            # Compute the mean and residuals
            _mean = np.sum(self.weights * self.dataset, axis=1)
            _residual = (self.dataset - _mean[:, None])
            # Compute the biased covariance
            self._data_covariance = np.atleast_2d(np.dot(_residual * self.weights, _residual.T))
            # Correct for bias (http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_covariance)
            self._data_covariance /= (1 - np.sum(self.weights ** 2))
            self._data_inv_cov = np.linalg.inv(self._data_covariance)

        self.covariance = self._data_covariance * self.factor**2
        self.inv_cov = self._data_inv_cov / self.factor**2
        self._norm_factor = np.sqrt(np.linalg.det(2*np.pi*self.covariance)) #* self.n

def to_angle_space(a, xhat, yhat, zhat):
    """Convert from distance-space to angle-space.

    Let the normalized, reference vector be :math:`\\hat{z}`, the common position to be point :math:`\\vec{O}`, and a point in
    question be :math:`\\vec{P}`. Let :math:`\\vec{a}=\\vec{P}-\\vec{O}`. Then the two angles giving the deviation of
    :math:`\\vec{P}` from :math:`\\hat{z}` are given by:

    .. math:: \\theta_1 = sign(\\vec{a} \\cdot \\hat{x}) \\cos^{-1} \\left(\\frac{(\\vec{a} - [\\vec{a} \\cdot \\hat{y}] \
    \\hat{y}) \\cdot \\hat{z}}{||(\\vec{a} - [\\vec{a} \\cdot \\hat{y}]\\hat{y}||}\\right)

    .. math:: \\theta_2 = sign(\\vec{a} \\cdot \\hat{y}) \\cos^{-1} \\left( \\frac{(\\vec{a} - [\\vec{a} \\cdot \\hat{x}] \
    \\hat{x}) \\cdot\\hat{z}}{||(\\vec{a} - [\\vec{a} \\cdot \\hat{x}]\\hat{x}||} \\right)

    Parameters
    ----------
    a : array, (nchan, 3)
        Vectors to be converted to angles. Emanating from single location ('lens')

    xhat : array, (3)
        Unit vector perp to yhat and zhat

    yhat : array, (3)
        Unit vector perp to xhat and zhat

    zhat : array, (3)
        Unit vector defining zero-angle position. Points from lens to some other point.

    Returns
    -------
    angles : array (nchan, 2)
        Angle that 'a' deviates from zhat in the xhat direction and in the yhat direction
    """
    nchan = a.shape[0]

    # Find angle in xhat, zhat plane
    a_dot_xhat = np.dot(a, xhat)
    a_dot_yhat = np.dot(a, yhat)
    a_dot_yhat_yhat = a_dot_yhat.reshape(nchan, 1) * yhat
    a_minus_a_dot_yhat_yhat = a - a_dot_yhat_yhat
    a_minus_a_dot_yhat_yhat_mag = np.linalg.norm(a_minus_a_dot_yhat_yhat, axis=1)
    arg1 = np.dot(a_minus_a_dot_yhat_yhat, zhat) / a_minus_a_dot_yhat_yhat_mag
    arg1[arg1 > 1.] = 1.                              # remove floating point errors
    ang1 = np.sign(a_dot_xhat) * np.arccos(arg1)      # (nchan)

    # Find angle in yhat, zhat plane
    a_dot_xhat = np.sum(a * xhat.reshape(1, 3), axis=1)
    a_dot_xhat_xhat = a_dot_xhat.reshape(nchan, 1) * xhat
    a_minus_a_dot_xhat_xhat = a - a_dot_xhat_xhat
    a_minus_a_dot_xhat_xhat_mag = np.linalg.norm(a_minus_a_dot_xhat_xhat, axis=1)
    arg2 = np.dot(a_minus_a_dot_xhat_xhat, zhat) / a_minus_a_dot_xhat_xhat_mag
    arg2[arg2 > 1.] = 1.                             # remove floating point errors
    ang2 = np.sign(a_dot_yhat) * np.arccos(arg2)     # (nchan)

    return np.array([ang1, ang2]).T     # (nchan, 2)

def project_image(axis=None,
                  lens=None,
                  data=None,
                  beam_pt=None,
                  beam_axis=None):
    """Given several lines of sight and an intensity per LOS, project an image on a plane perpendicular
    to the average LOS axis. USES ANGLES FOUND IN A GENERAL WAY

    Let the normalized, average LOS be :math:`\\hat{z}`, the lens position to be point :math:`\\vec{O}`, and a point in
    question be :math:`\\vec{P}`. Let :math:`\\vec{a}=\\vec{P}-\\vec{O}`.

    Parameters
    ----------
    axis : array (nchan, 3)
        Normalized axis vectors defining LOS

    lens : array (3)
        Common location for all LOS (aperture or lens)

    data : array (nchan)
        Data to be projected

    beam_pt : array (3)
        A point on the beam centerline

    beam_axis : array (3)
        Vector along beam centerline

    Returns
    -------
    x1 : array (100)
        Relative coordinates for grid_data (angle (deg) perpendicular to average LOS axis)

    x2 : array (101)
        Second set of relative coordinates for grid_data (angle (deg) perpendicular to average LOS axis)

    grid_data : array (100, 101)
        Data interpolated onto a uniform grid on a plane perpendicular to the average LOS axis

    beam_pt1 : array (2)
        Beam coordinates in 2-angle space for 1st beam point

    beam_pt2 : array (2)
        Beam coordinates in 2-angle space for 2nd beam point

    Todo
    ----
    * Choose x or yhat to consistantly be the one that is more parallel to beam axis. Could do by crossing
      beam_axis instead of any_vec.
    """
    # Average LOS axis
    zhat = axis.mean(0)
    zhat = zhat / np.linalg.norm(zhat)  # (3)

#    # Find any vector perp to zhat (by crossing w/ any non-colinear vector) to define the plane
#    any_vec = np.array([zhat[0] + 5., zhat[1], zhat[2]])   # 5. is arbitrary
#    xhat = np.cross(zhat, any_vec)
#    xhat = xhat / np.linalg.norm(xhat)
#
#    # Find second plane vector perp to first
#    yhat = np.cross(zhat, xhat)  # (3)

    # Get unit vector perp to beam
    yhat = np.cross(zhat, beam_axis) / np.linalg.norm(beam_axis)

    # Get unit vector along beam
    xhat = np.cross(yhat, zhat)

    # Get angles along xhat and yhat for LOS
    angs = to_angle_space(axis, xhat, yhat, zhat)    # (nchan, 2)

    # Interpolate data onto uniform grid of angles
    n1d = 100    # no. of grid points in each direction
    x1 = np.linspace(angs[:, 0].min(), angs[:, 0].max(), num = n1d)
    x2 = np.linspace(angs[:, 1].min(), angs[:, 1].max(), num = n1d + 1)
    x1_grid, x2_grid = np.meshgrid(x1, x2, indexing='ij')
    grid_data = interpolate.griddata(np.array([angs[:, 0], angs[:, 1]]).T, data, (x1_grid, x2_grid), fill_value=0.)

    # Pick two points on the beam centerline
    beam_vecs = np.zeros((2, 3))
    beam_vecs[0, :] = beam_pt - lens
    beam_vecs[1, :] = beam_pt + beam_axis * 10. - lens

    # Get angles along xhat and yhat for beam points
    beam_angs = to_angle_space(beam_vecs, xhat, yhat, zhat)   # (2, 2) = (2-pts, 2-space)

    # Define beam centerline in angle coordinates and move beam points to edge of angle grid
    beam_axis = np.squeeze(np.diff(beam_angs, axis=0))
    if beam_axis[0] != 0.:
        t1 = (x1.min() - beam_angs[0, 0]) / beam_axis[0]
        t2 = (x1.max() - beam_angs[0, 0]) / beam_axis[0]
    else:
        t1 = -np.inf
        t2 = np.inf
    if beam_axis[1] != 0.:
        t1_b = (x2.min() - beam_angs[0, 1]) / beam_axis[1]
        t2_b = (x2.max() - beam_angs[0, 1]) / beam_axis[1]
        t1 = np.max([t1, t1_b])
        t2 = np.min([t2, t2_b])
    beam_angs[0, :] = [beam_angs[0, 0] + beam_axis[0] * t1, beam_angs[0, 1] + beam_axis[1] * t1]
    beam_angs[1, :] = [beam_angs[1, 0] + beam_axis[0] * t2, beam_angs[1, 1] + beam_axis[1] * t2]

    return x1, x2, grid_data, beam_angs[0, :], beam_angs[1, :]

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
        Total number of spectral channels (lines of sight)

    lens_loc : 2D array
        Cartesian coords of all lenses in machine coords, (nchan, 3)

    Returns
    -------
    uniq_lens_indices : list
        Indeces to locate spectra for each unique lens location

    nlenses : int
        Number of unique len locations

    """
    nchan = lens_loc.shape[0]
    chan = np.arange(nchan)
    lens_list = [tuple(lens_loc[i,:]) for i in range(nchan)]
    lens_set = set(lens_list)

    uniq_lens_indices = []
    for lens in lens_set:
        wl = np.array([l == lens for l in lens_list])
        lens_ind = chan[wl]
        uniq_lens_indices.append(lens_ind)

    return uniq_lens_indices, len(uniq_lens_indices)


class Spectra:
    """ Spectra object that contains plot methods and parameters"""
    def __init__(self, nml):
        result_dir = nml["result_dir"]
        runid = nml["runid"]
        spec_file = os.path.join(result_dir, runid+'_spectra.h5')
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
#            self.projection_dist = tk.StringVar(value = 100.)

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

            # Get all LOS vectors for this lens
            lens_axis = self.lens_axis[ch, :]           # (this_nchan, 3), all LOS axes for this lens
            lens_loc = self.lens_loc[ch[0], :]          # (3), same for all in ch

            # Project all LOS data onto 2D grid perpendicular to average LOS, a distance projection_dist from the lens
            x1, x2, grid_spec, beam_pt1, beam_pt2 = project_image(axis=lens_axis,
                                                                  lens=lens_loc,
                                                                  data=spec,
                                                                  beam_pt=self.beam_src,
                                                                  beam_axis=self.beam_axis)

            # Plot contour
            c = ax.contourf(np.degrees(x1), np.degrees(x2), grid_spec.T, 50)
            cb = fig.colorbar(c)
            cb.ax.set_ylabel('[$Ph\ /\ (s\ sr\ m^2)$]')
            ax.set_title('Intensity\nLens at [{:4.0f},{:4.0f},{:4.0f}]'.format(lens_loc[0], lens_loc[1], lens_loc[2]))
            ax.set_xlabel('X1 [deg.]')
            ax.set_ylabel('X2 [deg.]')

            # Overplot beam centerline
#            beam_pt1 = np.degrees(beam_pt1)
#            beam_pt2 = np.degrees(beam_pt2)
#            ax.plot([beam_pt1[0], beam_pt2[0]], [beam_pt1[1], beam_pt2[1]], color = 'magenta')
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
        result_dir = nml["result_dir"]
        runid = nml["runid"]
        npa_file = os.path.join(result_dir, runid + '_npa.h5')
        wght_file = os.path.join(result_dir, runid + '_npa_weights.h5')
        neut_file = os.path.join(result_dir, runid + '_neutrals.h5')
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
        result_dir = nml["result_dir"]
        runid = nml["runid"]
        npa_wght_file = os.path.join(result_dir,runid+'_npa_weights.h5')
        fida_wght_file = os.path.join(result_dir,runid+'_fida_weights.h5')

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
        result_dir = nml["result_dir"]
        runid = nml["runid"]
        neut_file = os.path.join(result_dir,runid+'_neutrals.h5')
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
                        use_histogram = True

                        if use_histogram is True:
                            ###############################################
                            # MAINTAIN ORIGINAL NUMBER OF POINTS
                            ###############################################

                            # Use machine coords and they're not the same as beam coords (so must rebin)
                            ax.set_xlabel('X [cm]')
                            ax.set_ylabel('Y [cm]')

                            # Need to bin data onto mach regular grid before taking projections
                            fdens_hist = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(),
                                                        bins = (self.nx, self.ny), weights=self.fdens.flatten())
                            fdens = fdens_hist[0]

                            # Histogram returns edges of shape (nx+1). Convert to centers
                            xedges = fdens_hist[1]
                            yedges = fdens_hist[2]
                            dx = xedges[1] - xedges[0]
                            dy = yedges[1] - yedges[0]
                            x = xedges[0:-1] + dx / 2.
                            y = yedges[0:-1] + dy / 2.

                            hdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(),
                                                   bins = (xedges, yedges), weights=self.hdens.flatten())[0]
                            tdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(),
                                                   bins = (xedges, yedges), weights=self.tdens.flatten())[0]
                            halodens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(),
                                                      bins = (xedges, yedges), weights=self.halodens.flatten())[0]

                            # histogram2d sums weights, need mean
                            fdens = fdens / self.nz
                            hdens = hdens / self.nz
                            tdens = tdens / self.nz
                            halodens = halodens / self.nz

#                            print('original dx {}'.format((self.x_grid[1, 0, 0] - self.x_grid[0, 0, 0]) * (self.y_grid[0, 1, 0] - self.y_grid[0, 0, 0])))
#                            print('new dx {}'.format(dx * dy))

                            ###############################################
                            # MAINTAIN ORIGINAL RESOLUTION
                            ###############################################

#                            # Use machine coords and they're not the same as beam coords (so must rebin)
#                            ax.set_xlabel('X [cm]')
#                            ax.set_ylabel('Y [cm]')
#
#                            dx = np.abs(self.x_grid[1, 0, 0] - self.x_grid[0, 0, 0])
#                            dy = np.abs(self.y_grid[0, 1, 0] - self.y_grid[0, 0, 0])
#                            xmin, xmax = self.x_grid.min(), self.x_grid.max()
#                            ymin, ymax = self.y_grid.min(), self.y_grid.max()
#                            nx = np.floor((xmax - xmin) / dx) + 1
#                            ny = np.floor((ymax - ymin) / dy) + 1
#                            x = np.linspace(xmin, xmax, num=nx)
#                            y = np.linspace(ymin, ymax, num=ny)
#                            xedges = x - dx / 2.
#                            yedges = y - dx / 2.
#                            xedges = np.append(xedges, xedges[-1] + dx)
#                            yedges = np.append(yedges, yedges[-1] + dy)
#
#                            # Need to bin data onto mach regular grid before taking projections
#                            fdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(),
#                                                        bins = (xedges, yedges), weights=self.fdens.flatten())[0]
#                            hdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(),
#                                                   bins = (xedges, yedges), weights=self.hdens.flatten())[0]
#                            tdens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(),
#                                                   bins = (xedges, yedges), weights=self.tdens.flatten())[0]
#                            halodens = np.histogram2d(self.x_grid.flatten(), self.y_grid.flatten(),
#                                                      bins = (xedges, yedges), weights=self.halodens.flatten())[0]
#
#                            # histogram2d sums weights, need mean
#                            fdens = fdens / self.nz
#                            hdens = hdens / self.nz
#                            tdens = tdens / self.nz
#                            halodens = halodens / self.nz
                        else:
                            # Use KDE
                            x = np.linspace(self.x_grid.min(), self.x_grid.max(), self.nx)
                            y = np.linspace(self.y_grid.min(), self.y_grid.max(), self.ny)
                            xx, yy = np.meshgrid(x, y)
                            fdens = np.zeros((3, np.max(self.fdens.shape)))  # make (ndims, ndata)
                            fdens[0, 0:self.nx] = self.fdens[:, 0, 0]
                            fdens[1, 0:self.ny] = self.fdens[0, :, 0]
                            fdens[2, 0:self.nz] = self.fdens[0, 0, :]
                            pdf = gaussian_kde(fdens, weights=fdens)
                            fdens = pdf((np.ravel(xx), np.ravel(yy)))
                            fdens = np.reshape(fdens, xx.shape)
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

#            ttk.Label(self.imaging_frame, text = 'Projection Distance (cm)').pack()
#            ttk.Entry(self.imaging_frame, textvariable = self.spec.projection_dist, state = tk.NORMAL, width = 10).pack()
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
        self.namelistfile = askopenfilename(filetypes=[('Namelist Files','*.dat')])
        self.nml = self.read_nml(self.namelistfile)
        self.spec = Spectra(self.nml)
        self.npa = NPA(self.nml)
        self.neut = Neutrals(self.nml)
        self.wght = Weights(self.nml)

if __name__ == '__main__':
    root = tk.Tk()
    Viewer(root)
    root.mainloop()

