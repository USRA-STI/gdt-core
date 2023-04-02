# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.colors import colorConverter
from matplotlib.patches import Polygon
from scipy.spatial.transform import Rotation

from .defaults import *
from gdt.core.data_primitives import EnergyBins
from gdt.core.coords import SpacecraftFrame

__all__ = ['earth_line', 'earth_points', 'effective_area', 'errorband', 
           'galactic_plane', 'histo', 'histo_errorbars', 'histo_filled', 
           'lightcurve_background', 'response_matrix', 'saa_polygon',
           'selection_line', 'selections',  'sky_annulus','sky_circle', 
           'sky_heatmap', 'sky_line', 'sky_point', 'sky_polygon',
           'spectrum_background']

# ---------- Lightcurve and Spectra ----------#
def histo(bins, ax, color='C0', edges_to_zero=False, **kwargs):
    """Plot a rate histogram either lightcurves or count spectra.
    
    Args:
        bins (:class:`~gdt.core.data_primitives.TimeBins` or \
              :class:`~gdt.core.data_primitives.EnergyBins`):
            The lightcurve or count spectrum histograms
        ax (:class:`matplotlib.axes`): The axis on which to plot.
        color (str, optional): The color of the histogram. Default is 'C0'
        edges_to_zero (bool, optional):
            If True, then the farthest edges of the histogram will drop to zero.
            Default is True.
        **kwargs: Other plotting options
    
    Returns:
        (list of matplotlib.lines.Line2D)
    """
    bin_segs = bins.contiguous_bins()
    refs = []
    for seg in bin_segs:
        edges = np.concatenate(
            ([seg.lo_edges[0]], seg.lo_edges, [seg.hi_edges[-1]]))
        if isinstance(seg, EnergyBins):
            rates_attr = seg.rates_per_kev
        else:
            rates_attr = seg.rates
        if edges_to_zero:
            rates = np.concatenate(([0.0], rates_attr, [0.0]))
        else:
            rates = np.concatenate(
                ([rates_attr[0]], rates_attr, [rates_attr[-1]]))

        p = ax.step(edges, rates, where='post', color=color, **kwargs)
        refs.append(p)
    return refs


def histo_errorbars(bins, ax, color='C0', **kwargs):
    """Plot errorbars for lightcurves or count spectra.
    
    Args:
        bins (:class:`~gdt.core.data_primitives.TimeBins` or \
              :class:`~gdt.core.data_primitives.EnergyBins`):
            The lightcurve or count spectrum histograms
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The color of the errorbars. Default is 'C0'
        **kwargs: Other plotting options
    
    Returns:
        (list of matplotlib.container.ErrorbarContainer)
    """
    bin_segs = bins.contiguous_bins()
    refs = []
    for seg in bin_segs:
        if isinstance(seg, EnergyBins):
            rates_attr = seg.rates_per_kev
            uncert_attr = seg.rate_uncertainty_per_kev
        else:
            rates_attr = seg.rates
            uncert_attr = seg.rate_uncertainty           
        p = ax.errorbar(seg.centroids, rates_attr, uncert_attr,
                        capsize=0, fmt='none', color=color, **kwargs)
        refs.append(p)
    return refs


def histo_filled(bins, ax, color=DATA_SELECTED_COLOR,
                 fill_alpha=DATA_SELECTED_ALPHA, **kwargs):
    """Plot a filled histogram.
    
    Args:
        bins (:class:`~gdt.core.data_primitives.TimeBins` or \
              :class:`~gdt.core.data_primitives.EnergyBins`):
            The lightcurve or count spectrum histograms
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The color of the filled histogram
        fill_alpha (float, optional): The alpha of the fill
        **kwargs: Other plotting options
    
    Returns:
        (list of matplotlib.lines.Line2D and \
         matplotlib.collections.PolyCollection)
    """
    refs = histo(bins, ax, color=color, zorder=3, **kwargs)
    if isinstance(bins, EnergyBins):
        rates_attr = bins.rates_per_kev
    else:
        rates_attr = bins.rates        
    zeros = np.zeros(bins.size + 1)
    rates = np.append(rates_attr, rates_attr[-1])
    edges = np.append(bins.lo_edges, bins.hi_edges[-1])
    b1 = ax.plot((edges[0], edges[0]), (zeros[0], rates[0]), color=color,
                 zorder=2)
    b2 = ax.plot((edges[-1], edges[-1]), (zeros[-1], rates[-1]), color=color,
                 zorder=2)
    f = ax.fill_between(edges, zeros, rates, color=color, step='post',
                        alpha=fill_alpha, zorder=4, **kwargs)
    
    refs.extend([b1, b2, f])
    return refs


def selection_line(xpos, ax, **kwargs):
    """Plot a selection line.
    
    Args:
        xpos (float): The position of the selection line
        ax (:class:`matplotlib.axes`): The axis on which to plot
        **kwargs: Other plotting options
    
    Returns:
        (matplotlib.lines.Line2D)
    """
    ylim = ax.get_ylim()
    ref = ax.plot([xpos, xpos], ylim, **kwargs)
    return ref


def selections(bounds, ax, **kwargs):
    """Plot selection bounds.
    
    Args:
        bounds (list of tuples): List of selection bounds
        ax (:class:`matplotlib.axes`): The axis on which to plot
        **kwargs: Other plotting options
    
    Returns:
        (2-tuple of matplotlib.lines.Line2D)
    """
    refs1 = []
    refs2 = []
    for bound in bounds:
        p = selection_line(bound[0], ax, **kwargs)
        refs1.append(p[0])
        p = selection_line(bound[1], ax, **kwargs)
        refs2.append(p[0])

    return (refs1, refs2)


def errorband(x, y_upper, y_lower, ax, **kwargs):
    """Plot an error band.
    
    Args:
        x (np.array): The x values
        y_upper (np.array): The upper y values of the error band
        y_lower (np.array): The lower y values of the error band
        ax (:class:`matplotlib.axes`): The axis on which to plot
        **kwargs: Other plotting options
    
    Returns:
        (list of matplotlib.collections.PolyCollection)    
    """
    refs = ax.fill_between(x, y_upper.squeeze(), y_lower.squeeze(), **kwargs)
    return refs


def lightcurve_background(backrates, ax, cent_color=None, err_color=None,
                          cent_alpha=None, err_alpha=None, **kwargs):
    """Plot a lightcurve background model with an error band.
    
    Args:
        backrates (:class:`~gdt.background.primitives.BackgroundRates`):
            The background rates object integrated over energy. If there is more
            than one remaining energy channel, the background will be integrated
            over the remaining energy channels.
        ax (:class:`matplotlib.axes`): The axis on which to plot
        cent_color (str): Color of the centroid line
        err_color (str): Color of the errorband
        cent_alpha (float): Alpha of the centroid line
        err_alpha (float): Alpha of the errorband
        **kwargs: Other plotting options
    
    Returns:
        (list of matplotlib.lines.Line2D and \
         matplotlib.collections.PolyCollection)   
    """
    if backrates.num_chans > 1:
        backrates = backrates.integrate_energy()
    times = backrates.time_centroids
    rates = backrates.rates
    uncert = backrates.rate_uncertainty
    p2 = errorband(times, rates + uncert, rates - uncert, ax, alpha=err_alpha,
                   color=err_color, linestyle='-', **kwargs)
    p1 = ax.plot(times, rates, color=cent_color, alpha=cent_alpha,
                 **kwargs)
    refs = [p1, p2]
    return refs


def spectrum_background(backspec, ax, cent_color=None, err_color=None,
                        cent_alpha=None, err_alpha=None, **kwargs):
    """Plot a count spectrum background model with an error band.
    
    Args:
        backspec (:class:`~gdt.background.primitives.BackgroundSpectrum`):
            The background rates object integrated over energy. If there is more
            than one remaining energy channel, the background will be integrated
            over the remaining energy channels.
        ax (:class:`matplotlib.axes`): The axis on which to plot
        cent_color (str): Color of the centroid line
        err_color (str): Color of the errorband
        cent_alpha (float): Alpha of the centroid line
        err_alpha (fl): Alpha of the errorband
        **kwargs: Other plotting options
    
    Returns:
        (list of matplotlib.lines.Line2D and \
         matplotlib.collections.PolyCollection)   
    """
    rates = backspec.rates_per_kev
    uncert = backspec.rate_uncertainty_per_kev
    edges = np.append(backspec.lo_edges, backspec.hi_edges[-1])
    # plot the centroid of the model
    p1 = ax.step(edges, np.append(rates, rates[-1]), where='post',
                 color=cent_color, alpha=cent_alpha, **kwargs)

    # construct the stepped errorband to fill between
    energies = np.array(
        (backspec.lo_edges, backspec.hi_edges)).T.flatten()
    upper = np.array((rates + uncert, rates + uncert)).T.flatten()
    lower = np.array((rates - uncert, rates - uncert)).T.flatten()
    
    # some change in matplotlib around v3.3.3 necessitates this next line.
    # basically, plotting a fill_between when neighboring x-values are identical
    # is not handled (crashes), whereas in older versions, it was handled
    # successfully and as expected.  this line should be superfluous because it
    # tells fill_between to fill between every single y pairs, which is what it
    # should do by default, so ¯\_(ツ)_/¯   
    where = np.ones(energies.size, dtype=bool)
    p2 = errorband(energies, upper, lower, ax, color=err_color,
                   alpha=err_alpha, where=where, **kwargs)
    refs = [p1, p2]
    return refs


# ---------- DRM ----------#
def response_matrix(phot_energies, chan_energies, matrix, ax, cmap='Greens',
                    num_contours=100, norm=None, **kwargs):
    """Make a filled contour plot of a response matrix.
    
    Args:
        phot_energies (np.array): The incident photon energy bin centroids
        chan_energies (np.array): The recorded energy channel centroids
        matrix (np.array): The effective area matrix corresponding to the photon bin and 
                           energy channels
        ax (matplotlib.axes): The axis on which to plot
        cmap (str, optional): The color map to use. Default is 'Greens'
        num_contours (int, optional): The number of contours to draw. These will 
                                      be equally spaced in log-space. Default 
                                      is 100
        norm (matplotlib.colors.Normalize or similar, optional):
            The normalization used to scale the colormapping to the heatmap 
            values. This can be initialized by Normalize, LogNorm, SymLogNorm, 
            PowerNorm, or some custom normalization.
        **kwargs: Other keyword arguments to be passed to 
                  matplotlib.pyplot.contourf
    
    Returns:
        (matplotlib.collections.QuadMesh)
    """
    mask = (matrix > 0.0)
    levels = np.geomspace(matrix[mask].min(), matrix.max(), num_contours)
    image = ax.contourf(phot_energies, chan_energies, matrix, levels=levels,
                        cmap=cmap, norm=norm)
    return image


def effective_area(bins, ax, color='C0', orientation='vertical', **kwargs):
    """Plot a histogram of the effective area of an instrument response.
    
    Args:
        bins: (:class:`~gdt.core.data_primitives.Bins`): The histogram of 
                                                       effective area
        ax (matplotlib.axes): The axis on which to plot
        **kwargs (optional): Other plotting options
    
    Returns:
        (list of matplotlib.lines.line2D)
    """
    edges = np.concatenate(
        ([bins.lo_edges[0]], bins.lo_edges, [bins.lo_edges[-1]]))
    counts = np.concatenate(([0.0], bins.counts, [0.0]))
    if orientation == 'horizontal':
        p = ax.step(counts, edges, where='post', color=color, **kwargs)
    else:
        p = ax.step(edges, counts, where='post', color=color, **kwargs)
    return [p]


# ---------- Earth and Orbital ----------#
def saa_polygon(lat_saa, lon_saa, proj, color='darkred', alpha=0.4, **kwargs):
    """Plot the SAA polygon on the Earth.
    
    Args:
        lat_saa (np.array): Array of latitude points
        lon_saa (np.array): Array of longitude points
        proj (Cartopy Projection): The Cartopy projection
        color (str, optional): The color of the polygon
        alpha (float, optional): The alpha opacity of the interior of the polygon
        kwargs (optional): Other plotting keywords
    
    Returns:
        (matplotlib.patches.Polygon)
    """
    edge = colorConverter.to_rgba(color, alpha=1.0)
    face = colorConverter.to_rgba(color, alpha=alpha)

    # plot the polygon
    xy = list(zip(lon_saa, lat_saa))
    poly = Polygon(xy, edgecolor=edge, facecolor=face, transform=proj, **kwargs)
    return poly


def earth_line(lat, lon, proj, color='C0', alpha=0.4, **kwargs):
    """Plot a line on the Earth (e.g. orbit).
    
    Args:
        lat (np.array): Array of latitudes
        lon (np.array): Array of longitudes
        proj (GeoAxesSubplot): The Cartopy projection 
        color (str, optional): The color of the line
        alpha (float, optional): The alpha opacity of the line
        kwargs (optional): Other plotting keywords
    
    Returns:
        (list of matplotlib.lines.Line2D)
    """
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    lon[(lon > 180.0)] -= 360.0
    path = np.vstack((lon, lat))
    isplit = np.nonzero(np.abs(np.diff(path[0])) > 5.0)[0]
    segments = np.split(path, isplit + 1, axis=1)

    refs = []
    for segment in segments:
        refs.append(proj.plot(segment[0], segment[1], color=color, alpha=alpha,
                              **kwargs))
    return refs


def earth_points(lat, lon, proj, color='C0', alpha=1.0, size=10, **kwargs):
    """Plot a point or points on the Earth.
    
    Args:
        lat (np.array): Array of latitudes
        lon (np.array): Array of longitudes
        proj (GeoAxesSubplot): The Cartopy projection 
        color (str, optional): The color of the point(s)
        alpha (float, optional): The alpha opacity of the point(s)
        size (float, optional): The size of the point(s). Default is 10.
        kwargs (optional): Other plotting keywords.
    
    Returns:
        (list of matplotlib.collections.PathCollection)
    """
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    lon[(lon > 180.0)] -= 360.0
    ref = proj.scatter(lon, lat, color=color, alpha=alpha, s=size, **kwargs)
    return [ref]


# ---------- Sky and Fermi Inertial Coordinates ----------#

def circle(center_x, center_y, radius, num_points=100):
    """Compute the points of a circle on a sphere as defined by the center of
    the circle and the angular radius.
    
    Args:
        center_x (float): The azimuthal center of the circle, in radians.
        center_y (float): The polar center of the circle, in radians.
        radius (float): The angular radius of the circle, in radians.
        num_points (int, optional): The number of points on the circle. 
                                    Default is 100.

    Returns:
        (list of np.array, list of np.array)
    """
    ra = 2.0 * np.pi * np.arange(num_points) / float(num_points)
    dec = np.array([np.pi/2.0 - radius] * ra.size)
    
    # rotation vectors for the z and y axes to the circle center
    rot1 = Rotation.from_euler('z', -center_x).as_matrix()
    rot2 = Rotation.from_euler('y', -center_y).as_matrix()
    
    # uses the rotation vectors to rotate the points on the circle to the
    # frame defined by the center of the circle
    coord = SkyCoord(ra, dec, unit='rad')
    xyz = SkyCoord(*np.dot(np.dot(rot2, rot1).squeeze().T, coord.cartesian.xyz),
                   representation_type='cartesian', frame='icrs')
    ra = xyz.spherical.lon.to('rad').value
    dec = xyz.spherical.lat.to('rad').value
    
    mask = ra > 2.0*np.pi
    ra[mask] = ra[mask] - 2.0*np.pi
    return ra, dec
    

def split_on_meridian(phi, theta):
    """Split a set of phi, theta points that span the meridian into contiguous
    segments.
    
    Args:
        phi (np.array): The azimuthal coordinate, in radians
        theta (np.array): The polar coordinate, in radians

    Returns:
        (list of np.array, list of np.array)
    """
    # if difference between consecutive phi points is > pi, then the
    # meridian splits the line/polygon
    isplit = np.nonzero(np.abs(np.diff(phi, prepend=phi[-1])) > np.pi)[0]
    
    # if only one split is detected, the line/polygon is at the pole
    if len(isplit) == 1:
        idx = phi.argsort()
        phi = phi[idx]
        theta = theta[idx]
                    
        # detect which pole we are wrapping around and assign an array of
        # thetas to be at that peole. special case (second case) is if we are 
        # precisely at the equator.
        if (np.pi/2.0 - theta.max()) < np.abs(-np.pi/2.0 - theta.min()):
            theta_edges = np.array([np.pi/2.0] * theta.size)
        elif ( (np.pi/2.0 - theta.max()) / \
               np.abs(-np.pi/2.0 - theta.min()) ) <=  (1.0 + 1e-3):
            north_size = int(np.ceil(theta.size/2))
            south_size = theta.size - north_size
            north_edges = np.array([np.pi/2.0] * north_size)
            south_edges = np.array([-np.pi/2.0] * south_size)
            theta_edges = np.concatenate((north_edges, south_edges))
        else:
            theta_edges = np.array([-np.pi/2.0] * theta.size)
        
        # assign a reversed array of phi, and set the start and end points to
        # the meridian
        phi_edges = phi[::-1]
        phi_edges[0] = 2.0*np.pi
        phi_edges[-1] = 0.0
        
        phi = [np.concatenate((phi, phi_edges))]
        theta = [np.concatenate((theta, theta_edges))]
    
    # if two splits, the line/polygon is not at the pole
    elif len(isplit) == 2:
        # split into two chunks
        phi = np.split(phi, isplit)
        theta = np.split(theta, isplit)
        if phi[0].size < phi[-1].size:
            phi = phi[1:]
            theta = theta[1:]
        else:
            phi = phi[:-1]
            theta = theta[:-1]
        
        for i in range(len(phi)):
            # detect which side of the meridian the chunk is on and assign an
            # array of phis to be on that side of the meridian
            if 2.0*np.pi - phi[i].max() < phi[i].min():
                phi_edges = np.array([2.0*np.pi] * phi[i].size)
            else:
                phi_edges = np.array([0.0] * phi[i].size)
            
            # assign a reversed array of theta
            theta_edges = theta[i][::-1]
            
            phi[i] = np.concatenate((phi[i], phi_edges))
            theta[i] = np.concatenate((theta[i], theta_edges))
    else:
        phi = [phi]
        theta = [theta]
    
    return (phi, theta)


def sky_point(x, y, ax, flipped=True, frame='equatorial', **kwargs):
    """Plot a point on the sky.
        
    Args:
        x (float): The azimuthal coordinate, in degrees
        y (float): The polar coordinate, in degrees
        ax (matplotlib.axes): Plot axes object
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', or 'spacecraft'.
                               Default is 'equatorial'
        **kwargs: Other plotting options
    
    Returns:
        (matplotlib.collections.PathCollection)
    """
    theta = np.array(np.deg2rad(y))
    phi = np.array(np.deg2rad(x - 180.0))
    if frame == 'spacecraft':
        flipped = False
        phi -= np.pi
        if phi < -np.pi:
            phi += 2 * np.pi
    elif frame == 'galactic':
        phi -= np.pi
        phi[phi < -np.pi] += 2 * np.pi   
    else:
        pass 

    if flipped:
        phi = -phi
    point = ax.scatter(phi, theta, **kwargs)
    return point


def sky_line(x, y, ax, flipped=True, frame='equatorial', **kwargs):
    """Plot a line on a sky map, wrapping at the meridian.
        
    Args:
        x (float): The azimuthal coordinates, in degrees
        y (float): The polar coordinates, in degrees
        ax (matplotlib.axes): Plot axes object
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', or 'spacecraft'.
                               Default is 'equatorial'
        **kwargs: Other plotting options
    
    Returns:
        (list of matplotlib.collections.PathCollection)
    """
    theta = np.deg2rad(y)
    phi = np.deg2rad(x - 180.0)
    if frame == 'spacecraft':
        flipped = False
        phi -= np.pi
        phi[phi < -np.pi] += 2 * np.pi
    elif frame == 'galactic':
        phi -= np.pi
        phi[phi < -np.pi] += 2 * np.pi   
    else:
        pass 

    if flipped:
        phi = -phi
    seg = np.vstack((phi, theta))

    # here is where we split the segments at the meridian
    isplit = np.nonzero(np.abs(np.diff(seg[0])) > np.pi)[0]
    subsegs = np.split(seg, isplit+1, axis=1)

    # plot each path segment
    segrefs = []
    for seg in subsegs:
        ref = ax.plot(seg[0], seg[1], **kwargs)
        segrefs.append(ref)

    return segrefs


def sky_circle(center_x, center_y, radius, ax, flipped=True, frame='equatorial',
               face_color=None, face_alpha=None, edge_color=None,
               edge_alpha=None, **kwargs):
    """Plot a circle on the sky.
    
    Args:
        center_x (float): The azimuthal center, in degrees
        center_y (float): The polar center, in degrees
        radius (float): The ROI radius in degrees
        ax (matplotlib.axes): Plot axes object
        flipped (bool, optional):  If True, the azimuthal axis is flipped, 
                                   following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', or 'spacecraft'.
                               Default is 'equatorial'
        face_color (str, optional): The color of the circle fill
        face_alpha (float, optional): The alpha of the circle fill
        edge_color (str, optional): The color of the circle edge
        edge_alpha (float, optional): The alpha of the circle edge
        **kwargs: Other plotting options

    Returns:
        (list of matplotlib.patches.Polygon)
    """
    zen_true = False
    theta = np.deg2rad(90.0 - center_y)
    phi = np.deg2rad(center_x)
    if frame == 'spacecraft':
        flipped = False
        phi -= np.pi
        if phi < -np.pi:
            phi += 2 * np.pi
    elif frame == 'galactic':
        phi -= np.pi
        if phi < -np.pi:
            phi += 2 * np.pi        
    else:
        pass

    rad = np.deg2rad(radius)

    # The native matplotlib functions don't cut and display the circle polygon
    # correctly on map projections
    
    num_pts = int(rad * 500.0)
    phi_pts, theta_pts = circle(phi, theta, rad, num_points=num_pts)
    phi_pts, theta_pts = split_on_meridian(phi_pts, theta_pts)
    
    # plot each polygon section
    edge = colorConverter.to_rgba(edge_color, alpha=edge_alpha)
    face = colorConverter.to_rgba(face_color, alpha=face_alpha)
    patches = []
    for i in range(len(phi_pts)):
        phi_pts[i] -= np.pi
        if flipped:
            phi_pts[i] *= -1.0
        pts = np.vstack((phi_pts[i], theta_pts[i])).T
        patch = ax.add_patch( plt.Polygon(pts, facecolor=face, edgecolor=edge, \
                             **kwargs))
        patches.append(patch)
    return patches


def sky_annulus(center_x, center_y, radius, width, ax, color='black',
                alpha=0.3,
                frame='equatorial', flipped=True, **kwargs):
    """Plot an annulus on the sky defined by its center, radius, and width.
    
    Args:
        center_x (float): The azimuthal center, in degrees
        center_y (float): The polar center, in degrees
        radius (float): The radius in degrees, defined as the angular distance 
                        from the center to the middle of the width of the 
                        annulus band
        width (float): The width of the annulus in degrees
        ax (matplotlib.axes): Plot axes object
        color (string, optional): The color of the annulus. Default is black.
        alpha (float, optional): The opacity of the annulus. Default is 0.3
        frame (str, optional): Either 'equatorial', 'galactic', or 'spacecraft'.
                               Default is 'equatorial'
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following the equatorial convention
        **kwargs: Other plotting options

    Returns:
        (list of matplotlib.patches.Polygon)
    """
    zen_true = False
    edge = colorConverter.to_rgba(color, alpha=1.0)
    face = colorConverter.to_rgba(color, alpha=alpha)

    inner_radius = np.deg2rad(radius - width / 2.0)
    outer_radius = np.deg2rad(radius + width / 2.0)
    center_theta = np.deg2rad(90.0 - center_y)
    center_phi = np.deg2rad(center_x)
    if frame == 'spacecraft':
        flipped = False
        center_phi = np.deg2rad(center_x - 180.0)
    elif frame == 'galactic':
        center_phi = np.deg2rad(center_x + 180.0)
    else:
        pass

    # get the plot points for the inner and outer circles
    num_pts = int(outer_radius * 500.0)
    phi_inner, theta_inner = circle(center_phi, center_theta, inner_radius,
                                    num_points=num_pts)
    phi_inner, theta_inner = split_on_meridian(phi_inner, theta_inner)
    inner = [np.vstack((p, t)).T for p, t in zip(phi_inner, theta_inner)]
    
    phi_outer, theta_outer = circle(center_phi, center_theta, outer_radius,
                                    num_points=num_pts)
    phi_outer, theta_outer = split_on_meridian(phi_outer, theta_outer)
    outer = [np.vstack((p, t)).T for p, t in zip(phi_outer, theta_outer)]

    x1 = []
    y1 = []
    x2 = []
    y2 = []
    polys = []
    # plot the inner circle
    for section in inner:
        section[:, 0] -= np.pi
        if flipped:
            section[:, 0] *= -1.0
        polys.append(plt.Polygon(section, ec=edge, fill=False, **kwargs))
        ax.add_patch(polys[-1])
        x1.extend(section[:, 0])
        y1.extend(section[:, 1])

    # plot the outer circle
    for section in outer:
        section[:, 0] -= np.pi
        if flipped:
            section[:, 0] *= -1.0
        polys.append(plt.Polygon(section, ec=edge, fill=False, **kwargs))
        ax.add_patch(polys[-1])
        x2.extend(section[:, 0])
        y2.extend(section[:, 1])

    # organize and plot the fill between the circles
    # organize and plot the fill between the circles
    x1.append(x1[0])
    y1.append(y1[0])
    x2.append(x2[0])
    y2.append(y2[0])
    x2 = x2[::-1]
    y2 = y2[::-1]
    xs = np.concatenate((x1, x2))
    ys = np.concatenate((y1, y2))
    f = ax.fill(np.ravel(xs), np.ravel(ys), facecolor=face, zorder=1000)
    
    polys.append(f)
    return polys


def sky_polygon(x, y, ax, face_color=None, edge_color=None, edge_alpha=1.0,
                face_alpha=0.3, flipped=True, frame='equatorial', **kwargs):
    """Plot single polygon on a sky map, wrapping at the meridian.
        
    Args:
        x (float): The azimuthal coordinates, in degrees
        y (float): The polar coordinates, in degrees
        ax (matplotlib.axes): Plot axes object
        face_color (str, optional): The color of the polygon fill
        face_alpha (float, optional): The alpha of the polygon fill
        edge_color (str, optional): The color of the polygon edge
        edge_alpha (float, optional): The alpha of the polygon edge
        frame (str, optional): Either 'equatorial', 'galactic', or 'spacecraft'.
                               Default is 'equatorial'
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following the equatorial convention
        **kwargs: Other plotting options
        
    Returns:
        (list of matplotlib.lines.Line2D and matplotlib.patches.Polygon)
    """
    refs = sky_line(x, y, ax, color=edge_color, alpha=edge_alpha, frame=frame,
                    flipped=flipped, **kwargs)

    theta = np.deg2rad(y)
    phi = np.deg2rad(x - 180.0)
    
    def split_func(phi, theta):
        # this is needed to determine when contour spans the full x range
        if (phi.min() == -np.pi) and (phi.max() == np.pi):
            if theta.mean() > 0.0:  # fill in positive dec
                theta = np.insert(theta, 0, np.pi / 2.0)
                theta = np.append(theta, np.pi / 2.0)
                phi = np.insert(phi, 0, -np.pi)
                phi = np.append(phi, np.pi)
            else:  # fill in negative dec
                theta = np.insert(theta, 0, -np.pi / 2.0)
                theta = np.append(theta, -np.pi / 2.0)
                phi = np.insert(phi, 0, np.pi)
                phi = np.append(phi, -np.pi)
        return (phi, theta)
    
    if frame == 'spacecraft':
        flipped = False
        phi1 = phi - np.pi
        theta1 = np.copy(theta)
        phi1, theta1 = split_func(phi1, theta1)
        f1 = ax.fill(phi1, theta1, color=face_color, alpha=face_alpha, **kwargs)

        phi2 = phi + np.pi
        theta2 = np.copy(theta)
        phi2, theta2 = split_func(phi2, theta2)
        f2 = ax.fill(phi2, theta2, color=face_color, alpha=face_alpha, **kwargs)
        refs.extend([f1, f2])

    elif frame == 'galactic':
        phi1 = phi - np.pi
        theta1 = np.copy(theta)
        if flipped:
            phi1 = -phi1
        phi1, theta1 = split_func(phi1, theta1)
        f1 = ax.fill(phi1, theta1, color=face_color, alpha=face_alpha, **kwargs)

        phi2 = phi + np.pi
        theta2 = np.copy(theta)
        if flipped:
            phi2 = -phi2
        phi2, theta2 = split_func(phi2, theta2)
        f2 = ax.fill(phi2, theta2, color=face_color, alpha=face_alpha, **kwargs)
        refs.extend([f1, f2])

    else:
        if flipped:
            phi = -phi
        f = ax.fill(phi, theta, color=face_color, alpha=face_alpha, **kwargs)
        refs.append(f)

    return refs


def galactic_plane(ax, flipped=True, frame='equatorial', zen=True, 
                   outer_color='dimgray', inner_color='black', line_alpha=0.5, 
                   center_alpha=0.75, **kwargs):
    """Plot the galactic plane on the sky.
        
    Args:
        ax (matplotlib.axes): Plot axes object
        outer_color (str, optional): The color of the outer line
        inner_color (str, optional): The color of the inner line
        line_alpha (float, optional): The alpha of the line
        center_alpha (float, optional): The alpha of the center
        frame (str or :class:`~gdt.core.coords.SpacecraftFrame`, optional): 
            If a string, then can either be 'equatorial' or 'galactic'.
            Otherwise, it is the spacecraft frame definition. Defaults is 
            'equatorial'.
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following the equatorial convention
        zen (bool, optional): Only used if the frame is a spacecraft frame. 
                              If set to true, plots the polar axis in zenith,
                              False plots in elevation. Default is True.
        **kwargs: Other plotting options
    
    Returns:
        (list of matplotlib.lines.Line2D and 
         matplotlib.collections.PathCollection)
    """
    x = np.arange(0, 360, dtype=float)
    y = np.zeros_like(x)
    
    if isinstance(frame, str):
        if frame == 'equatorial':
            gc = SkyCoord(l=x * u.degree, b=y * u.degree, frame='galactic')
            x, y = (gc.gcrs.ra.deg, gc.gcrs.dec.deg)
        elif frame == 'galactic':
            pass
    elif isinstance(frame, SpacecraftFrame):
        flipped = False
        gc = SkyCoord(l=x * u.degree, b=y * u.degree, frame='galactic')
        gc = gc.transform_to(frame)
        frame = 'spacecraft'
        # counterintuitive, but it's based on how matplotlib specifies the 
        # plotting axes
        if zen:
            x, y = (gc.az.deg, gc.el.deg)
        else:
            x, y = (gc.az.deg, 90.0-gc.el.deg)
    else:
        raise TypeError('If frame is not a str, it must be a ' \
                        'SpacecraftFrame')
        
        
    line1 = sky_line(x, y, ax, color=outer_color, linewidth=3,
                     alpha=line_alpha, flipped=flipped, frame=frame)
    line2 = sky_line(x, y, ax, color=inner_color, linewidth=1,
                     alpha=line_alpha, flipped=flipped, frame=frame)

    # plot Galactic center
    pt1 = sky_point(x[0], y[0], ax, marker='o', c=outer_color, s=100,
                    alpha=center_alpha, edgecolor=None, flipped=flipped,
                    frame=frame)
    pt2 = sky_point(x[0], y[0], ax, marker='o', c=inner_color, s=20,
                    alpha=center_alpha, edgecolor=None, flipped=flipped,
                    frame=frame)

    return [line1, line2, pt1, pt2]


def sky_heatmap(x, y, heatmap, ax, cmap='RdPu', norm=None, flipped=True,
                frame='equatorial', **kwargs):
    """Plot a heatmap on the sky as a colormap gradient.
       
    Args:
        x (np.array): The azimuthal coordinate array of the heatmap grid
        y (np.array): The polar coordinate array of the heatmap grid
        heatmap (np.array): The heatmap array, of shape (x.size, y.size)
        ax (matplotlib.axes): Plot axes object
        cmap (str, optional): The colormap. Default is 'RdPu'
        norm (matplotlib.colors.Normalize or similar, optional):
            The normalization used to scale the colormapping to the heatmap 
            values. This can be initialized by Normalize, LogNorm, SymLogNorm, 
            PowerNorm, or some custom normalization.  Default is 
            PowerNorm(gamma=0.3).
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following the equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', or 'spacecraft'.
                               Default is 'equatorial'
        **kwargs: Other plotting options
    
    Returns:
        (matplotlib.collections.QuadMesh)
    """
    theta = np.deg2rad(y)
    phi = np.deg2rad(x - 180.0)
    
    if frame == 'spacecraft':
        flipped = False
        phi1 = phi - np.pi
        image1 = ax.pcolormesh(phi1, theta, heatmap, rasterized=True, cmap=cmap,
                               norm=norm)
        
        phi2 = phi + np.pi
        image2 = ax.pcolormesh(phi2, theta, heatmap, rasterized=True, cmap=cmap,
                               norm=norm)
        return [image1, image2]
    elif frame == 'galactic':
        phi1 = phi - np.pi
        if flipped:
            phi1 = -phi1
        image1 = ax.pcolormesh(phi1, theta, heatmap, rasterized=True, cmap=cmap,
                               norm=norm)
        
        phi2 = phi + np.pi
        if flipped:
            phi2 = -phi2
        image2 = ax.pcolormesh(phi2, theta, heatmap, rasterized=True, cmap=cmap,
                               norm=norm)
        return [image1, image2]
    else:
        if flipped:
            phi = -phi
        image = ax.pcolormesh(phi, theta, heatmap, rasterized=True, cmap=cmap,
                              norm=norm)
        return image
    

    return image

