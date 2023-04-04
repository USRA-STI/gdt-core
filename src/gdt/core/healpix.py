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
from scipy.stats import norm
from astropy.coordinates import SkyCoord

import healpy as hp
from matplotlib.pyplot import contour as Contour
from matplotlib.patches import Polygon
from gdt.core.plot.lib import circle as sky_circle

__all__ = ['HealPix', 'HealPixEffectiveArea', 'HealPixLocalization']

class HealPix():
    """Base class for HEALPix files
    """
    def __init__(self):
        self._hpx = []
        self._trigtime = None
        self._nside = 0

    @property
    def npix(self):
        """(int): Number of pixels in the HEALPix map"""
        return len(self._hpx)

    @property
    def nside(self):
        """(int): The HEALPix resolution"""
        return hp.npix2nside(self.npix)

    @property
    def pixel_area(self):
        """(float): The area of each pixel in square degrees"""
        return 4.0 * 180.0 ** 2 / (np.pi * self.npix)

    @property
    def trigtime(self):
        """(float): The reference time"""
        return self._trigtime

    def convolve(self, model, *args, **kwargs):
        """Convolve the map with a model kernel.  The model can be a Gaussian
        kernel or any mixture of Gaussian kernels. Uses `healpy.smoothing 
        <https://healpy.readthedocs.io/en/latest/generated/healpy.sphtfunc.smoothing.html>`_.
        
        An example of a model kernel with a 50%/50% mixture of two Gaussians,
        one with a 1-deg width, and the other with a 3-deg width::
            
            def gauss_mix_example():
                sigma1 = np.deg2rad(1.0)
                sigma2 = np.deg2rad(3.0)
                frac1 = 0.50
                return ([sigma1, sigma2], [frac1])
        
        Args: 
            model (<function>): The function representing the model kernel
            *args: Arguments to be passed to the model kernel function
        
        Returns:
            (:class:`HealPix`)
        """
        # evaluate model
        sigmas, fracs = model(*args)

        # determine number of gaussians, and ensure that they match the 
        # number of fractional weights
        num_sigmas = len(sigmas)
        if len(fracs) != num_sigmas:
            if len(fracs) + 1 != num_sigmas:
                raise ValueError(
                    'Number of mixture fraction parameters is incorrect')
            fracs.append(1.0 - np.sum(fracs))

        # for each gaussian, apply the smoothing at the prescribed weight
        new_hpx = np.zeros(self._hpx.shape)
        for i in range(num_sigmas):
            new_hpx += fracs[i] * hp.smoothing(self._hpx, sigma=sigmas[i])
                
        return type(self).from_data(new_hpx, trigtime=self.trigtime, **kwargs)

    @classmethod
    def from_data(cls, hpx_arr, trigtime=None, filename=None, **kwargs):
        """Create a HealPix object from healpix arrays
        
        Args:
            hpx_arr (np.array): The HEALPix array
            trigtime (float, optional): The reference time for the map
            filename (str, optional): The filename
        
        Returns:        
            (:class:`HealPix`)
        """
        obj = cls()
        try:
            iter(hpx_arr)
            hpx_arr = np.asarray(hpx_arr)
        except:
            raise TypeError('hpx_arr must be an array')
        obj._hpx = hpx_arr
        
        if trigtime is not None:
            try:
                trigtime = float(trigtime)
            except:
                raise TypeError('trigtime must be a non-negative float')
            if trigtime < 0.0:
                raise ValueError('trigtime must be non-negative')
        obj._trigtime = trigtime
                
        obj._filename = filename
        
        return obj

    @classmethod
    def multiply(cls, healpix1, healpix2, primary=0, output_nside=128, 
                 **kwargs):
        """Multiply two HealPix maps and return a new HealPix object
        
        Args:
            healpix1 (:class:`HealPix`): One of the HEALPix maps to multiply
            healpix2 (:class:`HealPix`): The other HEALPix map to multiply
            primary (int, optional): If 0, use the first map metadata, 
                                     or if 1, use the second map metadata. 
                                     Default is 0.
            output_nside (int, optional): The nside of the multiplied map. 
                                          Default is 128.
        Returns
            (:class:`HealPix`)
        """
        if not isinstance(healpix1, HealPix) or \
           not isinstance(healpix2, HealPix):
            raise TypeError('healpix1 and healpix2 must be HealPix objects')
        
        if primary != 0 and primary != 1:
            raise ValueError('primary must be either 0 or 1')
        
        # if different resolutions, upgrade the lower res, then multiply
        if healpix1.nside > healpix2.nside:
            hpx = healpix1._hpx * hp.ud_grade(healpix2._hpx,
                                                nside_out=healpix1.nside)
        elif healpix1.nside < healpix2.nside:
            hpx = healpix2._hpx * hp.ud_grade(healpix1._hpx,
                                                nside_out=healpix2.nside)
        else:
            hpx = healpix1._hpx * healpix2._hpx

        # output resolution and normalize
        hpx = hp.ud_grade(hpx, output_nside)

        # copy trigtime
        if primary == 0:
            trigtime = healpix1.trigtime
        else:
            trigtime = healpix2.trigtime
        
        return cls.from_data(hpx, trigtime=trigtime, **kwargs)

    @staticmethod
    def _dec_to_theta(dec):
        return np.deg2rad(90.0 - dec)

    @staticmethod
    def _phi_to_ra(phi):
        return np.rad2deg(phi)

    @staticmethod
    def _ra_to_phi(ra):
        return np.deg2rad(ra % 360.0)

    @staticmethod
    def _theta_to_dec(theta):
        return np.rad2deg(np.pi / 2.0 - theta)

    def _ang_to_pix(self, ra, dec):
        # convert RA/Dec to healpixels
        theta = self._dec_to_theta(dec)
        phi = self._ra_to_phi(ra)
        pix = hp.ang2pix(self.nside, theta, phi)
        return pix

    def _mesh_grid(self, num_phi, num_theta):
        try:
            numpts_phi = int(num_phi)
            numpts_theta = int(num_theta)
        except:
            raise TypeError('num_phi and num_theta must be positive integers')
        if num_phi <= 0 or num_theta <= 0:
            raise ValueError('num_phi and num_theta must be positive')
        
        # create the mesh grid in phi and theta
        theta = np.linspace(np.pi, 0.0, num_theta)
        phi = np.linspace(0.0, 2 * np.pi, num_phi)
        phi_grid, theta_grid = np.meshgrid(phi, theta)
        grid_pix = hp.ang2pix(self.nside, theta_grid, phi_grid)
        return (grid_pix, phi, theta)

    def __repr__(self):
        s = '<{0}: NSIDE={1}>'.format(self.__class__.__name__, self.nside)
        return s


class HealPixEffectiveArea(HealPix):
    """Class for effective area HEALPix files 
    """
    pass

    @property
    def eff_area(self):
        """(np.array): The effective area HEALPix array"""
        return self._hpx
    
    def effective_area(self, az, zen):
        """Calculate the effective area at a given point.  This
        function interpolates the map at the requested point rather than
        providing the vale at the nearest pixel center.
        
        Args:
            az (float): The Azimuth
            zen (float): The Zenith
        
        Returns:
            (float)
        """
        phi = self._ra_to_phi(az)
        theta = np.deg2rad(zen)
        eff_area = hp.get_interp_val(self.eff_area, theta, phi)
        return eff_area  
    
    def make_grid(self, numpts_az=360, numpts_zen=180):
        """Return the effective area mapped to a grid.
        
        Args:
            numpts_az (int, optional): The number of grid points along the Az 
                                       axis. Default is 360.
            numpts_zen (int, optional): The number of grid points along the Zen 
                                        axis. Default is 180.

        Returns: 
            3-tuple containing:
            
            - *np.array*: The effective area array with shape
                      (``numpts_zen``, ``numpts_az``)
            - *np.array*: The Az grid points
            - *np.array*: The Zen grid points
        """        
        grid_pix, phi, theta = self._mesh_grid(numpts_az, numpts_zen)
        return (self.eff_area[grid_pix], self._phi_to_ra(phi), np.rad2deg(theta))
    
    @classmethod
    def from_cosine(cls, az, zen, eff_area_normal, coeff=1.0, nside=64, 
                    filename=None, **kwargs):
        """Create an object with effective area that has a cosine angular
        depndence.
        
        Args:
            az (float): The azimuth of the detector normal
            zen (float): The zenith of the detector normal
            eff_area_normal (float): The effective area at the detector normal
            coeff (float, optional): The cosine coefficient to use. Default is 1
            nside (int, optional): The nside of the HEALPix to make. Default is 64
            filename (str, optional): The filename
        
        Returns:
            (:class:`HealPixEffectiveArea`)
        """
        try:
            az = float(az)
            zen = float(zen)
            eff_area_normal = float(eff_area_normal)
            coeff = float(coeff)
        except:
            raise TypeError('az, zen, eff_area_normal, and coeff must be floats')
    
        az = az % 360.0
        if zen < 0.0 or zen > 180.0:
            raise ValueError('zen must be between 0 and 180')

        if eff_area_normal < 0:
            raise ValueError('eff_area_normal must be non-negative')
        if coeff <= 0.0:
            raise ValueError('coeff must be positive')
       
        coord = SkyCoord(az, 90-zen, unit='deg')
        
        npix = hp.nside2npix(nside)
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        azs = cls._phi_to_ra(phi)
        zens = np.rad2deg(theta)
        coords = SkyCoord(azs, 90.0-zens, unit='deg')
        
        angles = coord.separation(coords)
        eff_area = eff_area_normal * np.cos(coeff * angles.to('rad')).value
        mask = coeff * angles.value > 90.0
        eff_area[mask] = 0.0
        return cls.from_data(eff_area, filename=filename, **kwargs)

    @classmethod
    def from_uniform(cls, eff_area, nside=64, filename=None, **kwargs):
        """Create an object with uniform effective area.
        
        Args:
            eff_area (float): The effective area.
            nside (int, optional): The nside of the HEALPix to make.
            filename (str, optional): The filename
        
        Returns:
            (:class:`HealPixEffectiveArea`)
        """
        try:
            eff_area = float(eff_area)
        except:
            raise TypeError('eff_area must be a non-negative float')
        if eff_area < 0:
            raise ValueError('eff_area must be non-negative')
        
        npix = hp.nside2npix(nside)
        hpx_arr = np.full(npix, eff_area)
        return cls.from_data(hpx_arr, filename=filename, **kwargs)

    @classmethod
    def sum(cls, healpix1, healpix2, primary=0, output_nside=128, **kwargs):
        """Add two HealPixEffectiveArea maps and return a new object
        
        Args:
            healpix1 (:class:`HealPixEffectiveArea`): One of the HEALPix maps 
                                                      to sum
            healpix2 (:class:`HealPixEffectiveArea`): The other HEALPix map 
                                                      to sum
            primary (int, optional): If 0, use the first map metadata, 
                                     or if 1, use the second map metadata. 
                                     Default is 0.
            output_nside (int, optional): The nside of the multiplied map. 
                                          Default is 128.
        Returns
            (:class:`HealPixEffectiveArea`)
        """
        if not isinstance(healpix1, HealPixEffectiveArea) or \
           not isinstance(healpix2, HealPixEffectiveArea):
            raise TypeError('healpix1 and healpix2 must be ' \
                            'HealPixEffectiveArea objects')
        
        if primary != 0 and primary != 1:
            raise ValueError('primary must be either 0 or 1')
        
        # if different resolutions, upgrade the lower res, then multiply
        if healpix1.nside > healpix2.nside:
            hpx = healpix1._hpx + hp.ud_grade(healpix2._hpx,
                                                nside_out=healpix1.nside)
        elif healpix1.nside < healpix2.nside:
            hpx = healpix2._hpx + hp.ud_grade(healpix1._hpx,
                                                nside_out=healpix2.nside)
        else:
            hpx = healpix1._hpx + healpix2._hpx

        # output resolution and normalize
        hpx = hp.ud_grade(hpx, output_nside)

        # copy trigtime info
        if primary == 0:
            trigtime = healpix1.trigtime
        else:
            trigtime = healpix2.trigtime
        
        return cls.from_data(hpx, trigtime=trigtime, **kwargs)

        
class HealPixLocalization(HealPix):
    """Class for localization HEALPix files 
    """
    def __init__(self):
        super().__init__()
        self._sig = None
    
    @property
    def centroid(self):
        """(float, float): The RA, Dec of the centroid"""
        pix = np.argmax(self.prob)
        theta, phi = hp.pix2ang(self.nside, pix)
        return (self._phi_to_ra(phi), self._theta_to_dec(theta))

    @property
    def prob(self):
        """(np.array): The HEALPix array for the probability/pixel"""
        return self._hpx
    
    @property
    def sig(self):
        """(np.array): The HEALPix array for the significance of each pixel"""
        return self._sig
    
    def area(self, clevel):
        """Calculate the sky area contained within a given confidence region
        
        Args:
            clevel (float): The localization confidence level (valid range 0-1)
        
        Returns:
            (float)
        """
        if clevel < 0.0 or clevel > 1.0:
            raise ValueError('clevel must be between 0 and 1')
        numpix = np.sum((1.0 - self.sig) <= clevel)
        return numpix * self.pixel_area

    def confidence(self, ra, dec):
        """Calculate the localization confidence level for a given point. 
        This function interpolates the map at the requested point rather than
        providing the value at the nearest pixel center.
        
        Args:
            ra (float): The RA
            dec (float): The Dec
        
        Returns:
            (float)
        """
        phi = self._ra_to_phi(ra)
        theta = self._dec_to_theta(dec)
        return 1.0 - hp.get_interp_val(self.sig, theta, phi)

    def confidence_region_path(self, clevel, numpts_ra=360, numpts_dec=180):
        """Return the bounding path for a given confidence region.
        
        Args:
            clevel (float): The localization confidence level (valid range 0-1)
            numpts_ra (int, optional): The number of grid points along the RA 
                                       axis. Default is 360.
            numpts_dec (int, optional): The number of grid points along the Dec 
                                        axis. Default is 180.
        
        Returns:
            ([(np.array, np.array), ...]): A list of RA, Dec points, where each
                                           item in the list is a continuous 
                                           closed path.
        """
        if clevel < 0.0 or clevel > 1.0:
            raise ValueError('clevel must be between 0 and 1')

        # create the grid and integrated probability array
        grid_pix, phi, theta = self._mesh_grid(numpts_ra, numpts_dec)
        sig_arr = 1.0 - self.sig[grid_pix]
        ra = self._phi_to_ra(phi)
        dec = self._theta_to_dec(theta)

        # use matplotlib contour to produce a path object
        contour = Contour(ra, dec, sig_arr, [clevel])

        # get the contour path, which is made up of segments
        paths = contour.collections[0].get_paths()

        # extract all the vertices
        pts = [path.vertices for path in paths]

        # unfortunately matplotlib will plot this, so we need to remove
        for c in contour.collections:
            c.remove()

        return pts

    def probability(self, ra, dec, per_pixel=False):
        """Calculate the localization probability at a given point.  This
        function interpolates the map at the requested point rather than
        providing the vale at the nearest pixel center.
        
        Args:
            ra (float): The RA
            dec (float): The Dec
            per_pixel (bool, optional): 
                If True, return probability per pixel, otherwise return 
                probability per square degree. Default is False.
        
        Returns:
            (float)
        """
        phi = self._ra_to_phi(ra)
        theta = self._dec_to_theta(dec)
        prob = hp.get_interp_val(self.prob, theta, phi)
        if not per_pixel:
            prob /= self.pixel_area
        return prob

    def prob_array(self, numpts_ra=360, numpts_dec=180, sqdegrees=True,
                   sig=False):
        """Return the localization probability mapped to a grid on the sky
        
        Args:
            numpts_ra (int, optional): The number of grid points along the RA 
                                       axis. Default is 360.
            numpts_dec (int, optional): The number of grid points along the Dec 
                                        axis. Default is 180.
            sqdegrees (bool, optional): 
                If True, the prob_array is in units of probability per square 
                degrees, otherwise in units of probability per pixel. 
                Default is True
            sig (bool, optional): Set True to retun the significance map on a 
                                  grid instead of the probability. 
                                  Default is False.

        Returns: 
            3-tuple containing:
            
            - *np.array*: The probability (or significance) array with shape \
                      (``numpts_dec``, ``numpts_ra``)
            - *np.array*: The RA grid points
            - *np.array*: The Dec grid points
        """        
        grid_pix, phi, theta = self._mesh_grid(numpts_ra, numpts_dec)

        if sig:
            sqdegrees = False
            prob_arr = self.sig[grid_pix]
        else:
            prob_arr = self.prob[grid_pix]
        if sqdegrees:
            prob_arr /= self.pixel_area
        return (prob_arr, self._phi_to_ra(phi), self._theta_to_dec(theta))

    def region_probability(self, healpix, prior=0.5):
        r"""The probability that the HealPix localization is associated with
        another HealPixLocalization map.  This is calculated against the null 
        hypothesis that the two HealPix maps are unassociated:
        
        :math:`P(A | \mathcal{I}) = 
        \frac{P(\mathcal{I} | A) \ P(A)}
        {P(\mathcal{I} | A) \ P(A) + P(\mathcal{I} | \neg A) \ P(\neg A)}`
        
        where
        
        * :math:`P(\mathcal{I} | A)` is the integral over the overlap of the two 
          maps once the Earth occultation has been removed for *this* map.
        * :math:`P(\mathcal{I} | \neg A)` is the integral over the overlap of
          *this* map with a uniform distribution on the sky (i.e. the probability 
          the localization is associated with a random point on the sky)
        * :math:`P(A)` is the prior probability that *this* localization is 
          associated with the *other* HEALPix map.
        
        Args:
            healpix (:class:`HealPixLocalization`): The healpix map for which to 
                                                    calculate the spatial 
                                                    association
            prior (float, optional): The prior probability that the localization
                                     is associated with the source. 
                                     Default is 0.5
        
        Returns:  
            (float)
        """
        if (prior < 0.0) or (prior > 1.0):
            raise ValueError('Prior probability must be within 0-1, inclusive')
        # convert uniform prob/sr to prob/pixel
        u = 1.0 / (4.0 * np.pi)

        # ensure maps are the same resolution
        probmap1 = self.prob
        probmap2 = healpix.prob
        if self.nside > healpix.nside:
            probmap2 = hp.ud_grade(probmap2, nside_out=self.nside)
            probmap2 = self._assert_prob(probmap2)
            u *= hp.nside2resol(self.nside) ** 2
        elif self.nside < healpix.nside:
            probmap1 = hp.ud_grade(probmap1, nside_out=healpix.nside)
            probmap1 = self._assert_prob(probmap1)
            u *= hp.nside2resol(healpix.nside) ** 2
        else:
            u *= hp.nside2resol(self.nside) ** 2

        # alternative hypothesis: they are related
        alt_hyp = (probmap1 * probmap2).sum()
        # null hypothesis: one of the maps is from an unassociated source
        # (uniform spatial probability)
        null_hyp = (probmap1 * u).sum()

        # since we have an exhaustive and complete list of possibilities, we can
        # easily calculate the probability
        prob = (alt_hyp*prior) / ((alt_hyp*prior) + (null_hyp*(1.0-prior)))
        return prob

    def source_probability(self, ra, dec, prior=0.5):
        r"""The probability that the HealPix localization is associated with
        a known point location.  This is calculated against the null hypothesis
        that the HealPix localization originates from an unassociated random
        source that has equal probability of origination anywhere in the sky:
        
        :math:`P(A | \mathcal{I}) = 
        \frac{P(\mathcal{I} | A) \ P(A)}
        {P(\mathcal{I} | A) \ P(A) + P(\mathcal{I} | \neg A) \ P(\neg A)}`
        
        where
        
        * :math:`P(\mathcal{I} | A)` is the probability of the localization at
          the point source once
        * :math:`P(\mathcal{I} | \neg A)` is the probability per pixel assuming 
          a uniform distribution on the sky (i.e. the probability the 
          localization is associated with a random point on the sky)
        * :math:`P(A)` is the prior probability that the localization is 
          associated with the point source
        
        Args:
            ra (float): The RA of the known source location
            dec (float): The Dec of the known source location
            prior (float, optional): The prior probability that the localization
                                     is associated with the source. 
                                     Default is 0.5
        
        Returns:        
            (float)
        """
        if (prior < 0.0) or (prior > 1.0):
            raise ValueError('Prior probability must be within 0-1, inclusive')
        # convert uniform prob/sr to prob/pixel
        u = 1.0 / (4.0 * np.pi)
        u *= hp.nside2resol(self.nside) ** 2

        # the pixel probability of the skymap at the location of the point source
        p = self.probability(ra, dec, per_pixel=True)

        # null hypothesis is that they are not associated, therefore the sky map
        # is result of some source that has uniform probability on the sky
        prob = (p*prior) / ((p*prior) + (u*(1.0-prior)))
        return prob

    @classmethod
    def from_annulus(cls, center_ra, center_dec, radius, sigma, nside=None,
                     trigtime=None, filename=None, **kwargs):
        """Create a HealPixLocalization object of a Gaussian-width annulus.
        
        Args:
            center_ra (float): The RA of the center of the annulus
            center_dec (float): The Dec of the center of the annulus
            radius (float): The radius of the annulus, in degrees, measured to 
                            the center of the of the annulus
            sigma (float): The Gaussian standard deviation width of the annulus, 
                           in degrees
            nside (int, optional): The nside of the HEALPix to make. By default,
                                   the nside is automatically determined by the 
                                   ``sigma`` width.  Set this argument to 
                                   override the default. 
            trigtime (float, optional): The reference time for the map
            filename (str, optional): The filename
                
        Return:
            (:class:`HealPixLocalization`)
        """
        try:
            center_ra = float(center_ra)
            center_dec = float(center_dec)
            radius = float(radius)
            sigma = float(sigma)
        except:
            raise TypeError('center_ra, center_dec, radius, and sigma must be' \
                            ' floats')
        center_ra = center_ra % 360.0
        if center_dec < -90.0 or center_dec > 90.0:
            raise ValueError('center_dec must be between -90 and 90')
        if radius < 0:
            raise ValueError('radius must be positive')
        if sigma < 0:
            raise ValueError('sigma must be positive')
        
        # Automatically calculate appropriate nside by taking the closest nside
        # with an average resolution that matches 0.2*sigma
        if nside is None:
            nsides = 2**np.arange(15)
            pix_res = hp.nside2resol(nsides, True)/60.0
            idx = np.abs(pix_res-sigma/5.0).argmin()
            nside = nsides[idx]
        
        # get everything in the right units
        center_phi = cls._ra_to_phi(center_ra)
        center_theta = cls._dec_to_theta(center_dec)
        radius_rad = np.deg2rad(radius)
        sigma_rad = np.deg2rad(sigma)

        # number of points in the circle based on the approximate arclength 
        # and resolution
        res = hp.nside2resol(nside)
        
        # calculate normal distribution about annulus radius with sigma width
        x = np.linspace(0.0, np.pi, int(10.0*np.pi/res))
        pdf = norm.pdf(x, loc=radius_rad, scale=sigma_rad)

        # cycle through annuli of radii from 0 to 180 degree with the 
        # appropriate amplitude and fill the probability map
        probmap = np.zeros(hp.nside2npix(nside))
        for i in range(x.size):
            # no need to waste time on pixels that will have ~0 probability...
            if pdf[i]/pdf.max() < 1e-10:
                continue
            
            # approximate arclength determines number of points in each annulus
            arclength = 2.0*np.pi*x[i]
            numpts = int(np.ceil(arclength/res))*10
            circ = sky_circle(x[i], center_phi, center_theta, num_points=numpts)
            theta = np.pi / 2.0 - circ[1]
            phi = circ[0]
            
            # convert to pixel indixes and fill the map
            idx = hp.ang2pix(nside, theta, phi)
            probmap[idx] = pdf[i]
            mask = (probmap[idx] > 0.0)
            probmap[idx[~mask]] = pdf[i]
            probmap[idx[mask]] = (probmap[idx[mask]] + pdf[i])/2.0

        obj = cls.from_data(probmap, trigtime=trigtime, filename=filename, 
                            **kwargs)
        return obj

    @classmethod
    def from_data(cls, prob_arr, trigtime=None, filename=None, **kwargs):
        """Create a HealPixLocalization object from a HEALPix probability array.
        
        Args:
            prob_arr (np.array): The HEALPix array
            trigtime (float, optional): The reference time for the map
            filename (str, optional): The filename
        
        Returns:        
            (:class:`HealPixLocalization`)
        """
        obj = super().from_data(prob_arr, trigtime=trigtime, filename=filename,
                                **kwargs)
        obj._hpx = obj._assert_prob(obj._hpx)
        obj._sig = obj._assert_sig(1.0 - cls._credible_levels(obj.prob))
        return obj

    @classmethod
    def from_gaussian(cls, center_ra, center_dec, sigma, nside=None, 
                      trigtime=None, filename=None, **kwargs):
        """Create a HealPixLocalization object of a Gaussian
        
        Args:
            center_ra (float): The RA of the center of the Gaussian
            center_dec (float): The Dec of the center of the Gaussian
            sigma (float): The Gaussian standard deviation, in degrees
            nside (int, optional): The nside of the HEALPix to make. By default,
                                   the nside is automatically determined by the 
                                   `sigma` of the Gaussian.  Set this argument 
                                   to override the default. 
            trigtime (float, optional): The reference time for the map
            filename (str, optional): The filename
        
        Returns:
            (:class:`HealPixLocalization`)
        """
        try:
            center_ra = float(center_ra)
            center_dec = float(center_dec)
            sigma = float(sigma)
        except:
            raise TypeError('center_ra, center_dec, and sigma must be floats')
        center_ra = center_ra % 360.0
        if center_dec < -90.0 or center_dec > 90.0:
            raise ValueError('center_dec must be between -90 and 90')
        if sigma < 0:
            raise ValueError('sigma must be positive')

        # Automatically calculate appropriate nside by taking the closest nside
        # with an average resolution that matches 0.2*sigma
        if nside is None:
            nsides = 2**np.arange(15)
            pix_res = hp.nside2resol(nsides, True)/60.0
            idx = np.abs(pix_res-sigma/10.0).argmin()
            nside = nsides[idx]
        
        # get everything in the right units
        center_phi = cls._ra_to_phi(center_ra)
        center_theta = cls._dec_to_theta(center_dec)
        sigma_rad = np.deg2rad(sigma)

        # point probability
        npix = hp.nside2npix(nside)
        probmap = np.zeros(npix)
        probmap[hp.ang2pix(nside, center_theta, center_phi)] = 1.0

        # then smooth out using appropriate gaussian kernel
        probmap = hp.smoothing(probmap, sigma=sigma_rad)

        obj = cls.from_data(probmap, trigtime=trigtime, filename=filename, 
                            **kwargs)
        return obj

    @classmethod
    def from_vertices(cls, ra_pts, dec_pts, nside=64, trigtime=None, 
                      filename=None, **kwargs):
        """Create a HealPixLocalization object from a list of RA, Dec vertices.
        The probability within the vertices will be distributed uniformly and
        zero probability outside the vertices.
        
        Args:
            ra_pts (np.array): The array of RA coordinates
            dec_pts (np.array): The array of Dec coordinates
            nside (int, optional): The nside of the HEALPix to make. Default is 64.
            trigtime (float, optional): The reference time for the map
            filename (str, optional): The filename
        
        Returns:
            (:class:`HealPixLocalization`)
        """
        ra_pts = np.asarray(ra_pts) % 360.0
        dec_pts = np.asarray(dec_pts)
        if np.any((dec_pts < -90.0) | (dec_pts > 90.0)):
            raise ValueError('all dec_pts must be between -90 and 90')
                
        poly = Polygon(np.vstack((ra_pts, dec_pts)).T, closed=True)

        npix = hp.nside2npix(nside)
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        ra = cls._phi_to_ra(phi)
        dec = cls._theta_to_dec(theta)
        mask = poly.contains_points(np.vstack((ra, dec)).T)

        probmap = np.zeros(npix)
        probmap[mask] = 1.0
        obj = cls.from_data(probmap, trigtime=trigtime, filename=filename, 
                            **kwargs)
        return obj

    @staticmethod
    def _credible_levels(p):
        """Calculate the credible levels of a probability array using a greedy
        algorithm.
    
        Args:
            p (np.array): The probability array
    
        Returns:    
             (np.array)
        """
        p = np.asarray(p)
        p_flat = p.flatten()
        idx = np.argsort(p_flat)[::-1]
        clevels = np.empty_like(p_flat)
        clevels[idx] = np.cumsum(p_flat[idx])
        return clevels.reshape(p.shape)

    def _assert_prob(self, prob):
        # ensure that the pixels have valid probability:
        # each pixel must be > 0 and sum == 1.
        prob[prob < 0.0] = 0.0
        prob /= prob.sum()
        return prob

    def _assert_sig(self, sig):
        # ensure that the pixels have valid significance:
        # each pixel must have significance [0, 1]
        if sig is not None:
            sig[sig < 0.0] = 0.0
            sig[sig > 1.0] = 1.0
        return sig

    def __repr__(self):
        s = '<{0}: \n'.format(self.__class__.__name__)
        s += ' NSIDE={0}; trigtime={1};\n'.format(self.nside, self.trigtime)
        s += ' centroid={}>'.format(self.centroid)
        return s
