# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Developed by: Joshua Wood
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Developed by: Jacob Smith
#               University of Alabama in Huntsville
#               Center for Space Plasma and Aeronomic Research
#
# This software is not subject to EAR.
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
import matplotlib.lines
import matplotlib.pyplot as plt
import matplotlib.patches

from astropy.utils import isiterable
from rich.progress import track
from astropy import units as u

from .plot.lightcurve import Lightcurve
from .binning.unbinned import bin_by_time
from .binning.binned import rebin_by_edge_index, combine_by_factor

__all__ = ['TriggerAlgorithm', 'Trigger', 'SlidingWindowMethod']

class TriggerAlgorithm:
    """Class for holding information related to a trigger algorithm

    Parameters:
        timescale (int): The duration in integer milliseconds of the window
                         over which we count events contributing to a trigger
        offset (int):    The time offset in integer milliseconds between the
                         start time of the trigger window and the beginning of
                         instrument data
        channels (list|tuple): The energy channel range given as
                         (channel_start, channel_stop) for counts contributing
                         to the trigger window
        threshold (float): The signficance above background for including
                           a detector in the trigger
    """

    def __init__(self, timescale: int, offset: int, channels: list | tuple, threshold: float):

        # apply sanity checks
        if not isinstance(timescale, int):
            raise ValueError("Timescale must be an integer number of milliseconds")
        if not isinstance(offset, int):
            raise ValueError("Offset must be an integer number of milliseconds")
        if channels[1] < channels[0]:
            raise ValueError("channels[1] must be >= channels[0]")

        self._timescale = timescale
        self._offset = offset
        self._channels = list(channels)
        self._threshold = threshold

    @property
    def timescale(self):
        """(int): The trigger window duration in integer milliseconds"""
        return self._timescale

    @property
    def offset(self):
        """(int): The offset of the first trigger window relative to data start time in integer milliseconds"""
        return self._offset

    @property
    def channels(self):
        """(list): The energy channel range contributing to the trigger"""
        return self._channels

    @property
    def threshold(self):
        """(float): The significance level above background for including a detector in the trigger"""
        return self._threshold

    def timescale_in_unit(self, unit):
        """Method for retrieving the trigger timescale with an astropy unit
 
        Args:
            unit (astropy.units): The desired unit for returning the trigger timescale
                                  (e.g. astropy.units.second)

        Returns:
            (astropy.units)
        """
        return (self._timescale * u.ms).to(unit)

    def offset_in_unit(self, unit):
        """Method for retrieving the trigger offset with an astropy unit
 
        Args:
            unit (astropy.units): The desired unit for returning the trigger offset
                                  (e.g. astropy.units.second)

        Returns:
            (astropy.units)
        """
        return (self._offset * u.ms).to(unit)

    def __repr__(self):
        return "TriggerAlgorithm(timescale {:4d} ms, offset {:4d} ms, channels {}, threshold {:.2f} sigma)".format(
            self._timescale, self._offset, self._channels, self._threshold)


class Trigger:
    """Class for storing trigger information

    Parameters:
        alg_num (int):            Trigger algorithm number
        alg (TriggerAlgorithm):   Class defining the algorithm properties
        time (float):             Trigger time in seconds
        sig (np.array):           Significance over background in each detector
        triggered_det (np.array): Array with values of True for detectors
                                  contributing to the trigger, False otherwise
        detector_names (list):    List of detector names associated with
                                  entries in triggered_det
    """

    def __init__(self, alg_num: int, alg: TriggerAlgorithm, time: float,
                 sig: np.array, triggered_det: np.array, detector_names: list):
        self.alg_num = alg_num
        self.alg = alg
        self.time = time
        self.sig = sig
        self.triggered_det = triggered_det
        self.detector_names = detector_names

    @property
    def triggered_det_num(self):
        """(list): List of detector numbers contributing to this trigger"""
        return [i for i in range(self.triggered_det.size) if self.triggered_det[i]]

    @property
    def triggered_det_names(self):
        """(list): List of detector names contributing to this trigger"""
        return [self.detector_names[i] for i in self.triggered_det_num]

    def __repr__(self):
        """(str): String representation of trigger information"""
        return 'Trigger <algorithm num {}, timescale {} ms, time {} s, triggered det {}>'.format(
            self.alg_num, self.alg.timescale, self.time, self.triggered_det_names)


class SlidingWindowMethod:
    """Class for applying a sliding window trigger which identifies
    transients by looking for windows where n detectors to exceed
    a given significance threshold above background.

    Parameters:
        algorithms (dict):       Dictionary defining the set of trigger
                                 algorithms applied in the trigger method
        background_window (int): Length of time in integer milliseconds
                                 used to compute background counts
        background_offset (int): Time offset in integer milliseconds between
                                 the end of the trigger and background windows
        resolution (int):        The time resolution in integer milliseconds
                                 of binned PHAII data used for the trigger
                                 method. All other timescales (windows,
                                 offsets) should be multiples of this value
        channel_edges (list):    List of channel edges used to define energy
                                 bins for PHAII data
        det_names (list):        Names of required detector numbers for
                                 the trigger method
        n (int):                 Minimum number of detectors needed for trigger
        verbose (float):         Show progress bars when True
    """

    def __init__(self, algorithms: dict, background_window: int,
                 background_offset: int, resolution: int,
                 channel_edges: list, det_names: list,
                 n: int = 1, verbose: bool = True):

        self.n = n
        self.verbose = verbose
        self.resolution = resolution
        self.algorithms = algorithms
        self.channel_edges = channel_edges
        self.det_names = det_names

        if background_window % self.resolution != 0 or background_window <= 0:
            raise ValueError("Background window must be an integer multiple of resolution in ms")
        if background_offset % self.resolution != 0:
            raise ValueError("Background offset must be an integer multiple of resolution in ms")

        self.background_window = background_window
        self.background_offset = background_offset

        self.background_window_bins = self.background_window // self.resolution
        self.background_offset_bins = self.background_offset // self.resolution

        self.phaiis = None

    @property
    def algorithms(self):
        """(dict): Dictionary defining the set of trigger algorithms"""
        return self._algorithms

    @algorithms.setter
    def algorithms(self, algorithms: dict):
        """Method to set the algorithms used when applying the trigger method.

        Note: This performs sanity checks to ensure the algorithms are passed
        as a dictionary and are matched to the minimum resolution used to bin
        data files.

        Args:
            algorithms (dict): Dictionary defining the set of trigger
                               algorithms applied in the trigger method.
        """
        if not isinstance(algorithms, dict):
            raise ValueError("algorithms should be a dictionary object with the format int:TriggerAlgorithm()")
        for alg_num, alg in algorithms.items():
            if alg.timescale % self.resolution != 0:
                raise ValueError("Trigger algorithm timescale must be an integer multiple of resolution in ms")
            if alg.offset % self.resolution != 0:
                raise ValueError("Trigger algorithm offset must be an integer multiple of resolution in ms")

        self._algorithms = algorithms

    def apply_trigger(self, holdoff: float = None, debug: bool = False):
        """Applies the set of trigger algorithms to prepared data

        Args:
            holdoff (float, optional): Trigger holdoff time in seconds.
                                       Set to 300 sec to replicate GCN
                                       trigger notices.
            debug (bool, optional):    Show debugging information when True

        Returns:
            (list of Trigger)
        """
        if self.phaiis is None:
            raise ValueError("Need to run prepare_data() method on TTE data before applying the trigger")

        # pre-compute algorithm values in terms of the phaii bin resolution
        bins = {}
        for alg_num, alg in self.algorithms.items():
            bins[alg_num] = {"window": alg.timescale // self.resolution,
                             "offset": alg.offset // self.resolution}

        # loop over time bins and apply all trigger algorithms whose
        # time windows end at the current time bin. This allows us to
        # to compute the background counts once for all algorithms
        # that need to be evaluated.
        triggers = []
        description = "Applying algorithms"
        nbins = self.phaiis[0].data.size[0] # total number of time bins
        for i in track(range(nbins), description=description, disable=not self.verbose):

            # i represents the index of the last bin in the trigger window.
            # the stop index when summing window bins should be this +1.
            window_stop_index = i + 1

            # calculate background window indices for this time, which starts at
            # self.background_offset before the end of the trigger window
            background_stop_index = window_stop_index - self.background_offset_bins
            background_start_index = background_stop_index - self.background_window_bins

            # skip cases where we can't fully form the background window
            if background_start_index < 0:
                continue

            # compute background counts and exposure
            # note: we'll update background counts during the algorithm loop if the energy channels change
            background_channels = list(self.algorithms.values())[0].channels
            background_counts, background_exposure = [], []
            for phaii in self.phaiis:
                background_counts.append(phaii.data.counts[background_start_index:background_stop_index,
                                                           background_channels[0]:background_channels[1]+1].sum())
                background_exposure.append(phaii.data.exposure[background_start_index:background_stop_index].sum())

            # loop to check trigger window counts against background for each algorithm
            for alg_num, alg in self.algorithms.items():
 
                window_bins = bins[alg_num]['window']
                offset_bins = bins[alg_num]['offset']
                window_start_index = window_stop_index - window_bins

                # skip times with out-of-phase or incomplete trigger window
                if ((window_start_index - offset_bins) % window_bins) != 0 or i < window_bins:
                    continue

                # update background channels if they differ from current alg
                if alg.channels != background_channels:
                    background_channels = alg.channels
                    background_counts = [phaii.data.counts[background_start_index:background_stop_index,
                                                           background_channels[0]:background_channels[1]+1].sum()
                                         for phaii in self.phaiis]

                # optional debugging info
                if debug:
                    t0 = self.phaiis[0].data.tstart[0]
                    tstart = self.phaiis[0].data.tstart[window_start_index]
                    tstop = self.phaiis[0].data.tstop[window_stop_index - 1]
                    bkg_tstart = self.phaiis[0].data.tstart[background_start_index]
                    bkg_tstop = self.phaiis[0].data.tstop[background_stop_index - 1]

                    print(f"Evaluating Alg #{alg_num}: Timescale {alg.timescale*0.001:.3f} s, Offset {alg.offset*0.001:.3f} s")
                    print(f"    Window ({tstart-t0:.3f} s, {tstop-t0:.3f}) s, Delta {tstop-tstart:.3f} s")
                    print(f"    Background Window ({bkg_tstart-t0:.3f}, {bkg_tstop-t0:.3f}) s, Background Delta {bkg_tstop-bkg_tstart:.3f} s, Offset {tstop - bkg_tstop:.3f} s")

                # sum phaii bins to determine counts and exposure within trigger window
                counts, exposure = [], []
                for phaii in self.phaiis:
                    counts.append(phaii.data.counts[window_start_index:window_stop_index,
                                                    alg.channels[0]:alg.channels[1]+1].sum())
                    exposure.append(phaii.data.exposure[window_start_index:window_stop_index].sum())

                # compute significance under a simple Gaussian approximation.
                # this ignores more rigorous statistical formulations like
                # Li & Ma ApJ 272 (1983). This is Ok since the thresholds
                # are empirically tuned to limit the false positive rate.
                scaled_background = np.array(exposure) / np.array(background_exposure) * np.array(background_counts)
                excess = np.array(counts) - scaled_background
                sig = excess / np.sqrt(scaled_background)

                triggered_det = sig > alg.threshold
                if triggered_det.sum() >= self.n:
                    trigtime = self.phaiis[0].data.tstop[window_stop_index - 1]
                    triggers.append(Trigger(alg_num, alg, trigtime, sig, triggered_det, self.det_names))
                    if debug:
                       print("Found ", str(triggers[-1]))

        if holdoff:
            triggers = self.apply_holdoff(triggers, holdoff)

        return triggers

    def apply_holdoff(self, triggers: list, holdoff: float):
        """Apply a holdoff which ignores triggers occurring
        within the holdoff timescale of a prior trigger.

        Args:
            triggers (list): List of Trigger objects
            holdoff (float): Trigger holdoff time in seconds

        Returns:
            (list)
        """
        # ensure we're sorted by trigger time
        triggers.sort(key=lambda x: x.time)

        # loop to select triggers outside the holdoff timescale of prior triggers
        selected_triggers = []
        for trigger in triggers:
            if len(selected_triggers):
                # replace existing trigger with the more significant one
                if trigger.time == selected_triggers[-1].time and \
                    selected_triggers[-1].sig.max() < trigger.sig.max():
                    selected_triggers[-1] = trigger
                # append new trigger
                if trigger.time > selected_triggers[-1].time + holdoff:
                    selected_triggers.append(trigger)
            else:
                # always append the first trigger
                selected_triggers.append(trigger)

        return selected_triggers

    def prepare_data(self, ttes: list, time_range: list = None):
        """Method to prepare TTE data into a PHAII format for use with
        trigger algorithms.

        Args:
            ttes (list): List of Tte objects
            time_range (list): List with [tstart, tstop]
        """
        if len(ttes) != len(self.det_names):
            raise ValueError("Mismatch between detector names and number of TTE files")

        # ensure time range is within the range of provided data
        if time_range is None:
            time_range = list(ttes[0].time_range)
        for tte in ttes:
            if tte.time_range[0] > time_range[0]:
                time_range[0] = tte.time_range[0]
            if tte.time_range[1] < time_range[1]:
                time_range[1] = tte.time_range[1]

        self.phaiis = []
        description = "Binning TTE into {} energy channels".format(len(self.channel_edges)-1)
        for tte in track(ttes, description=description, disable=not self.verbose):
            phaii = tte.to_phaii(bin_by_time, self.resolution * 0.001, time_ref=0, time_range=time_range)
            phaii = phaii.rebin_energy(rebin_by_edge_index, np.array(self.channel_edges))
            self.phaiis.append(phaii)

    def waterfall_plot(self, triggers: list, **kwargs):
        """Create a waterfall plot showing all trigger windows as a function
        of algorithm number versus time.

        Args:
            triggers (list):   List of triggers
            kwargs (optional): Optional arguments passed to rectangular
                               patches drawn for each trigger

        Returns:
            (list)
        """
        # default plotting options
        if not len(kwargs):
            kwargs = {'facecolor':'C0', 'edgecolor':'b', 'linewidth':3}

        rects = []
        for trigger in triggers:
            width = trigger.alg.timescale_in_unit(u.second).value
            height = 1
            xy = (trigger.time - width, trigger.alg_num - 0.5)
            rects.append(plt.Rectangle(xy, width, height, **kwargs))
            plt.gca().add_patch(rects[-1])

        plt.xlabel("Time [s]")
        plt.ylabel("Algorithm Number")
        plt.gca().autoscale_view()

        return rects

    def lightcurve_plot(self, trigger: Trigger, time_range: tuple = None, detectors: list = None,
                        figsize: tuple = None, resolution: int = None, legend: bool = True):
        """Plot detector light curves with trigger overlay.

        Args:
            trigger (Trigger):  The trigger to overlay
            time_range (tuple): Start and end of time range to plot in seconds
            detectors (list):   List of detector numbers to plot
            figsize (tuple):    Dimensions of the figure in inches
            resolution (int):   Lightcurve resolution in integer milliseconds
            legend (bool):      Show legend when True

        Returns:
            (tuple)
        """
        if self.phaiis is None:
            raise ValueError("Need to run prepare_data() method on TTE data before applying the trigger")

        # default margins
        top = 0.9
        bottom = 0.1

        if detectors is None:
            detectors = np.arange(trigger.triggered_det.size)
        elif len(detectors) == 0:
            raise ValueError("Detector size must be greater than zero")
        if figsize is None:
            topsize = 0.5
            bottomsize = 1.0 if legend else 0.5
            figsize = [12, 2 * len(detectors) + topsize + bottomsize]
            # fix margins to prevent excessive white space
            top = 1 - topsize / figsize[1]
            bottom = bottomsize / figsize[1]
        if resolution is None:
            resolution = trigger.alg.timescale
        if resolution < self.resolution:
            raise ValueError(f"Lightcurve resolution cannot be less than PHAII resolution ({self.resolution} ms)")
        if time_range is None:
            time_range = (self.phaiis[0].data.tstart[0], self.phaiis[0].data.tstop[-1])

        fig, axes = plt.subplots(len(detectors), 1, sharex=True, sharey=False,
                                 figsize=figsize)
        fig.subplots_adjust(hspace=0.06, wspace=0.15, top=top, bottom=bottom)

        try:
            axes[0]
        except TypeError:
            axes = [axes]

        [ax.tick_params(axis='both', which='both', direction='in') for ax in axes]
        axes[0].set_title(trigger)

        bin_factor = resolution // self.resolution
        channel_range = trigger.alg.channels

        # convert resolution to seconds for future calculations
        resolution *= 0.001

        lcplots = []
        for i, det_num in enumerate(detectors):

            phaii = self.phaiis[det_num]

            # match lightcurve start/stop to data range
            tstart = max(time_range[0], phaii.data.tstart[0])
            tstop = min(time_range[1], phaii.data.tstop[-1])

            # ensure start/stop is a multiple of resolution from the trigger time
            tstart = trigger.time + int((tstart - trigger.time)/resolution) * resolution
            tstop = trigger.time + int((tstop - trigger.time)/resolution) * resolution

            # match range to nearest bins to prevent float rounding
            # issues when rebinning lightcurves by factor
            tstart = phaii.data.tstart[np.fabs(tstart - phaii.data.tstart).argmin()]
            tstop = phaii.data.tstop[np.fabs(tstop - phaii.data.tstop).argmin()]

            # create lightcurve with the appropriate resolution and range
            lc = phaii.to_lightcurve(channel_range=channel_range, time_range=(tstart, tstop))
            lc = lc.rebin(combine_by_factor, bin_factor, tstart=tstart, tstop=tstop)
            lcplots.append(Lightcurve(data=lc, ax=axes[i]))

            # trigger window bounds and exposure
            trigger_stop = trigger.time
            trigger_start = trigger_stop - trigger.alg.timescale * 0.001

            start_bin = np.fabs(trigger_start - phaii.data.tstart).argmin()
            stop_bin = np.fabs(trigger_stop - phaii.data.tstop).argmin()

            exposure = phaii.data.exposure[start_bin:stop_bin+1].sum()

            # background window bounds and exposure / counts
            background_stop = trigger.time - self.background_offset * 0.001
            background_start = background_stop - self.background_window * 0.001

            start_bin = np.fabs(background_start - phaii.data.tstart).argmin()
            stop_bin = np.fabs(background_stop - phaii.data.tstop).argmin()

            background_counts = phaii.data.counts[start_bin:stop_bin+1, channel_range[0]:channel_range[1]+1].sum()
            background_exposure = phaii.data.exposure[start_bin:stop_bin+1].sum()

            # trigger threshold
            scaled_background = background_counts * exposure / background_exposure
            threshold = scaled_background + trigger.alg.threshold * np.sqrt(scaled_background)

            # plot the background fit and window
            bg, = axes[i].plot([background_start, background_stop], 2 * [background_counts / background_exposure], color='r')
            bg_span = axes[i].axvspan(background_start, background_stop, color='r', alpha=0.1)

            # plot the trigger threshold and window
            trig_thresh, = axes[i].plot([trigger_start - 0.5 * trigger.alg.timescale * 0.001,
                                         trigger_stop + 0.5 * trigger.alg.timescale * 0.001],
                                        2 * [threshold / exposure], color='b')
            trig_win = axes[i].axvspan(trigger_start, trigger_stop, color='b', alpha=0.2)

            if not trigger.triggered_det[det_num]:
                axes[i].set_facecolor('#BEC0BF')
            if i != len(detectors) - 1:
                axes[i].set_xlabel("")
            axes[i].minorticks_on()

            # detector label
            x = 0.01
            y = 0.87
            label = self.det_names[i] + int(trigger.triggered_det[det_num]) * " (triggered)"
            axes[i].annotate(label, xycoords='axes fraction', xy=(x, y), zorder=1000)

            # energy range label
            emin = phaii.ebounds[channel_range[0]].emin
            emax = phaii.ebounds[channel_range[1]].emax
            axes[i].annotate('{0:3.0f}-{1:3.0f} keV'.format(emin, emax), xy=(1-x,y),
                             xycoords='axes fraction', horizontalalignment='right',
                             zorder=1000)

        if legend:
            handles = [
                matplotlib.patches.Rectangle((0,0), 1, 1, color='r', alpha=0.1),
                matplotlib.lines.Line2D([0, 1], [0, 1], color='r'),
                matplotlib.patches.Rectangle((0,0), 1, 1, color='b', alpha=0.2),
                matplotlib.lines.Line2D([0, 1], [0, 1], color='b')]
            labels = [
                "Background Window",
                "Background Fit",
                "Trigger Window",
                "Trigger Threshold"]
            legend_ax = fig.add_axes((0, 0, 1, bottom))
            legend_ax.set_axis_off()
            legend_ax.legend(handles, labels,
                loc=(0.15, 0.0), ncols=4, frameon=False, fontsize=12)

        return fig, axes, lcplots
