# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 425 3rd Street SW, Suite 950, Washington, DC 20024
#
# Copyright 2017-2024 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland, Adam Goldstein and Oliver J. Roberts
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


class Duration:
    def __init__(self, timebins_list, bkgds_list, duration_interval):
        # intializes with a list of TimeBins and BackgroundRates
        # this only checks that the inputs are valid
        self.timebins_list = timebins_list
        self.bkgds_list = bkgds_list
        self.duration_interval = duration_interval

    def quantiles(self, tparam, confidence):
        loconf = np.quantile(tparam, ((1 + confidence) / 2))
        uppconf = np.quantile(tparam, (1 - ((1 + confidence) / 2)))
        return loconf, uppconf

    def findtparams(self, array, dur_per, tbins):
        paramloc = np.where(array <= dur_per * np.max(array))
        param = np.max(tbins[paramloc[0]])
        return param

    def findtparams_err(self, num_sims, array, dur_per, tbins):
        list_lo = []
        for ii in range(num_sims):
            param2 = self.findtparams(array[ii], dur_per, tbins)
            list_lo.append(param2)
        return list_lo

    def error_prop(self, a, b):
        # error propagation, multiplication.
        c = np.sqrt((a ** 2) + (b ** 2))
        return c
        
    def calculate(self, num_sims, confidence):
        # Performs the calculation. The user must specify the duration interval they
        # want to calculate (e.g. (0.05, 0.95)). The user can also define the number
        # of sims and confidence region for the uncertainty calculation. This
        # will return a 3-tuple: (value, - error, + error)

        list_src_cts = []
        list_src_centroids = []
        list_src_cts_err = []
        for xx in range(len(self.timebins_list)):
            s_rates = self.timebins_list[xx]
            s_rates_cts = s_rates.counts
            s_rates_ct_err = s_rates.count_uncertainty
            s_rates_centroids = s_rates.centroids
            list_src_cts.append(s_rates_cts)
            list_src_centroids.append(s_rates_centroids)
            list_src_cts_err.append(s_rates_ct_err)
        # break
        timebins_y = np.sum(list_src_cts, axis=0)
        timebins_x = list_src_centroids[0]

        list_bkg_cts = []
        list_bkg_cts_err = []
        for x in range(len(self.bkgds_list)):
            b_rates = self.bkgds_list[x]
            b_rates_cts = b_rates.counts
            b_rates_cts_err = b_rates.count_uncertainty
            list_bkg_cts.append(b_rates_cts)
            list_bkg_cts_err.append(b_rates_cts_err)
            # break
        br = np.sum(list_bkg_cts, axis=0)

        workingrate = timebins_y.T[:, ] - br.T[:, ]
        cumflsum = workingrate.cumsum(axis=1)
        plotter_full = (timebins_x, cumflsum.T[:, ])

        t_lower = self.findtparams(plotter_full[1], self.duration_interval[0], timebins_x)
        t_higher = self.findtparams(plotter_full[1], self.duration_interval[1], timebins_x)
        t_diff = (t_higher - t_lower)

        # Errors  #################################################################################
        p_source_err_list = np.random.poisson(lam=(np.abs(timebins_y)), size=(num_sims, len(timebins_x)))

        err_prop = []
        for i, name in enumerate(list_bkg_cts_err):
            br_ind_err = (list_bkg_cts_err[i] ** 2)
            err_prop.append(br_ind_err)
        br_ind_errs = np.sum(err_prop, axis=0)
        br_err = np.sqrt(br_ind_errs)

        mu11, sigma11 = br[:, 0], br_err  # mean and standard deviation
        p_bkd_err_list = np.random.normal(mu11, sigma11, size=(num_sims, len(timebins_x)))

        diff_rate = p_source_err_list - p_bkd_err_list
        cuflux = diff_rate.cumsum(axis=1)
        arr = cuflux.T[:, ]

        arr2 = arr.T
        timebins2 = timebins_x.T

        f_err_lower = self.findtparams_err(num_sims, arr2, self.duration_interval[0], timebins2)
        f_err_higher = self.findtparams_err(num_sims, arr2, self.duration_interval[1], timebins2)

        #  Final Numbers #################################################################################

        tdiff_upp_lo_err = t_higher - self.quantiles(f_err_higher, confidence)[0]
        tdiff_upp_hi_err = t_higher - self.quantiles(f_err_higher, confidence)[1]

        tdiff_low_lo_err = t_lower - self.quantiles(f_err_lower, confidence)[0]
        tdiff_low_hi_err = t_lower - self.quantiles(f_err_lower, confidence)[1]

        tdiff_err_lo = self.error_prop(tdiff_upp_lo_err, tdiff_low_lo_err) * -1
        tdiff_err_hi = self.error_prop(tdiff_upp_hi_err, tdiff_low_hi_err)

        return t_diff, tdiff_err_lo, tdiff_err_hi
