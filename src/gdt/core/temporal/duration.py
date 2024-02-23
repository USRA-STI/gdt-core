# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 425 3rd Street SW, Suite 950, Washington, DC 20024
#
# Copyright 2017-2024 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland, Adam Goldstein and Oliver Roberts
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
        # intializes with a list of TimeBins, BackgroundRates and specified duration interval.
        # The user must specify the duration interval they want to calculate (e.g. (0.05, 0.95)).
        # This only checks that the inputs are valid
        self.timebins_list = timebins_list
        self.bkgds_list = bkgds_list
        self.duration_interval = duration_interval

    def calculate(self, num_sims, confidence):
        # Performs the calculation. The user can also define the number
        # of sims and confidence region for the uncertainty calculation. This
        # will return a 3-tuple: (value, - error, + error)

        ####################  Calculation  ###########################

        list_src_cts = []
        list_src_centroids = []
        for xx in range(len(self.timebins_list)):
            s_rates = self.timebins_list[xx]
            s_rates_cts = s_rates.counts
            s_rates_centroids = s_rates.centroids
            list_src_cts.append(s_rates_cts)
            list_src_centroids.append(s_rates_centroids)
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
        # plt.plot(plotter_full[0],plotter_full[1])

        f_lower = np.where(plotter_full[1] <= self.duration_interval[0] * np.max(plotter_full[1]))
        t_lower = np.max(timebins_x[f_lower[0]])
        # t_lower = round(t_lower2, 2)

        f_higher = np.where(plotter_full[1] <= self.duration_interval[1] * np.max(plotter_full[1]))
        t_higher = np.max(timebins_x[f_higher[0]])
        # t_higher = round(t_higher2, 2)

        t_diff = (t_higher - t_lower)
        # t_diff = round(f_diff, 4)

        ####################  Errors  ###########################

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
        timbins2 = timebins_x.T

        f_err_lower = []
        for i in range(num_sims):
            h_err_lower = np.where(arr2[i] <= self.duration_interval[0] * np.max(arr2[i]))
            g_err_lower = np.max(timbins2[h_err_lower[0]])
        f_err_lower.append(g_err_lower)

        f_err_higher = []
        for i in range(num_sims):
            h_err_higher = np.where(arr2[i] <= self.duration_interval[1] * np.max(arr2[i]))
            g_err_higher = np.max(timbins2[h_err_higher[0]])
            f_err_higher.append(g_err_higher)

        ####################  Final Numbers  ###########################

        tupp_err_loconf = np.quantile(f_err_higher,((1+confidence)/2))
        tupp_err_hiconf = np.quantile(f_err_higher,(1-((1+confidence)/2)))

        tlow_err_loconf = np.quantile(f_err_lower,((1+confidence)/2))
        tlow_err_hiconf = np.quantile(f_err_lower,(1-((1+confidence)/2)))

        tdiff_upp_lo_err = t_higher - tupp_err_loconf
        tdiff_upp_hi_err = t_higher - tupp_err_hiconf

        tdiff_low_lo_err = t_lower - tlow_err_loconf
        tdiff_low_hi_err = t_lower - tlow_err_hiconf

        fdiff_err_lo = np.sqrt((tdiff_upp_lo_err ** 2) + (tdiff_low_lo_err ** 2))
        fdiff_err_hi = np.sqrt((tdiff_upp_hi_err ** 2) + (tdiff_low_hi_err ** 2))
        tdiff_err_lo = fdiff_err_lo * -1
        tdiff_err_hi = fdiff_err_hi

        return t_diff, tdiff_err_lo, tdiff_err_hi