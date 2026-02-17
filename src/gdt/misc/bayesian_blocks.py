from gdt.core.data_primitives import TimeBins

import numpy as np

# Copyright (c) 2011-2022, Astropy Developers
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice, this list of conditions
#    and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions
#    and the following disclaimer in the documentation and/or other materials provided with the
#    distribution.
#  * Neither the name of the Astropy Team nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# From Israel Martinez:
# Much of this code was copied from the bayesian blocks implementation
# in astropy https://github.com/astropy/astropy/blob/2db2f820eb51c95fcb3e187328c8cbd99ecd24df/astropy/stats/bayesian_blocks.py
# With the following differences:
# 1. The following issue was fixed https://github.com/astropy/astropy/issues/14017
# 2. It can handle bins with 0 events (Scargle's convention is to always have at least 1)
# 3. The code was simplified to handle the specific case of binned data
# 4. In returns change point indices, instead of edges.

def bayesian_blocks(lc, p0=0.05, gamma=None, ncp_prior=None):
    """
    Perform a binned bayesian blocks identification (Scargle, 2013),
    taking the exposure into account.

    Args:
        lc (TimeBins): Lighcurve
        p0 (float): False positive rate
        gamma (float): Alternatively, provide the gamma parameter 
            (slope in the prior). Take precedence over ``p0``
        ncp_prior (float): Specify the prior in the number of bins,
            taking precedence over ``p0`` and ``gamma``.

    Return:
        bb_indices (array): Bin indices of change points. A new lightcurve can be 
            computed as lc.rebin(rebin_by_edge_index, bb_indices)
    """
    
    # Compute prior
    if ncp_prior is None:
        if gamma is not None:
            ncp_prior = -np.log(gamma)
        elif p0 is not None:
            #Eq. 21 in Scargle (2013) (log missing)
            ncp_prior = 4 - np.log(73.53 * p0 * (lc.size**-0.478))
        else:
            raise RuntimeError("Specify either p0, gamma or ncp_prior")
        
    # ----------------------------------------------------------------
    # Start with first data cell; add one cell at each iteration
    # ----------------------------------------------------------------
    exposure_cumsum = np.append(np.cumsum(lc.exposure[::-1])[::-1], 0)
    counts_cumsum = np.append(np.cumsum(lc.counts[::-1])[::-1], 0)
    
    best = np.zeros(lc.size, dtype=float)
    last = np.zeros(lc.size, dtype=int)

    for R in range(lc.size):

        # evaluate fitness function. Eq. 19 from Scargle 2013
        T_k = exposure_cumsum[: R + 1] - exposure_cumsum[R + 1]
        N_k = counts_cumsum[: R + 1] - counts_cumsum[R+1] 

        # When N_k = 0, fit_vec is nan, but it should be 0
        fit_vec_log = np.zeros(N_k.size) #Prevent uninitialized values
        np.log(N_k / T_k, out = fit_vec_log, where = N_k != 0)
        
        fit_vec = N_k * fit_vec_log

        A_R = fit_vec - ncp_prior
        A_R[1:] += best[:R]

        i_max = np.argmax(A_R)
        last[R] = i_max
        best[R] = A_R[i_max]

    # ----------------------------------------------------------------
    # Now find changepoints by iteratively peeling off the last block
    # ----------------------------------------------------------------
    change_points = np.zeros(lc.size, dtype=int)
    i_cp = lc.size
    ind = lc.size
    while i_cp > 0:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    if i_cp == 0:
        change_points[i_cp] = 0
    change_points = change_points[i_cp:]

    return change_points
