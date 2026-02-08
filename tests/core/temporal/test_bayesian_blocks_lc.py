import unittest

import numpy as np
from gdt.core.data_primitives import TimeBins
from gdt.core.temporal import BayesianBlocksLightcurve

from gdt.core.plot.lightcurve import Lightcurve

import matplotlib.pyplot as plt

nbins = 100
time_lo = np.arange(nbins, dtype=float)
time_hi = time_lo + 1.
time_mid = time_lo + .5
exposure = np.ones(nbins)
mean_counts = 100 + .1 * time_mid * time_mid  # Bkg
mean_counts[30:40] += 100
rng = np.random.default_rng(seed=123)  # Make test deterministic
counts = rng.poisson(mean_counts).astype(int)
lc = TimeBins(counts, time_lo, time_hi, exposure)

class TestBayesianBlocksLightcurve(unittest.TestCase):

    def test_compute_bayesian_blocks(self):

        bb_lc = BayesianBlocksLightcurve(lc)
        bb_lc.compute_bayesian_blocks()

        self.assertEqual(bb_lc.duration(1.), 9.0)
        self.assertEqual(bb_lc.signal_range.tstart, 30.0)
        self.assertEqual(bb_lc.signal_range.tstop, 40.0)

if __name__ == '__main__':
    unittest.main()
