#
#     Authors: William Cleveland (USRA),
#              Adam Goldstein (USRA) and
#              Daniel Kocevski (NASA)
#
#     Portions of the code are Copyright 2020 William Cleveland and
#     Adam Goldstein, Universities Space Research Association
#     All rights reserved.
#
#     Written for the Fermi Gamma-ray Burst Monitor (Fermi-GBM)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import os
import matplotlib.pyplot as plt
import astropy.units as u

from gdt.core.detector import Detectors

class ImageFileMixin:
    this_dir = os.path.dirname(__file__)
    image_file = os.path.join(this_dir, "test.png")

    def tearDown(self):
        plt.close('all')
        try:
            os.remove(self.image_file)
        except:
            pass


class ExampleDetectors(Detectors):
    det0 = ('Det0', 0,  45.0 * u.deg,  45.0 * u.deg)
    det1 = ('Det1', 1, 270.0 * u.deg, 135.0 * u.deg)
