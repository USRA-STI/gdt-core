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
import math
from astropy.io.fits.card import Card

def _sci_notation(value: float):
    exponent = int(math.floor(math.log10(value)))
    mantissa = float(value / math.pow(10, exponent))
    return mantissa, exponent


def _create_card(key: str, value: str, comment: str) -> Card:
    """create card with a string image."""
    raw_str = f'{key:8s}= {value:>20s}'
    if comment is not None:
        raw_str += f' / {comment}'
    return Card.fromstring(raw_str)


def fixed_card(key: str, value: float, places: int = 5, comment: str = None) -> Card:
    """Create a header card with a fixed decimal floating-point value.

    Args:
        key (str): key index value of the header card.
        value (float): The value being stored within the header card.
        places (int): The number of decimal places represented by the value.
        comment (str): The comment string to add to the header card.

    Returns:
        fits.Card: The header card with a fixed floating-point value.
    """
    fmt = fr'{{:.{places}f}}'
    return _create_card(key, fmt.format(value), comment)


def scientific_card(key: str, value: float, places: int = 5, comment: str = None, *,
                    use_d: bool = False) -> Card:
    """Create a header card with a float-point value in scientific notation.

    Args:
        key (str): key index value of the header card.
        value (float): The value being stored within the header card.
        places (int): The number of decimal places represented by the value.
        comment (str): The comment string to add to the header card.
        use_d (bool): If true, use 'D' as seperator otherwise use 'E'. Default is False.
    Returns:
        Card: The header card with an exponential floating-point value.
    """
    val_str = f'{{:.{places}f}}{{:s}}{{:d}}'
    mantissa, exponent = _sci_notation(value)
    sep = 'D' if use_d else 'E'
    return _create_card(key, val_str.format(mantissa, sep, exponent), comment)
