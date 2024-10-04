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
from astropy.io import fits


def _create_card(key: str, value: str, comment: str) -> fits.Card:
    """create card with a string image."""
    raw_str = f'{key:8s}= {value:20s}'
    if comment is not None:
        raw_str += f' / {comment}'
    return fits.Card.fromstring(raw_str)


def fixed_card(key: str, value: float, places: int = 5, comment: str = None) -> fits.Card:
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


def exponential_card(key: str, value: float, places: int = 5, use_double: bool = False,
                     comment: str = None) -> fits.Card:
    """Create a header card with a fixed decimal floating-point value.

    Args:
        key (str): key index value of the header card.
        value (float): The value being stored within the header card.
        places (int): The number of decimal places represented by the value.
        use_double (bool): Whether to use 'D' instead of 'E' to represent a double precision value.
        comment (str): The comment string to add to the header card.

    Returns:
        fits.Card: The header card with an exponential floating-point value.
    """
    fmt = fr'{{:.{places}e}}'
    mantissa, exponent = fmt.format(value).split('e')
    # converting the exponential string to integer to remove zero-padding
    exp_val = int(exponent)
    if use_double:
        value_str = f'{mantissa}D{exp_val}'
    else:
        value_str = f'{mantissa}E{exp_val}'
    return _create_card(key, value_str, comment)
