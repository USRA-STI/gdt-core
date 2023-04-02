#  CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
#  Contract No.: CA 80MSFC17M0022
#  Contractor Name: Universities Space Research Association
#  Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
#  Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
#  Developed by: William Cleveland and Adam Goldstein
#                Universities Space Research Association
#                Science and Technology Institute
#                https://sti.usra.edu
#
#  Developed by: Daniel Kocevski
#                National Aeronautics and Space Administration (NASA)
#                Marshall Space Flight Center
#                Astrophysics Branch (ST-12)
#
#   Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
#   in compliance with the License. You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software distributed under the License
#  is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
#  implied. See the License for the specific language governing permissions and limitations under the
#  License.
#
import datetime as dt
from typing import Union

import numpy as np

from astropy.time import Time, TimeDelta
import astropy.units as u


def truncate_hour(time: Time):
    """Truncates the given time object to the current hour.

    Args:
        time (astropy.time.Time): Time object that will be truncated.

    Returns:
        (astropy.time.Time)
    """
    d = time.utc.to_datetime()
    return Time(dt.datetime(d.year, d.month, d.day, d.hour), format='datetime', scale='utc')


def time_range(start: Time, stop: Union[Time, TimeDelta, u.Quantity, int, float],
               step: Union[TimeDelta, u.Quantity, int, float] = 1.0) -> Time:
    """Creates an array of Time values within the specified range.
    The values will be >= start times and < then the stop time (given or 
    calculated).

    Args:
        start (astropy.time.Time): Time value at the beginning of the range
        stop (astropy.time.Time, astropy.units.Quantity, int, or float): 
            Either a time value signifying the end of the range, the duration of 
            the time range as a Quantity with appropriate unit, or an int or 
            float representing the number of seconds.
        step (datetime.TimeDelta, astropy.units.Quantity, int, or float, optional):
            Either a TimeDelta value, Quantity with appropriate unit, or a int 
            or float representing the number of seconds.

    Returns:
        (astropy.time.Time)
    """
    # Preserve the format and scale of the start value
    start_format = start.format
    start_scale = start.scale

    # Determine the stop value
    if isinstance(stop, Time):
        stop_value = stop.unix_tai
    elif isinstance(stop, TimeDelta):
        stop_value = (start + stop).unix_tai
    elif isinstance(stop, u.Quantity):
        stop_value = (start + TimeDelta(stop)).unix_tai
    else:
        stop_value = (start + TimeDelta(stop, format='sec')).unix_tai

    # Determine the increment value
    if isinstance(step, TimeDelta):
        increment = step.sec
    elif isinstance(step, u.Quantity):
        increment = step.to(u.s).value
    else:
        increment = step

    # create the array
    result = Time(np.arange(start.unix_tai, stop_value, increment), format='unix_tai', scale=start_scale)
    result.format = start_format
    return result


def hours_range_from(hours: int, *, from_time=Time.now()):
    """Returns a collection of hours starting from the given hours but not 
    including the current hour of from time.

    Args:
        hours (int): Number of hours from the given ``from_time``
        from_time (astropy.time.Time, optional): The time from which the hours 
                                                 range is based (defaults to 
                                                 current time).

    Returns:
        (list of astropy.time.Time)
    """
    hr = truncate_hour(from_time)
    return time_range(hr - TimeDelta(hours * 3600, format='sec'), hr, 3600)
