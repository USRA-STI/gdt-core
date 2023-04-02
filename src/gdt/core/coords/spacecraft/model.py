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
from abc import ABC, abstractmethod
from astropy.timeseries import TimeSeries
from .frame import SpacecraftFrame


class SpacecraftFrameModelMixin(ABC):  # pragma: no cover - Simple mixin with not testable code.
    """Mixin that prototypes the retrieval and saving of spacecraft frame information."""

    @abstractmethod
    def get_spacecraft_frame(self) -> SpacecraftFrame:
        """Retrieves the spacecraft frame from the model.

        Returns:
            (:class:`~gdt.core.coords.SpacecraftFrame`)
        """
        pass

    def set_spacecraft_frame(self, frame: SpacecraftFrame):
        """Saves the spacecraft frame to the model.

        Args:
            frame (:class:`~gdt.core.coords.SpacecraftFrame`): 
                The object to be saved to the model.
        """
        raise NotImplementedError('Model is read-only')


class SpacecraftStatesModelMixin(ABC):  # pragma: no cover - Simple mixin with no testable code
    """Mixin that prototypes the retrieval and saving of spacecraft state information."""

    def get_spacecraft_states(self) -> TimeSeries:
        """Retrieves the spacecraft state information from the model.

        Returns:
            (astropy.timeseries.TimeSeries)
        """
        pass

    def set_spacecraft_states(self, series: TimeSeries):
        """Saves the spacecraft state to the model.

        Args:
            series (astropy.timeseries.TimeSeries): Object containing spacecraft 
                                                    states to be saved to model.
        """
        raise NotImplementedError('Model is read-only')
