#!/usr/bin/env python3

# Copyright (C) 2021 Stamatis Vretinaris
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <https://www.gnu.org/licenses/>.

"""This module provides functions to calculate errors between simulations."""

import logging

import numpy as np
import pandas as pd

import ektome.metektome.simulation as sim
from ektome.exceptions import SimulationNotFoundError

# Configure logging
logger = logging.getLogger(__name__)


class Error:
    """Calculates and reports errors between a vanilla and an excision simulation.

    Attributes:
        vnl_sim_name: Name of the vanilla simulation.
        exc_sim_name: Name of the corresponding excision simulation.
        dim: Dimensionality of the data.
    """

    def __init__(self, vnl_sim_name: str, dim: int = 3) -> None:
        """Initializes the Error object by loading both simulations.

        Args:
            vnl_sim_name: Name of the vanilla simulation folder.
            dim: Dimension of the data.

        Raises:
            SimulationNotFoundError: If either simulation cannot be loaded.
        """
        self.vnl_sim_name = vnl_sim_name
        # Assuming naming convention: vanilla_... -> excision_...
        suffix = vnl_sim_name.split("vanilla")[-1]
        self.exc_sim_name = f"excision{suffix}"
        self.dim = dim

        try:
            self.vanilla = sim.Simulation(self.vnl_sim_name, self.dim)
            self.excision = sim.Simulation(self.exc_sim_name, self.dim)
        except SimulationNotFoundError as exc:
            logger.error(f"Could not load simulations for error comparison: {exc}")
            raise

        self.q = self.vanilla.mp / self.vanilla.mm
        self.ex_r = self.excision.ex_r
        self.par_b = self.vanilla.par_b

        self._calculate_error_psi()
        self._calculate_error_psi_theoretical()

        self.psimax = self.calculate_max_with_mask(self.error_psi)
        self.psimaxt = self.calculate_max_with_mask(self.error_psi_t)

    def _calculate_error_psi(self) -> None:
        """Calculates the relative error in the psi field."""
        self.error_psi = abs(self.vanilla.psi - self.excision.psi) / self.vanilla.psi

    def _calculate_error_psi_theoretical(self) -> None:
        """Calculates a theoretical error estimate for the psi field."""
        temp = self.q / (4.0 * self.vanilla.par_b - self.q)
        self.error_psi_t = temp / self.vanilla.psi

    def calculate_max_with_mask(self, var) -> float:
        """Calculates the maximum value of a field within a masked region.

        The mask excludes the region inside the excision radius.

        Args:
            var: The field data (typically from kuibit).

        Returns:
            The maximum value found.
        """
        max_vals = []
        for _ref_level, _comp_index, unif_grid in var:
            coords = unif_grid.coordinates_from_grid()
            data = unif_grid.data
            
            # Create mask: True where we WANT to keep data
            if self.dim == 3:
                x, y, z = coords
                # Distance from puncture at (-par_b, 0, 0)
                dist_sq = (x + self.par_b)**2 + y**2 + z**2
            else:
                x, y = coords
                dist_sq = (x + self.par_b)**2 + y**2

            mask = dist_sq >= (self.ex_r**2)
            
            # Apply mask to data
            masked_data = data[mask]
            
            if masked_data.size > 0:
                max_vals.append(np.nanmax(masked_data))

        if not max_vals:
            logger.warning(f"No valid data points found after masking for {self.vnl_sim_name}")
            return 0.0

        return float(np.max(max_vals))

    def error_report(self) -> pd.DataFrame:
        """Generates a summary report of the errors and simulation parameters.

        Returns:
            A pandas DataFrame with one row of summary data.
        """
        report = {
            "q": self.q,
            "b": self.par_b,
            "ex_r": self.ex_r,
            "s1x": self.vanilla.s1x,
            "s1y": self.vanilla.s1y,
            "s1z": self.vanilla.s1z,
            "s1": self.vanilla.s1,
            "s2x": self.vanilla.s2x,
            "s2y": self.vanilla.s2y,
            "s2z": self.vanilla.s2z,
            "s2": self.vanilla.s2,
            "max_error_psi": self.psimax,
            "max_error_psi_theoretical": self.psimaxt,
        }
        return pd.DataFrame([report])
