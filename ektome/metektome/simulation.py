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

"""This module provides a class to represent and load simulation data."""

import logging
import re
from pathlib import Path

import numpy as np
from kuibit.simdir import SimDir

import ektome.globals as glb
from ektome.exceptions import SimulationNotFoundError

# Configure logging
logger = logging.getLogger(__name__)


class Simulation:
    """Represents a simulation, providing methods to load parameters and data.

    Attributes:
        sim_name: Name of the simulation.
        dim: Dimensionality of the data (2 or 3).
        sim_dir: Path to the simulation data.
        param_file: Path to the parameter file.
        metadata_file: Path to the metadata file (TwoPunctures.bbh).
    """

    def __init__(self, sim_name: str, dim: int = 3) -> None:
        """Initializes the Simulation object.

        Args:
            sim_name: Name of the simulation folder.
            dim: Dimension of the simulation data.

        Raises:
            SimulationNotFoundError: If the simulation directory or metadata is missing.
        """
        self.sim_name = sim_name
        self.dim = dim
        self.sim_dir = glb.SIMULATIONS_PATH / self.sim_name
        self.param_file = glb.PARFILES_PATH / f"{self.sim_name}.par"
        self.metadata_file = self.sim_dir / "TwoPunctures.bbh"

        if not self.sim_dir.exists():
            raise SimulationNotFoundError(f"Simulation directory not found: {self.sim_dir}")
        if not self.metadata_file.exists():
            # Fallback check if it's in a different location or named differently
            logger.warning(f"Metadata file {self.metadata_file} not found.")

        self._read_params()
        self._calculate_norms()
        
        if self.dim == 2:
            self._load_data_2d()
        else:
            self._load_data_3d()

    def _get_param_value(self, param_name: str, file_path: Path) -> float:
        """Parses a numerical value for a parameter from a file.

        Args:
            param_name: The name of the parameter to look for.
            file_path: The file to search in.

        Returns:
            The parsed float value.

        Raises:
            ValueError: If the parameter is not found or cannot be parsed.
        """
        if not file_path.exists():
            raise FileNotFoundError(f"File {file_path} not found while reading {param_name}")

        pattern = re.compile(rf"{param_name}\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")
        
        with open(file_path) as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    return float(match.group(1))
        
        raise ValueError(f"Parameter '{param_name}' not found in {file_path}")

    def _calculate_norms(self) -> None:
        """Calculates spin norms for both punctures."""
        self.s1 = np.sqrt(self.s1x**2 + self.s1y**2 + self.s1z**2)
        self.s2 = np.sqrt(self.s2x**2 + self.s2y**2 + self.s2z**2)

    def _read_params(self) -> None:
        """Reads simulation parameters from metadata and parameter files."""
        try:
            # ADM masses and separation from metadata
            self.mm = self._get_param_value("adm-mass2", self.metadata_file)
            self.mp = self._get_param_value("adm-mass1", self.metadata_file)
            self.par_b = self._get_param_value("separation", self.metadata_file)
            self.ex_r = self._get_param_value("excision", self.metadata_file)

            # Grid spacing from parameter file
            self.dx = self._get_param_value("dx", self.param_file)
            self.dy = self._get_param_value("dy", self.param_file)
            self.dz = self._get_param_value("dz", self.param_file)

            # Spins from metadata
            # Note: Mapping might be reversed as per legacy comments
            self.s1x = self._get_param_value("spin2x", self.metadata_file)
            self.s1y = self._get_param_value("spin2y", self.metadata_file)
            self.s1z = self._get_param_value("spin2z", self.metadata_file)
            self.s2x = self._get_param_value("spin1x", self.metadata_file)
            self.s2y = self._get_param_value("spin1y", self.metadata_file)
            self.s2z = self._get_param_value("spin1z", self.metadata_file)
        except (ValueError, FileNotFoundError) as exc:
            logger.error(f"Error reading parameters for {self.sim_name}: {exc}")
            # Assign defaults or re-raise
            raise SimulationNotFoundError(f"Missing required parameters in {self.sim_name}") from exc

    def _load_data_3d(self) -> None:
        """Loads 3D gravitational field data using kuibit."""
        logger.info(f"Loading 3D data for {self.sim_name}...")
        sim = SimDir(str(self.sim_dir))
        self.u = sim.gf.xyz.fields.puncture_u[0]
        self.psi = sim.gf.xyz.fields.my_psi[0]

    def _load_data_2d(self) -> None:
        """Loads 2D gravitational field data using kuibit."""
        logger.info(f"Loading 2D data for {self.sim_name}...")
        sim = SimDir(str(self.sim_dir))
        self.u = sim.gf.xy.fields.puncture_u[0]
        self.psi = sim.gf.xy.fields.my_psi[0]

    def _sphere_dist_sq(self, x, y, z) -> float:
        """Calculates distance squared from the first puncture."""
        return (x + self.par_b)**2 + y**2 + z**2

    def _circle_dist_sq(self, x, y) -> float:
        """Calculates distance squared from the first puncture in 2D."""
        return (x + self.par_b)**2 + y**2
