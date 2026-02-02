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

"""This module provides global variables which will be used throughout the package."""

import os
from pathlib import Path

# Schwarzshild radius for M = 1 in isotropic coords
MIN_R: float = 0.5

# Base project path - allow override via Environment Variable
PROJECT_PATH: Path = Path(os.getenv("EKTOME_PROJECT_PATH", Path.cwd())).resolve()

# Derived paths using pathlib
SIMULATIONS_PATH: Path = PROJECT_PATH / "simulations"
RESULTS_PATH: Path = PROJECT_PATH / "results"
FIGURES_PATH: Path = RESULTS_PATH / "figures"
PARFILES_PATH: Path = PROJECT_PATH / "parfiles"
SUBFILES_PATH: Path = PROJECT_PATH / "subfiles"
METADATA_PATH: Path = PROJECT_PATH / "metadata.txt"
CONFIG_PATH: Path = PROJECT_PATH / "submit.ini"

# Ensure essential directories exist (this will be used by spine.py)
REQUIRED_DIRECTORIES = [
    SIMULATIONS_PATH,
    RESULTS_PATH,
    FIGURES_PATH,
    FIGURES_PATH / "u",
    FIGURES_PATH / "psi",
    FIGURES_PATH / "error",
    PARFILES_PATH,
    SUBFILES_PATH,
]
