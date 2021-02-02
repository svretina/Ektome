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

min_r = 0.5 # Schwarzshild radius for M = 1 in isotropic coords
proj_path = os.getcwd()
simulations_path = f"{proj_path}/simulations"
results_path = f"{proj_path}/results"
figures_path = f"{proj_path}/results/figures"
parfiles_path = f"{proj_path}/parfiles"
subfiles_path = f"{proj_path}/subfiles"
metadata_path = f"{proj_path}/metadata.txt"
config_path = f"{proj_path}/submit.ini"
