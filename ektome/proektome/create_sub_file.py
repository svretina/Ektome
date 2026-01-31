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

"""This module provides helper functions to write submit files."""

from typing import Any, Dict

import ektome.globals as glb


def create_sub_dict(base_name: str) -> Dict[str, Any]:
    """Creates a dictionary containing the HTCondor submission information.

    Args:
        base_name: Name string that serves as a base to name parameter files.

    Returns:
        A dictionary with the submission info.
    """
    par_file_path = glb.PARFILES_PATH / f"{base_name}.par"
    
    # How much memory to request
    try:
        b_val = int(base_name.split("_")[2].split("b")[1])
    except (IndexError, ValueError):
        b_val = 0

    if 100 < b_val <= 1000:
        memory = "15000MB"
    elif b_val > 1000:
        memory = "12000MB"
    else:
        memory = "8000MB"  # Default fallback

    sim_dir = glb.SIMULATIONS_PATH / base_name
    
    sub_file_dict: Dict[str, Any] = {
        "executable": "/work/stamatis.vretinaris/opt/Cactus/exe/cactus_sim",
        "arguments": str(par_file_path),
        "universe": "vanilla",
        "output": str(sim_dir / "out"),
        "error": str(sim_dir / "err"),
        "log": str(sim_dir / "log"),
        "accounting_group": "cbc.prod.initial_data",
        "request_cpus": 5,
        "getenv": "True",
        "request_memory": memory,
        "Rank": "Memory",
        "on_exit_hold": "( ExitCode != 0 )",
        "kill_sig": "15",
    }

    return sub_file_dict
