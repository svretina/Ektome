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


"""This module provides helper functions to write submit files.
The functions available are:
- :py:func:`~.create_sub_file`: Writes a submition file.
"""

import ektome.globals as glb


def create_sub_dict(base_name):
    """Writes the submition file.

    :param simulation: Dictionary containing simulation info.
    :type simulation: dict
    :param: proj_path: Path of project directory.
    :type proj_path: str
    :param base_name: Name string that serves as a base to name parameter files.
    :type base_name: str
    :returns: A Dictionary with the submition info
    :rtype: dict
    """

    par_file_path = "/".join((glb.parfiles_path, base_name + ".par"))
    # How much memory to request
    # if simulation["par_b"] == 50:
    #     memory = "15000MB"
    # elif simulation["par_b"] == 55:
    #     memory = "20000MB"
    # elif simulation["par_b"] == 60:
    #     memory = "25000MB"
    # else:
    #     memory = "10000MB"

    sim_dir = "/".join((glb.simulations_path, base_name))
    sub_file_dict = {
        "executable": "/work/stamatis.vretinaris/opt/Cactus/exe/cactus_sim",
        "arguments": par_file_path,
        "universe": "vanilla",
        "output": sim_dir + "/out",
        "error": sim_dir + "/err",
        "log": sim_dir + "/log",
        "accounting_group": "cbc.prod.initial_data",
        "request_cpus": "Target.TotalSlotCPUs",
        "getenv": "True",
        "request_memory": 18000,
        "Rank": "Memory",
        "on_exit_hold": "( ExitCode != 0 )",
        "kill_sig": "15",
    }
    # "request_memory":f"max({{{memory}, Target.TotalSlotMemory}})",

    return sub_file_dict
