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

"""This module provides helper functions to manage submition of simulations.
The functions available are:
- :py:func:`~.read_config`: Reads a config file and outputs a dictionary with the values of the config file.
- :py:func:`~.create_dirs`: Checks if output directories exists and if not, creates them
- :py:func:`~.get_ini_file`: Searches the directory for a .ini config file and returns its path.

"""

import os
import sys
import math
import numpy as np
import configparser as cfg
import ektome.globals as glb
import itertools
import ektome.proektome.submit as sb
import ektome.proektome.binary as bnr
import ektome.proektome.simulation as simulation

def read_config():
    """Reads a config file and returns a dictionary with its
    contents.

    :returns: Dictionary with the values of the config file
              variable.
    :rtype: dict
    """

    parser = cfg.ConfigParser()
    try:
        parser.read(glb.config_path)
    except ValueError:
        parser.read(get_ini_file())
    config = {}
    for section in parser:
        temp = {}
        if section == "DEFAULT":
            continue
        for key in parser[section]:
            temp[key] = parser[section][key]
        config[section] = temp
    return config


def create_dirs():
    """Checks if output directories already exists and if not,
    it creates them.
    """

    if not os.path.exists(f"{glb.simulations_path}"):
        os.mkdir(f"{glb.simulations_path}")
    else:
        pass
    if not os.path.exists(f"{glb.parfiles_path}"):
        os.mkdir(f"{glb.parfiles_path}")
    else:
        pass
    if not os.path.exists(f"{glb.results_path}"):
        os.mkdir(f"{glb.results_path}")
    else:
        pass
    if not os.path.exists(f"{glb.figures_path}"):
        os.mkdir(f"{glb.figures_path}")
    else:
        pass
    if not os.path.exists(f"{glb.figures_path}/u"):
        os.mkdir(f"{glb.figures_path}/u")
    else:
        pass
    if not os.path.exists(f"{glb.figures_path}/psi"):
        os.mkdir(f"{glb.figures_path}/psi")
    else:
        pass
    if not os.path.exists(f"{glb.figures_path}/error"):
        os.mkdir(f"{glb.figures_path}/error")
    else:
        pass


def get_ini_file():
    """Searches for a .ini config file and returns its path.

    :returns: Path of the .ini config file
    :rtype: str
    """

    for file in os.listdir(glb.proj_path):
        if file.endswith(".ini"):
            return file
        continue
    return None


def clear_submit_metadata():
    if os.path.exists(glb.metadata_path):
        os.remove(glb.metadata_path)

def my_arange(start,step,end):
    N = int((end - start) / step + 1)
    tmp = np.linspace(start, end, N)
    return np.round(tmp, 3)


def create_config_arrays(cfg_dict):
    """Creates numpy arrays spaced according to the steps
    provided in the config file.

    :param cfg_dict: Dictionary containing the config information
    :type cfg_dict: dict
    :returns: Dictionary with the numpy arrays
    :rtype: dict
    """
    config_array = {}
    for section in cfg_dict:
        try:
            temp_array = cfg_dict[section]["array"].split(",")
            content = []
            for i in temp_array:
                if ":" in i:
                    start, step, end = [float(j) for j in i.split(":")]
                    array = my_arange(start,step,end)
                    content.append(array)
                else:
                    content.append(np.array([float(i)]))
            config_array[section] = np.sort(np.concatenate(content))


        except (ValueError, KeyError):
            config_array[section] = np.nan
    return config_array


def create_simulation_dict_and_submit(cfg_arr):
    """Creates a dictionary with all the parameters for the
    simulation.

    :param cfg_arr: Dictionary with the config info.
    :type cfg_arr: dict.
    :returns: Dictionary with simulation info
    """
    spin1 = list(itertools.product(cfg_arr["minus_spin_x"],
                                   cfg_arr["minus_spin_y"],
                                   cfg_arr["minus_spin_z"]))

    spin2 = list(itertools.product(cfg_arr["plus_spin_x"],
                                   cfg_arr["plus_spin_y"],
                                   cfg_arr["plus_spin_z"]))

    for q in cfg_arr["mass_ratio"]:
        exr = q / 2.0
        binary = bnr.Binary(q)
        for s1 in spin1:
            if np.linalg.norm(s1) >= 1:
                print("|Spin| > 1")
                continue
            for s2 in spin2:
                if np.linalg.norm(s2) >= 1:
                    print("|Spin| > 1")
                    continue
                for n_orb in cfg_arr["number_of_orbits"]:
                    b = binary.semimajor(n_orb)
                    p1, p2 = binary.quasicircular_inspiral(q,
                                                           2*b,
                                                           s1,
                                                           s2)
                    sim = simulation.Simulation(q=q,b=b,
                                                px1=p1[0],
                                                py1=p1[1],
                                                pz1=p1[2],
                                                sx1=s1[0],
                                                sy1=s1[1],
                                                sz1=s1[2],
                                                px2=p2[0],
                                                py2=p2[1],
                                                pz2=p2[2],
                                                sx2=s2[0],
                                                sy2=s2[1],
                                                sz2=s2[2],
                                                exr=exr)
                    sb.submit(sim)


if __name__ == "__main__":
    # Read the config file
    config_dict = read_config()
    config_arr = create_config_arrays(config_dict)
    # Create output directories
    create_dirs()
    # Create simulation dictionary and submit it
    create_simulation_dict_and_submit(config_arr)
