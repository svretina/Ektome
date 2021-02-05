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

import configparser as cfg
import os
import sys
import numpy as np
import ektome.proektome.submit as sb
import ektome.globals as glb
import math


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
    temp = {}
    for section in parser:
        if section == "DEFAULT":
            continue
        for key in parser[section]:
            temp[key] = parser[section][key]
        config[section] = temp
        temp = {}
    return config


def create_dirs():
    """Checks if output directories already exists and if not,
    it creates them.
    """

    if not os.path.exists(f"{glb.proj_path}/simulations"):
        os.mkdir(f"{glb.proj_path}/simulations")
    else:
        pass
    # if not os.path.exists(f"{glb.proj_path}/subfiles"):
    #     os.mkdir(f"{glb.proj_path}/subfiles")
    # else:
    #     pass
    if not os.path.exists(f"{glb.proj_path}/parfiles"):
        os.mkdir(f"{glb.proj_path}/parfiles")
    else:
        pass
    if not os.path.exists(f"{glb.proj_path}/results"):
        os.mkdir(project_path + "/results")
    else:
        pass
    if not os.path.exists(f"{glb.proj_path}/results/figures"):
        os.mkdir(f"{glb.proj_path}/results/figures")
    else:
        pass
    if not os.path.exists(f"{glb.proj_path}/results/figures/u"):
        os.mkdir(f"{glb.proj_path}/results/figures/u")
    else:
        pass
    if not os.path.exists(f"{glb.proj_path}/results/figures/psi"):
        os.mkdir(f"{glb.proj_path}/results/figures/psi")
    else:
        pass
    if not os.path.exists(f"{glb.proj_path}/results/figures/error"):
        os.mkdir(f"{glb.proj_path}/results/figures/error")
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
            if "array" in cfg_dict[section].keys():
                temp_array = cfg_dict[section]["array"].split(",")
                config_array[section] = np.array([float(x) for x in temp_array])
            else:
                start = float(cfg_dict[section]["start"])
                end = float(cfg_dict[section]["end"])
                step = float(cfg_dict[section]["step"])
                N = int((end - start) / step + 1)
                tmp = np.linspace(start, end, N)
                config_array[section] = np.round(tmp, 3)

        except (ValueError, KeyError):
            config_array[section] = np.nan
    return config_array


def submit_simulation(simulation_dict):
    _ = sb.submit(simulation_dict)


def create_simulation_dict_and_submit(cfg_arr):
    """Creates a dictionary with all the parameters for the
    simulation.

    :param cfg_arr: Dictionary with the config info.
    :type cfg_arr: dict.
    :returns: Dictionary with simulation info
    """
    # Can be replaced with cartesian product of arrays
    # and loop over it.
    simulation = {}
    counter = 0
    for q in cfg_arr["mass_ratio"]:
        for b in cfg_arr["par_b"]:
            for sx1 in cfg_arr["plus_spin_x"]:
                for sy1 in cfg_arr["plus_spin_y"]:
                    for sz1 in cfg_arr["plus_spin_z"]:
                        for pz1 in cfg_arr["plus_momentum_z"]:
                            for py1 in cfg_arr["plus_momentum_y"]:
                                for px1 in cfg_arr["plus_momentum_x"]:
                                    for sx2 in cfg_arr["minus_spin_x"]:
                                        for sy2 in cfg_arr["minus_spin_y"]:
                                            for sz2 in cfg_arr["minus_spin_z"]:
                                                for pz2 in cfg_arr["minus_momentum_z"]:
                                                    for py2 in cfg_arr[
                                                        "minus_momentum_y"
                                                    ]:
                                                        for px2 in cfg_arr[
                                                            "minus_momentum_x"
                                                        ]:
                                                            max_r = q / 2.0
                                                            if not math.isnan(
                                                                cfg_arr["excision"]
                                                            ):
                                                                step = float(
                                                                    cfg["excision"][
                                                                        "step"
                                                                    ]
                                                                )
                                                                radii = np.arange(
                                                                    glb.min_r,
                                                                    max_r + step,
                                                                    step,
                                                                )
                                                            else:
                                                                radii = [max_r]
                                                            for ex_r in radii:
                                                                if math.isnan(
                                                                    cfg_arr["error"]
                                                                ):
                                                                    simulation["q"] = q
                                                                    simulation[
                                                                        "par_b"
                                                                    ] = b
                                                                    simulation[
                                                                        "px1"
                                                                    ] = px1
                                                                    simulation[
                                                                        "py1"
                                                                    ] = py1
                                                                    simulation[
                                                                        "pz1"
                                                                    ] = pz1
                                                                    simulation[
                                                                        "sx1"
                                                                    ] = sx1
                                                                    simulation[
                                                                        "sy1"
                                                                    ] = sy1
                                                                    simulation[
                                                                        "sz1"
                                                                    ] = sz1
                                                                    simulation[
                                                                        "px2"
                                                                    ] = px2
                                                                    simulation[
                                                                        "py2"
                                                                    ] = py2
                                                                    simulation[
                                                                        "pz2"
                                                                    ] = pz2
                                                                    simulation[
                                                                        "sx2"
                                                                    ] = sx2
                                                                    simulation[
                                                                        "sy2"
                                                                    ] = sy2
                                                                    simulation[
                                                                        "sz2"
                                                                    ] = sz2

                                                                    simulation[
                                                                        "ex_r"
                                                                    ] = ex_r
                                                                    submit_simulation(
                                                                        simulation
                                                                    )
                                                                    sys.exit()
                                                                    counter = (
                                                                        counter + 1
                                                                    )


if __name__ == "__main__":
    # Read the config file
    config_dict = read_config()
    config_arr = create_config_arrays(config_dict)
    # Create output directories
    create_dirs()
    # Create simulation dictionary and submit it
    create_simulation_dict_and_submit(config_arr)


# Support for error arays ( run for different error tollerances)
# elif (( $error > $q/(4.0*$b - $q) ));then
# echo "Error too large"
# echo "Defaulting to EH of m+"
# err=$(echo "scale=3;$q/(4.0*$b - $q) "| bc )
# ./submit_bone.sh -q $q -b $b -px $px -py $py -e $err -exr $ex_r
# echo "========================="
# elif (( $error < $q/(4.0*$b - 1.0) ));then
# echo "Error too small"
# echo "Defaulting to EH of m-"
# err=$(echo "scale=3;$q/(4.0*$b - 1.0) "| bc )
# ./submit_bone.sh -q $q -b $b -px $px -py $py -e $err -exr $ex_r
# echo "========================="
# else
#     bash submit_bone.sh $q $b $p $error
# echo "error===================="
