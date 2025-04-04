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

"""This module provides functions to post process the simulation data."""

import os
import sys
import pandas as pd
from collections import Counter
from multiprocessing import Pool
import multiprocessing
import ektome.globals as glb
import ektome.metektome.error as e
import ektome.metektome.simulation as sim
import ektome.metektome.plotting_modules as pm
import time


def get_processed_sims(folders):
    processed_folders = []
    for folder in folders:
        sim_path = f"{glb.simulations_path}/{folder}"
        if os.path.exists(f"{sim_path}/psi.png"):
            processed_folders.append(folder)
    return processed_folders


def get_unfinished_sims(folders):
    unfinished_folders = []
    for folder in folders:
        sim_path = f"{glb.simulations_path}/{folder}"
        if not os.path.exists(f"{sim_path}/TwoPunctures.bbh"):
            unfinished_folders.append(folder)
    return unfinished_folders


def get_folders(recalc=False):
    print("Getting folders to process...")
    sim_folders = os.listdir(glb.simulations_path)
    vanilla_folders = [x for x in sim_folders if x.startswith("vanilla")]
    if recalc:
        unfinished = get_unfinished_sims(vanilla_folders)
        if not unfinished:
            exclude_folders = []
            folders = vanilla_folders
        else:
            exclude_folders = unfinished
            folders = list(
                (
                    Counter(vanilla_folders) - Counter(exclude_folders)
                ).elements()
            )
    else:
        processed = get_processed_sims(vanilla_folders)
        unfinished = get_unfinished_sims(vanilla_folders)
        if not processed and not unfinished:
            exclude_folders = []
            folders = vanilla_folders
        else:
            exclude_folders = get_processed_sims(
                vanilla_folders
            ) + get_unfinished_sims(vanilla_folders)

            folders = list(
                (
                    Counter(vanilla_folders) - Counter(exclude_folders)
                ).elements()
            )
    return folders


def get_error_dataframe(folder):
    print(folder)
    try:
        err = e.Error(folder, dim=3)
        info = err.error_report()
        # info.to_csv(f"error_data3Dcython.csv",
        #             mode='a', header=False, index=False)
        return info
    except FileNotFoundError:
        return pd.DataFrame()


def parallel_loop_over_folders(save=True, recalc=False, n_cpus=None):
    folders = get_folders(recalc)
    print("Starting folder parsing...")
    with Pool(processes=n_cpus) as pool:
        result = pool.map(get_error_dataframe, folders)

    error_info = pd.concat(result, ignore_ixndex=True)
    print("end")

    if save:
        if recalc:
            error_info.to_csv(
                f"{glb.results_path}/error_data3D.csv",
                mode="w",
                header=True,
                index=False,
            )
        else:
            error_info.to_csv(
                f"{glb.results_path}/error_data3D.csv",
                mode="a",
                header=False,
                index=False,
            )
    return error_info


def serial_loop_over_folders(save=True, recalc=False):
    folders = get_folders(recalc)
    error_info = pd.DataFrame()
    num_of_folders = len(folders)
    print("Serial loop")
    print(f"Number of files to process: {num_of_folders}")
    for idx, folder in enumerate(folders):
        print(f"{idx+1}/{num_of_folders}", folder)
        try:
            err = e.Error(folder, dim=3)
        except FileNotFoundError:
            continue
        info = err.error_report()
        # info.to_csv(f"{glb.results_path}/error_data3Dpythran.csv",
        #            mode='a', header=False, index=False)
        error_info = error_info.append(info, ignore_index=True)
    if save:
        if recalc:
            error_info.to_csv(
                f"{glb.results_path}/error_data3Dpythran_final.csv",
                mode="w",
                index=False,
            )
        else:
            error_info.to_csv(
                f"{glb.results_path}/error_data3Dpythran_final.csv",
                mode="a",
                header=False,
                index=False,
            )
    return error_info
