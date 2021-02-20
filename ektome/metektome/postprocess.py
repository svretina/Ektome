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

import ektome.globals as glb
import ektome.metektome.error as e
import ektome.metektome.simulation as sim
import ektome.metektome.plotting_modules as pm


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
    sim_folders = os.listdir(glb.simulations_path)
    vanilla_folders = [x for x in sim_folders if x.startswith("vanilla")]
    if recalc:
        exclude_folders = get_unfinished_sims(vanilla_folders)
    else:
        exclude_folders = get_processed_sims(vanilla_folders)
        + get_unfinished_sims()

    folders = list((Counter(vanilla_folders)
                   - Counter(exclude_folders)).elements())
    return folders

def get_error_dataframe(folder):
    print(folder)
    err = e.Error(folder)
    info = err.error_report()
    return info

def parallel_loop_over_folders(save=True, recalc=False):
    folders = get_folders(recalc)
    with Pool() as pool:
        result = pool.map(get_error_dataframe, folders)
    error_info = pd.concat(result, ignore_index=True)
    if save:
        if recalc:
            error_info.to_csv(f"{glb.results_path}/error_data.csv",
                          mode='w', header=True, index=False)
        else:
            error_info.to_csv(f"{glb.results_path}/error_data.csv",
                          mode='a', header=False, index=False)
    return error_info

def serial_loop_over_folders(save=True, plot=False, recalc=False):
    folders = get_folders(recalc)
    error_info = pd.DataFrame()
    for folder in folders:
        print(folder)
        err = e.Error(folder)
        info = err.error_report()
        error_info = error_info.append(info, ignore_index=True)
    if save:
        if recalc:
            error_info.to_csv(f"{glb.results_path}/error_data.csv",
                          mode='w', index=False)
        else:
            error_info.to_csv(f"{glb.results_path}/error_data.csv",
                          mode='a', header=False, index=False)
    return error_info


def error_analysis():
    error_data = pd.read_csv(f"{glb.results_path}/error_data.csv", header=0)
    # ????--------------------- v ----- check this
    # pm.plot_error_curve_momentum(error_data)
    pm.plot_error_curve_separation(error_data, "momentum")
