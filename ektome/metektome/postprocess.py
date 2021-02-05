#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os

import ektome.metektome.error as e
import ektome.metektome.plotting_modules as pm
import ektome.metektome.simulation as sim


def loop_over_folders(save=True, plot=True):
    sim_dir = os.getcwd() + "/simulations"
    folders = os.listdir(sim_dir)
    error_info = pd.DataFrame()

    for folder in folders:
        simulation = sim.Simulation(folder)
        if plot:
            pm.plot_psi_xy(simulation)
            pm.plot_u_xy(simulation)

        if "vanilla" in folder:
            err = e.Error(folder)

            info = err.error_report()
            error_info = error_info.append(info, ignore_index=True)

            if plot:
                pm.plot_error_psi_xy(err)
                pm.plot_error_u_xy(err)

    if save:
        results_dir = os.getcwd() + "/results"
        error_info.to_csv(results_dir + "/error_data.csv", index=False)
    return error_info


def error_analysis():
    results_dir = os.getcwd() + "/results"
    error_data = pd.read_csv(results_dir + "/error_data.csv", header=0)

    # print(error_data)
    pm.plot_error_curve_momentum(error_data)
    pm.plot_error_curve_separation(error_data)
