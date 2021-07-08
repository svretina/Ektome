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

"""This module provides functions to post process the csv data."""

import os
import numpy as np
import pandas as pd
import ektome.globals as glb

data = pd.read_csv(f"{glb.results_path}/error_data3Dpythran_final.csv")
data["q"] = data["q"].round(decimals=0)
qs = data["q"].unique()


def check_if_processed():
    prefix = "error_data_3D_q"
    bools = []
    for q in qs:
        q = int(q)
        if os.path.isfile(f"{prefix}{q}.csv"):
            bools.append(True)
        if os.path.isfile(f"{prefix}{q}_trend.csv"):
            bools.append(True)
    return np.prod(bools)


def seperate_csv_in_q():

    data["b"] = data["b"] / 2
    data["b"] = data["b"].round(decimals=0)
    data["s2x"] = data["s2x"] / (data["q"] ** 2)
    data["s2y"] = data["s2y"] / (data["q"] ** 2)
    data["s2z"] = data["s2z"] / (data["q"] ** 2)
    data["s2"] = data["s2"] / (data["q"] ** 2)
    data["diff"] = data["max_error_psi_theoretical"] - data["max_error_psi"]

    exclude_bad_sims(data, f"{glb.proj_path}/nans.dat")
    for q in qs:
        temp = data[data["q"] == q]
        counts = temp["b"].value_counts().to_frame()
        counts.reset_index(level=0, inplace=True)
        counts.columns = ["b", "counts"]

        bs1 = counts[counts["counts"] > 100]["b"].tolist()
        bs2 = counts[~(counts["counts"] > 100)]["b"].tolist()

        temp1 = temp[temp["b"].isin(bs1)]
        temp2 = temp[temp["b"].isin(bs2)]

        temp1.to_csv(f"{glb.results_path}/error_data_3D_q{int(q)}.csv", index=False)
        temp2.to_csv(
            f"{glb.results_path}/error_data_3D_q{int(q)}_trend.csv", index=False
        )


def exclude_bad_sims(df, filename):
    with open(filename, "r") as fl:
        lines = fl.readlines()

        for line in lines:
            info = sim_info_from_name(line)
            idx = df.index[
                (df["q"] == info["q"])
                & (df["s1x"] == info["s1x"])
                & (df["s1y"] == info["s1y"])
                & (df["s1z"] == info["s1z"])
                & (df["s2x"] == info["s2x"])
                & (df["s2y"] == info["s2y"])
                & (df["s2z"] == info["s2z"])
            ]

            df.drop(index=idx, inplace=True)
    return df


def sim_info_from_name(filename):
    fields = filename.split("_")
    mass_ratio = float(fields[1].split("q")[1])
    b = float(fields[2].split("b")[1])
    s1x = float(fields[3].split("sx1")[1])
    s1y = float(fields[4].split("sy1")[1])
    s1z = float(fields[5].split("sz1")[1])

    s2x = float(fields[6].split("sx2")[1])
    s2y = float(fields[7].split("sy2")[1])
    s2z = float(fields[8].split("sz2")[1])
    return dict(q=mass_ratio, b=b, s1x=s1x, s1y=s1y, s1z=s1z, s2x=s2x, s2y=s2y, s2z=s2z)
