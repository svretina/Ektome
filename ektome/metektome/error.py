#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import re, os
import matplotlib.cm as cm
import sys
from kuibit.simdir import SimDir

import ektome.metektome.simulation as sim


class Error:
    def __init__(self, vnl_sim_name):
        self.sim_name1 = vnl_sim_name
        self.sim_name2 = "excision" + vnl_sim_name.split("vanilla")[1]

        self.vanilla = sim.Simulation(self.sim_name1)
        self.excision = sim.Simulation(self.sim_name2)
        self.proj_dir = self.vanilla.proj_dir
        self.error_u = abs(self.vanilla.umesh - self.excision.umesh)
        self.error_u_masked = abs(
            self.vanilla.umesh_masked - self.excision.umesh_masked
        )

        self.error_psi = abs(self.vanilla.psimesh - self.excision.psimesh)
        self.error_psi_masked = abs(
            self.vanilla.psimesh_masked - self.excision.psimesh_masked
        )

        self.xmin = self.vanilla.xmin
        self.xmax = self.vanilla.xmax
        self.ymin = self.vanilla.ymin
        self.ymax = self.vanilla.ymax
        self.ex_r = self.excision.ex_r
        (
            self.indu,
            self.x_at_umax,
            self.y_at_umax,
            self.umax,
            self.umaxangle,
        ) = self.find_max("puncture_u")
        (
            self.indpsi,
            self.x_at_psimax,
            self.y_at_psimax,
            self.psimax,
            self.psimaxangle,
        ) = self.find_max("my_psi")
        self.figure_dir = self.vanilla.figure_dir

        self.base_name = self.vanilla.base_name

    def find_max(self, param):
        if param == "puncture_u":
            # print("puncture_u")
            var = self.error_u_masked.copy()
        elif param == "my_psi":
            # print("psi")
            var = self.error_psi_masked.copy()
        else:
            raise TypeError(
                "Invalid choice of parameter\
            in find max"
            )

        ind = np.unravel_index(np.nanargmax(var, axis=None), var.shape)

        angle = np.arctan2(self.excision.ymesh[ind], self.excision.xmesh[ind])
        return (
            ind,
            self.vanilla.xmesh[ind],
            self.vanilla.ymesh[ind],
            var[ind],
            angle,
        )

    def error_report(self):
        error_dict = {
            "x": self.x_at_psimax,
            "y": self.y_at_psimax,
            "q": self.vanilla.mp,
            "b": self.vanilla.par_b,
            "ex_r": self.ex_r,  # 1st BH
            "p1x": self.vanilla.p1x,
            "p1y": self.vanilla.p1y,
            "p1z": self.vanilla.p1z,
            "p1": self.vanilla.p1,
            "s1x": self.vanilla.s1x,
            "s1y": self.vanilla.s1y,
            "s1z": self.vanilla.s1z,
            "s1": self.vanilla.s1,
            # 2nd BH
            "p2x": self.vanilla.p2x,
            "p2y": self.vanilla.p2y,
            "p2z": self.vanilla.p2z,
            "p2": self.vanilla.p2,
            "s2x": self.vanilla.s2x,
            "s2y": self.vanilla.s2y,
            "s2z": self.vanilla.s2z,
            "s2": self.vanilla.s2,
            "angle": self.psimaxangle,
            "max_error_psi": self.psimax,
        }
        return pd.DataFrame.from_dict([error_dict])

    def plot_error_u_xy(self):
        plt.clf()
        plt.tick_params(direction="in")
        im = plt.imshow(
            self.error_u,
            interpolation="spline36",
            cmap="inferno_r",
            origin="lower",
            extent=[self.xmin, self.xmax, self.ymin, self.ymax],
            vmin=abs(self.error_u).min(),
            vmax=abs(self.error_u).max(),
        )
        # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
        plt.scatter(
            [-self.vanilla.par_b, self.vanilla.par_b],
            [0, 0],
            color="black",
            marker="x",
        )

        excision_sphere = plt.Circle(
            (-self.excision.par_b, 0),
            self.excision.ex_r,
            color="black",
            fill=False,
            ls="--",
        )

        rs_radius1 = plt.Circle(
            (-self.excision.par_b, 0),
            self.vanilla.mm / 2,
            color="black",
            fill=False,
        )

        rs_radius2 = plt.Circle(
            (self.excision.par_b, 0),
            self.vanilla.mp / 2,
            color="black",
            fill=False,
        )

        plt.gcf().gca().add_artist(excision_sphere)
        plt.gcf().gca().add_artist(rs_radius1)
        plt.gcf().gca().add_artist(rs_radius2)

        plt.scatter([self.x_at_umax], [self.y_at_umax], color="green", marker="+")

        plt.xlim(self.vanilla.xmin, self.vanilla.xmax)
        plt.ylim(self.vanilla.ymin, self.vanilla.ymax)
        plt.xlabel("x/M")
        plt.ylabel("y/M")
        plt.title(r"Error_u with Excision Radius: %.2f" % self.ex_r)
        plt.colorbar(im, orientation="horizontal")
        plt.savefig("%s-error-u.png" % self.vanilla.sim_dir)
        plt.savefig("%s-error-u.png" % self.excision.sim_dir)

    def plot_error_psi_xy(self):
        plt.clf()
        plt.tick_params(direction="in")
        im = plt.imshow(
            self.error_psi,
            interpolation="spline36",
            cmap="inferno_r",
            origin="lower",
            extent=[self.xmin, self.xmax, self.ymin, self.ymax],
            vmin=abs(self.error_psi).min(),
            vmax=abs(self.error_psi).max(),
        )
        # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
        plt.scatter(
            [-self.vanilla.par_b, self.vanilla.par_b],
            [0, 0],
            color="black",
            marker="x",
        )
        excision_sphere = plt.Circle(
            (-self.excision.par_b, 0),
            self.excision.ex_r,
            color="black",
            fill=False,
            ls="--",
        )

        rs_radius1 = plt.Circle(
            (-self.excision.par_b, 0),
            self.vanilla.mm / 2,
            color="black",
            fill=False,
        )

        rs_radius2 = plt.Circle(
            (self.excision.par_b, 0),
            self.vanilla.mp / 2,
            color="black",
            fill=False,
        )

        plt.gcf().gca().add_artist(excision_sphere)
        plt.gcf().gca().add_artist(rs_radius1)
        plt.gcf().gca().add_artist(rs_radius2)

        plt.scatter([self.x_at_psimax], [self.y_at_psimax], color="red", marker="+")

        plt.xlim(self.vanilla.xmin, self.vanilla.xmax)
        plt.ylim(self.vanilla.ymin, self.vanilla.ymax)
        plt.xlabel("x/M")
        plt.ylabel("y/M")
        plt.title(r"Error_$\psi$ with Excision Radius: %.2f" % self.ex_r)
        plt.colorbar(im, orientation="horizontal")
        plt.savefig("%s-error-psi.png" % self.vanilla.sim_dir)
        plt.savefig("%s-error-psi.png" % self.excision.sim_dir)
