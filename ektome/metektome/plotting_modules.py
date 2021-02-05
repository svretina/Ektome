#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os, sys
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe
import matplotlib as mpl

import seaborn as sns
import ektome.globals as glb
import ektome.metektome.osteo as o

# plt.rc('text', usetex=True)
plt.rc("font", family="Times New Roman")

from scipy.constants import golden_ratio


def get_fig_size(width=7, scale=1.0):
    # width = 3.36 # 242 pt
    base_size = np.array([1, 1 / scale / golden_ratio])
    fig_size = width * base_size
    return fig_size


sns.set_context("paper")
sns.set_style(style="whitegrid", rc=None)

mpl.rcParams["figure.dpi"] = 600
mpl.rcParams["figure.figsize"] = get_fig_size()
# mpl.rcParams['text.usetex'] = True
mpl.rc("font", **{"family": "serif", "serif": ["Times New Roman"]})
mpl.rcParams["lines.linewidth"] = 1.5
mpl.rcParams["axes.labelsize"] = 16
mpl.rcParams["xtick.labelsize"] = 13
mpl.rcParams["ytick.labelsize"] = 13
mpl.rcParams["legend.fontsize"] = 13
mpl.rcParams["image.cmap"] = "jet"

# Palette : https://coolors.co
{
    "Claret": "8b1e3f",
    "Dark Purple": "3c153b",
    "Eton Blue": "89bd9e",
    "Gold Crayola": "f0c987",
    "Cinnabar": "db4c40",
}
{
    "C1": "8b1e3f",
    "C2": "3c153b",
    "C3": "89bd9e",
    "C4": "f0c987",
    "C5": "db4c40",
}

# Path effects to define different borders around the detectors.
pe = [
    mpe.Stroke(linewidth=3, foreground="black"),
    mpe.Stroke(foreground="white", alpha=1),
    mpe.Normal(),
]


def plot_error_u_xy(error_obj):
    plt.clf()
    plt.tick_params(direction="in")
    im = plt.imshow(
        error_obj.error_u,
        interpolation="spline36",
        cmap="inferno_r",
        origin="lower",
        extent=[
            error_obj.xmin,
            error_obj.xmax,
            error_obj.ymin,
            error_obj.ymax,
        ],
        vmin=abs(error_obj.error_u).min(),
        vmax=abs(error_obj.error_u).max(),
    )
    # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
    plt.scatter(
        [-error_obj.vanilla.par_b, error_obj.vanilla.par_b],
        [0, 0],
        color="black",
        marker="x",
    )

    excision_sphere = plt.Circle(
        (-error_obj.excision.par_b, 0),
        error_obj.excision.ex_r,
        color="black",
        fill=False,
        ls="--",
    )

    rs_radius1 = plt.Circle(
        (-error_obj.excision.par_b, 0),
        error_obj.vanilla.mm / 2,
        color="black",
        fill=False,
    )

    rs_radius2 = plt.Circle(
        (error_obj.excision.par_b, 0),
        error_obj.vanilla.mp / 2,
        color="black",
        fill=False,
    )

    plt.gcf().gca().add_artist(excision_sphere)
    plt.gcf().gca().add_artist(rs_radius1)
    plt.gcf().gca().add_artist(rs_radius2)

    plt.scatter([error_obj.x_at_umax], [error_obj.y_at_umax], color="green", marker="+")

    plt.xlim(error_obj.vanilla.xmin, error_obj.vanilla.xmax)
    plt.ylim(error_obj.vanilla.ymin, error_obj.vanilla.ymax)
    plt.xlabel("x/M")
    plt.ylabel("y/M")
    plt.title(r"Error_u with Excision Radius: %.2f" % error_obj.ex_r)
    plt.colorbar(im, orientation="horizontal")

    savedir = error_obj.figure_dir + "/error"
    plt.savefig(f"{savedir}/{error_obj.base_name}-error-u.png")
    plt.savefig(f"{error_obj.vanilla.sim_dir}/error-u.png")
    plt.savefig(f"{error_obj.excision.sim_dir}/error-u.png")


def plot_error_psi_xy(error_obj):
    plt.clf()
    plt.tick_params(direction="in")
    im = plt.imshow(
        error_obj.error_psi,
        interpolation="spline36",
        cmap="inferno_r",
        origin="lower",
        extent=[
            error_obj.xmin,
            error_obj.xmax,
            error_obj.ymin,
            error_obj.ymax,
        ],
        vmin=abs(error_obj.error_psi).min(),
        vmax=abs(error_obj.error_psi).max(),
    )
    # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
    plt.scatter(
        [-error_obj.vanilla.par_b, error_obj.vanilla.par_b],
        [0, 0],
        color="black",
        marker="x",
    )
    excision_sphere = plt.Circle(
        (-error_obj.excision.par_b, 0),
        error_obj.excision.ex_r,
        color="black",
        fill=False,
        ls="--",
    )

    rs_radius1 = plt.Circle(
        (-error_obj.excision.par_b, 0),
        error_obj.vanilla.mm / 2,
        color="black",
        fill=False,
    )

    rs_radius2 = plt.Circle(
        (error_obj.excision.par_b, 0),
        error_obj.vanilla.mp / 2,
        color="black",
        fill=False,
    )

    plt.gcf().gca().add_artist(excision_sphere)
    plt.gcf().gca().add_artist(rs_radius1)
    plt.gcf().gca().add_artist(rs_radius2)

    plt.scatter(
        [error_obj.x_at_psimax],
        [error_obj.y_at_psimax],
        color="red",
        marker="+",
    )

    plt.xlim(error_obj.vanilla.xmin, error_obj.vanilla.xmax)
    plt.ylim(error_obj.vanilla.ymin, error_obj.vanilla.ymax)
    plt.xlabel("x/M")
    plt.ylabel("y/M")

    plt.title(rf"Error_$\psi$ with Excision Radius: {error_obj.ex_r}")
    plt.colorbar(im, orientation="horizontal")

    savedir = error_obj.figure_dir + "/error"
    plt.savefig(f"{savedir}/{error_obj.base_name}-error-psi.png")
    plt.savefig(f"{error_obj.vanilla.sim_dir}/error-psi.png")
    plt.savefig(f"{error_obj.excision.sim_dir}/error-psi.png")


def plot_u_xy(sim_obj):
    plt.clf()
    plt.tick_params(direction="in")
    im = plt.imshow(
        sim_obj.umesh,
        interpolation="spline36",
        cmap="inferno_r",
        origin="lower",
        extent=[sim_obj.xmin, sim_obj.xmax, sim_obj.ymin, sim_obj.ymax],
        vmin=abs(sim_obj.umesh).min(),
        vmax=abs(sim_obj.umesh).max(),
    )
    # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
    plt.scatter([-sim_obj.par_b, sim_obj.par_b], [0, 0], color="black", marker="x")
    excision_sphere = plt.Circle(
        (-sim_obj.par_b, 0), sim_obj.ex_r, color="black", fill=False, ls="--"
    )

    rs_radius1 = plt.Circle(
        (-sim_obj.par_b, 0), sim_obj.mm / 2, color="black", fill=False
    )

    rs_radius2 = plt.Circle(
        (sim_obj.par_b, 0), sim_obj.mp / 2, color="black", fill=False
    )

    xs, ys, psis = sim_obj.psi_value_at_excision()
    plt.scatter(xs, ys, color="green", marker="+")
    plt.gcf().gca().add_artist(excision_sphere)
    plt.gcf().gca().add_artist(rs_radius1)
    plt.gcf().gca().add_artist(rs_radius2)

    plt.xlim(sim_obj.xmin, sim_obj.xmax)
    plt.ylim(sim_obj.ymin, sim_obj.ymax)
    plt.xlabel("x/M")
    plt.ylabel("y/M")
    plt.title(r"u with Excision Radius: %.2f" % sim_obj.ex_r)
    plt.colorbar(im, orientation="horizontal")

    savedir = sim_obj.figure_dir + "/u"
    plt.savefig(f"{savedir}/{sim_obj.sim_name}.png")
    plt.savefig(f"{sim_obj.sim_dir}/u.png")


def plot_psi_xy(sim_obj):
    plt.clf()
    plt.tick_params(direction="in")
    im = plt.imshow(
        sim_obj.psimesh,
        interpolation="spline36",
        cmap="inferno_r",
        origin="lower",
        extent=[sim_obj.xmin, sim_obj.xmax, sim_obj.ymin, sim_obj.ymax],
        vmin=abs(sim_obj.psimesh).min(),
        vmax=abs(sim_obj.psimesh).max(),
    )
    # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
    plt.scatter([-sim_obj.par_b, sim_obj.par_b], [0, 0], color="black", marker="x")
    excision_sphere = plt.Circle(
        (-sim_obj.par_b, 0), sim_obj.ex_r, color="black", fill=False, ls="--"
    )

    rs_radius1 = plt.Circle(
        (-sim_obj.par_b, 0), sim_obj.mm / 2, color="black", fill=False
    )

    rs_radius2 = plt.Circle(
        (sim_obj.par_b, 0), sim_obj.mp / 2, color="black", fill=False
    )

    plt.gcf().gca().add_artist(excision_sphere)
    plt.gcf().gca().add_artist(rs_radius1)
    plt.gcf().gca().add_artist(rs_radius2)

    xs, ys, psis = sim_obj.psi_value_at_excision()
    plt.scatter(xs, ys, color="green", marker="+")

    plt.xlim(sim_obj.xmin, sim_obj.xmax)
    plt.ylim(sim_obj.ymin, sim_obj.ymax)
    plt.xlabel("x/M")
    plt.ylabel("y/M")
    plt.title(f"$\psi$ with Excision Radius: {sim_obj.ex_r}")
    plt.colorbar(im, orientation="horizontal")

    savedir = sim_obj.figure_dir + "/psi"

    plt.savefig(f"{savedir}/{sim_obj.sim_name}.png")
    plt.savefig(f"{sim_obj.sim_dir}/psi.png")


def plot_error_curve_separation(error_dict, mode):
    plt.figure()
    if not isinstance(error_dict, pd.core.frame.DataFrame):
        data = pd.read_csv(error_dict, header=0)
    else:
        data = error_dict

    # q = data["q"][0]
    # b = data["b"].values
    # theoretical = q/ ( 4.0 * b - q)
    if mode == "momentum":
        momentum = data[(data["s1"] == 0) & (data["s2"] == 0)]
        groups = ["p1", "p2", "p1x", "p1y", "p1z", "p2x", "p2y", "p2z"]
        for state, frame in momentum.groupby(by=groups):
            if state == (0.0, 0.0, -0.0, 0.0, 0.0, 0.0, -0.0, -0.0):
                continue
            frame = frame.sort_values(by=["b"])
            b = frame["b"].values
            q = frame["q"].values[0]
            theoretical = q / (4.0 * b - q)
            error_psi = frame["max_error_psi"].values
            dif = theoretical - error_psi
            # print(dif, frame["p1"].values[0], frame["p2"].values[0])
            # if np.any( dif < 0):
            #     print(frame[["p1","p2"]])
            plt.semilogy(
                b,
                error_psi,
                color="#888888",
                marker=".",
                alpha=0.4,
                label=rf"$p_1$={state[0]},\
                         $p_2$={state[1]}",
            )
    elif mode == "spin":
        spin = data[(data["p1"] == 0) & (data["p2"] == 0)]
        groups = ["s1", "s2", "s1x", "s1y", "s1z", "s2x", "s2y", "s2z"]
        for state, frame in spin.groupby(by=groups):
            if state == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0):
                continue
            frame = frame.sort_values(by=["b"])
            b = frame["b"].values
            q = frame["q"].values[0]
            theoretical = q / (4.0 * b - q)
            error_psi = frame["max_error_psi"].values
            dif = theoretical - error_psi
            # if np.any( dif < 0):
            #     print(frame[["p1","p2"]])
            plt.semilogy(
                b,
                error_psi,
                color="#888888",
                marker=".",
                alpha=0.4,
                label=rf"$p_1$={state[0]},\
                         $p_2$={state[1]}",
            )
    else:
        sys.exit()

    q = frame["q"].values[0]
    theoretical = q / (4.0 * b - q)
    plt.semilogy(
        b,
        theoretical,
        color="#cc0000",
        marker="+",
        lw=1.5,
        label=r"$p_1$=0.0, $p_2$=0.0",
    )
    # #"#cc0000",\
    # plt.plot(b,theoretical,"rx",label="p=0.0")
    plt.title(rf"Maximum Error in $\Psi$ for BHs with {mode}")
    plt.xlabel("par_b variable")
    plt.ylabel("Maximum Error")
    # plt.legend()
    results_dir = os.getcwd() + "/results/figures/error-analysis/"
    plt.tight_layout()
    plt.savefig(f"{results_dir}error_separation-{mode}.png")


def plot_error_curve_momentum(error_dict):
    plt.clf()
    if not isinstance(error_dict, pd.core.frame.DataFrame):
        data = pd.read_csv(error_dict, header=0)
    else:
        data = error_dict
    print(data)
    q = data["q"][0]
    b = data["b"].values
    theoretical = q / (4.0 * b - q)

    # plt.plot(b,theoretical,marker="x",label="p=0.0")
    markers = ["o", "+", "*", "1", "D", "X", ".", "^", "s", "H"]
    for state, frame in data.groupby(by=["b"]):
        counter = np.random.randint(0, len(markers))
        plt.scatter(
            frame["p1"],
            frame["max_error_psi"],
            marker=markers[counter],
            label=f"b={state}",
        )

    # plt.plot(b,theoretical,marker="x",label="p=0.0")

    plt.title(r"Maximum Error in $\Psi$")
    plt.xlabel("momentum")
    plt.ylabel("Maximum Error")
    plt.legend()
    results_dir = os.getcwd() + "/results/figures/error-analysis/"
    plt.savefig(f"{results_dir}error_momentum.png")
