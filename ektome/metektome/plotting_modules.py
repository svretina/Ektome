#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ektome.globals as glb
import tikzplotlib as tikz

fig_width_pt = 510.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]

params = {'backend': 'ps',
          'font.family': 'Times New Roman',
          'axes.labelsize': 12,
          'legend.fontsize': 10,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'text.usetex': False,
          'figure.figsize': fig_size}

plt.rcParams.update(params)


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

    xs, ys, _ = sim_obj.psi_value_at_excision()
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

    xs, ys, _ = sim_obj.psi_value_at_excision()
    plt.scatter(xs, ys, color="green", marker="+")

    plt.xlim(sim_obj.xmin, sim_obj.xmax)
    plt.ylim(sim_obj.ymin, sim_obj.ymax)
    plt.xlabel("x/M")
    plt.ylabel("y/M")
    plt.title(f"$\psi$ with Excision Radius: {sim_obj.ex_r}")
    plt.colorbar(im, orientation="horizontal")

    savedir = f"{sim_obj.figure_dir}/psi"

    plt.savefig(f"{savedir}/{sim_obj.sim_name}.png")
    plt.savefig(f"{savedir}.png")


def plot_error_curve_separation(error_dict):
    error_dict = f"{glb.results_path}/{error_dict}"
    plt.figure()
    plt.clf()
    if not isinstance(error_dict, pd.core.frame.DataFrame):
        data = pd.read_csv(error_dict, header=0)
    else:
        data = error_dict

    q = int(data['q'].unique()[0])
    data = data[(data['s1']<=0.9) | (data['s2']<=0.9)]
    # chi1 = data["s1z"]
    # chi2 = data["s2z"]
    mtot = 1 + data["q"]
    # chi_eff = ( 1 * chi1 + data["q"] * chi2)/ mtot
    data['b/m'] = data['b']/mtot
    data['b/m'] = data['b/m'].round(decimals=3)
    data.sort_values(by=['b'], inplace=True)

    plt.plot(data['b/m'], data['diff'])
    plt.title(f"mass ratio {q}")
    plt.savefig(f"{glb.results_path}/test_{q}.png")
    tikz.save(f"{glb.figures_path}/separation_{q}.tikz")

def plot_error_xeff(error_dict):
    error_dict = f"{glb.results_path}/{error_dict}"
    plt.figure()
    plt.clf()
    if not isinstance(error_dict, pd.core.frame.DataFrame):
        data = pd.read_csv(error_dict, header=0)
    else:
        data = error_dict

    q = int(data['q'].unique()[0])
    data = data[(data['s1']<=0.9) | (data['s2']<=0.9)]
    bs = data['b'].unique()

    markers = ['D','o']
    colors = ["#f1a340", "#998ec3" ]

    for index, b in enumerate(bs):
        temp_data = data[data['b'] == b]
        chi1 = temp_data["s1z"]
        chi2 = temp_data["s2z"]
        mtot = 1 + temp_data["q"]
        chi_eff = ( 1 * chi1 + temp_data["q"] * chi2)/ mtot

        plt.scatter(chi_eff, temp_data['diff'],
                    marker=markers[index],
                    edgecolor=colors[index],
                    facecolors='none',
                    label=f"b={temp_data['b'].unique()[0]/q}")

    plt.xlabel(r"$\chi_{eff}$")
    plt.title(f"q={q}")
    plt.grid(True, 'major', axis='y')
    plt.legend()
    plt.savefig(f"{glb.results_path}/test_xeff_{q}.png")
    tikz.save(f"{glb.figures_path}/xeff_{q}.tikz")
