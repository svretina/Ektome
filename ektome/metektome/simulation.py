#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from kuibit.simdir import SimDir
import ektome.globals as glb
import re


class Simulation:

    def __init__(self, sim_name):
        self.proj_dir = glb.proj_path
        self.sim_name = sim_name
        self.sim_dir = f"{glb.simulations_path}/{self.sim_name}"
        self.param_file = f"{glb.parfiles_path}/{self.sim_name}.par"
        self.metadata = f"{self.proj_dir}/simulations/{self.sim_name}/TwoPunctures.bbh"

        self.figure_dir = f"{self.proj_dir}/results/figures"

        self.base_name = "_".join(self.sim_dir.split("_")[1:])
        (
            self.mp,
            self.mm,
            self.par_b,
            self.ex_r,
            self.dx,
            self.dy,
            self.dz,
            self.p1x,
            self.p1y,
            self.p1z,
            self.p2x,
            self.p2y,
            self.p2z,
            self.s1x,
            self.s1y,
            self.s1z,
            self.s2x,
            self.s2y,
            self.s2z,
        ) = self.get_params()

        self.p1 = np.sqrt(
            self.p1x * self.p1x + self.p1y * self.p1y + self.p1z * self.p1z
        )

        self.p2 = np.sqrt(
            self.p2x * self.p2x + self.p2y * self.p2y + self.p2z * self.p2z
        )

        self.s1 = np.sqrt(
            self.s1x * self.s1x + self.s1y * self.s1y + self.s1z * self.s1z
        )

        self.s2 = np.sqrt(
            self.s2x * self.s2x + self.s2y * self.s2y + self.s2z * self.s2z
        )

        (
            self.xmesh,
            self.ymesh,
            self.umesh,
            self.psimesh,
            self.bounds,
            self.lenx,
            self.leny,
        ) = self.load_data()

        self.yvec = self.ymesh[:, 0]
        self.xvec = self.xmesh[0, :]

        self.xmin = self.bounds[0, 0]
        self.xmax = self.bounds[0, 1]
        self.ymin = self.bounds[1, 0]
        self.ymax = self.bounds[1, 1]
        self.umesh_masked = self.mask_data("puncture_u")
        self.psimesh_masked = self.mask_data("my_psi")
        self.umax_ind, self.umax = self.find_max("puncture_u")
        self.psimax_ind, self.psimax = self.find_max("my_psi")

    @staticmethod
    def get_line(string, fl):
        file = open(fl, "r")
        for line in file:
            if re.search(string, line):
                tmp = line
                break
        file.close()
        return tmp

    @staticmethod
    def get_nvalue(line):
        val = line.split("=")[1]
        return float(val)

    def read_param(self, param, fl):
        return self.get_nvalue(self.get_line(param, fl))

    def get_params(self):
        mm = self.read_param("bare-mass2", self.metadata)
        mp = self.read_param("bare-mass1", self.metadata)
        par_b = self.read_param("position1x", self.metadata)
        ex_r = self.read_param("excision", self.metadata)

        dx = self.read_param("dx", self.param_file)
        dy = self.read_param("dy", self.param_file)
        dz = self.read_param("dz", self.param_file)

        # Momenta
        p1x = self.read_param("momentum1x", self.metadata)
        p1y = self.read_param("momentum1y", self.metadata)
        p1z = self.read_param("momentum1z", self.metadata)
        p2x = self.read_param("momentum2x", self.metadata)
        p2y = self.read_param("momentum2y", self.metadata)
        p2z = self.read_param("momentum2z", self.metadata)
        # Spins
        s1x = self.read_param("spin1x", self.metadata)
        s1y = self.read_param("spin1y", self.metadata)
        s1z = self.read_param("spin1z", self.metadata)
        s2x = self.read_param("spin2x", self.metadata)
        s2y = self.read_param("spin2y", self.metadata)
        s2z = self.read_param("spin2z", self.metadata)
        return (
            mp,
            mm,
            par_b,
            ex_r,
            dx,
            dy,
            dz,
            p1x,
            p1y,
            p1z,
            p2x,
            p2y,
            p2z,
            s1x,
            s1y,
            s1z,
            s2x,
            s2y,
            s2z,
        )

    def load_data(self):
        sim = SimDir(self.sim_dir)

        u = sim.gf.xy.fields.puncture_u[0]
        psi = sim.gf.xy.fields.my_psi[0]

        # Use this with AMR
        # u   = u.merge_refinement_levels(resample=True)
        # psi = psi.merge_refinement_levels(resample=True)
        xmesh = u.coordinates()[0][0][0][:].T
        ymesh = u.coordinates()[1][0][0][:].T

        # filename = "/twopunctures-"+str(var)+".xy.asc-pp.txt"
        # data = pd.read_csv(self.sim_dir+filename,sep=",",comment='#',usecols=[9,10,12],names=["x","y","var"])
        # x = data['x'].values
        # y = data['y'].values

        xmin = xmesh.min()
        xmax = xmesh.max()

        ymin = ymesh.min()
        ymax = ymesh.max()

        lenx = xmesh.shape[1]
        leny = xmesh.shape[0]

        u = u[0][0].data.T
        psi = psi[0][0].data.T
        # var = data['var'].values

        # xmesh = x.reshape((leny,lenx))
        # ymesh = y.reshape((leny,lenx))
        # var = var.reshape((leny,lenx))
        bounds = np.array([[xmin, xmax], [ymin, ymax]])
        return xmesh, ymesh, u, psi, bounds, lenx, leny

    def mask_data(self, param):
        radius = self.ex_r
        x1d = self.xmesh.ravel()
        y1d = self.ymesh.ravel()
        if param == "puncture_u":
            var = self.umesh.copy()
        elif param == "my_psi":
            var = self.psimesh.copy()
        var1d = var.ravel()
        N = x1d.shape[0]
        for i in range(N):
            circle = (x1d[i] + self.par_b) * (x1d[i] + self.par_b) + y1d[i] * y1d[i]

            if circle < (radius * radius):
                var1d[i] = np.nan

            if x1d[i] > 0:
                var1d[i] = np.nan
        maskedmesh = var1d.reshape((self.leny, self.lenx))
        return maskedmesh

    def find_max(self, param):
        if param == "puncture_u":
            # print("puncture_u")
            var = self.umesh_masked
        elif param == "my_psi":
            # print("psi")
            var = self.psimesh_masked
        else:
            raise TypeError(
                "Invalid choice of parameter\
                in find max"
            )

        ind = np.unravel_index(np.nanargmax(var, axis=None), var.shape)

        return ind, var[ind]

    def psi_value_at_excision(self):
        x1d = self.xmesh.copy().ravel()
        y1d = self.ymesh.copy().ravel()

        var = self.psimesh.copy()
        var1d = var.ravel()
        N = x1d.shape[0]
        psis = []
        xs = []
        ys = []
        for i in range(N):
            circle = (x1d[i] + self.par_b) * (x1d[i] + self.par_b) + y1d[i] * y1d[i]
            if circle == (self.ex_r * self.ex_r):
                psis.append(var1d[i])
                xs.append(x1d[i])
                ys.append(y1d[i])
        return xs, ys, psis

    def plot_u_xy(self):
        plt.clf()
        plt.tick_params(direction="in")
        im = plt.imshow(
            self.umesh,
            interpolation="spline36",
            cmap="inferno_r",
            origin="lower",
            extent=[self.xmin, self.xmax, self.ymin, self.ymax],
            vmin=abs(self.umesh).min(),
            vmax=abs(self.umesh).max(),
        )
        # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
        plt.scatter([-self.par_b, self.par_b], [0, 0], color="black", marker="x")
        excision_sphere = plt.Circle(
            (-self.par_b, 0), self.ex_r, color="black", fill=False, ls="--"
        )

        rs_radius1 = plt.Circle(
            (-self.par_b, 0), self.mm / 2, color="black", fill=False
        )

        rs_radius2 = plt.Circle((self.par_b, 0), self.mp / 2, color="black", fill=False)

        xs, ys, _ = self.psi_value_at_excision()
        plt.scatter(xs, ys, color="green", marker="+")
        plt.gcf().gca().add_artist(excision_sphere)
        plt.gcf().gca().add_artist(rs_radius1)
        plt.gcf().gca().add_artist(rs_radius2)

        plt.xlim(self.xmin, self.xmax)
        plt.ylim(self.ymin, self.ymax)
        plt.xlabel("x/M")
        plt.ylabel("y/M")
        plt.title(r"u with Excision Radius: %.2f" % self.ex_r)
        plt.colorbar(im, orientation="horizontal")
        plt.savefig("%s-u.png" % self.sim_dir)

    def plot_psi_xy(self):
        plt.clf()
        plt.tick_params(direction="in")
        im = plt.imshow(
            self.psimesh,
            interpolation="spline36",
            cmap="inferno_r",
            origin="lower",
            extent=[self.xmin, self.xmax, self.ymin, self.ymax],
            vmin=abs(self.psimesh).min(),
            vmax=abs(self.psimesh).max(),
        )
        # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
        plt.scatter([-self.par_b, self.par_b], [0, 0], color="black", marker="x")
        excision_sphere = plt.Circle(
            (-self.par_b, 0), self.ex_r, color="black", fill=False, ls="--"
        )

        rs_radius1 = plt.Circle(
            (-self.par_b, 0), self.mm / 2, color="black", fill=False
        )

        rs_radius2 = plt.Circle((self.par_b, 0), self.mp / 2, color="black", fill=False)

        plt.gcf().gca().add_artist(excision_sphere)
        plt.gcf().gca().add_artist(rs_radius1)
        plt.gcf().gca().add_artist(rs_radius2)

        xs, ys, _ = self.psi_value_at_excision()
        plt.scatter(xs, ys, color="green", marker="+")

        plt.xlim(self.xmin, self.xmax)
        plt.ylim(self.ymin, self.ymax)
        plt.xlabel("x/M")
        plt.ylabel("y/M")
        plt.title(r" $\psi$ with Excision Radius: %.2f" % self.ex_r)
        plt.colorbar(im, orientation="horizontal")
        plt.savefig("%s-psi.png" % self.sim_dir)
