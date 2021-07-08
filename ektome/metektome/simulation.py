#!/usr/bin/env python3

import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from kuibit.simdir import SimDir
import ektome.globals as glb
import re


class Simulation:
    def __init__(self, sim_name, dim=2):
        self.proj_dir = glb.proj_path
        self.sim_name = sim_name
        self.sim_dir = f"{glb.simulations_path}/{self.sim_name}"
        self.param_file = f"{glb.parfiles_path}/{self.sim_name}.par"
        self.metadata = (
            f"{self.proj_dir}/simulations/{self.sim_name}/TwoPunctures.bbh"
        )

        self.figure_dir = f"{self.proj_dir}/results/figures"
        self.dim = dim
        self.base_name = "_".join(self.sim_dir.split("_")[1:])
        self._read_params()
        self._calculate_norms()
        if self.dim == 2:
            self._load_data_2D()
        else:
            self._load_data_3D()
        # self.umax = self.calculate_max_with_mask(self.u)
        # self.psimax = self.calculate_max_with_mask(self.psi)

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

    def _calculate_norms(self):
        # self.p1 = np.sqrt(self.p1x * self.p1x +
        #                   self.p1y * self.p1y +
        #                   self.p1z * self.p1z)
        # self.p2 = np.sqrt(self.p2x * self.p2x +
        #                   self.p2y * self.p2y +
        #                   self.p2z * self.p2z)
        self.s1 = np.sqrt(
            self.s1x * self.s1x + self.s1y * self.s1y + self.s1z * self.s1z
        )
        self.s2 = np.sqrt(
            self.s2x * self.s2x + self.s2y * self.s2y + self.s2z * self.s2z
        )

    def _read_params(self):
        # bare-mass 1   => m_plus   (+)
        # base-mass 2   => m_minus  (-)
        self.mm = self.read_param("adm-mass2", self.metadata)
        self.mp = self.read_param("adm-mass1", self.metadata)
        self.par_b = self.read_param("separation", self.metadata)
        self.ex_r = self.read_param("excision", self.metadata)

        self.dx = self.read_param("dx", self.param_file)
        self.dy = self.read_param("dy", self.param_file)
        self.dz = self.read_param("dz", self.param_file)

        # Careful !
        # in TwoPunctures.bbh:
        # 1 cooresponds to the (+) puncture
        # 2 coorespoinds to the (-) puncture
        # in our package:
        # 1 cooresponds to the (-) puncture
        # 2 cooresponds to the (+) puncture

        self.s1x = self.read_param("spin2x", self.metadata)
        self.s1y = self.read_param("spin2y", self.metadata)
        self.s1z = self.read_param("spin2z", self.metadata)
        self.s2x = self.read_param("spin1x", self.metadata)
        self.s2y = self.read_param("spin1y", self.metadata)
        self.s2z = self.read_param("spin1z", self.metadata)

    def _load_data_3D(self):
        sim = SimDir(self.sim_dir)
        # u[integration][refinement_level][components]
        self.u = sim.gf.xyz.fields.puncture_u[0]
        self.psi = sim.gf.xyz.fields.my_psi[0]

    def _load_data_2D(self):
        sim = SimDir(self.sim_dir)
        # u[integration][refinement_level][components]
        self.u = sim.gf.xy.fields.puncture_u[0]
        self.psi = sim.gf.xy.fields.my_psi[0]

    def _circle(self, x, y):
        return (x + self.par_b) * (x + self.par_b) + y * y

    def _sphere(self, x, y, z):
        return (x + self.par_b) * (x + self.par_b) + y * y + z * z

    def calculate_max_with_mask_3D(self, var):
        maxs = []
        for _ref_level, _comp_index, unif_grid in var:
            x, y, z = unif_grid.coordinates_from_grid()
            mask = np.ones(unif_grid.data.shape)
            for j in range(y.shape[0]):
                for i in range(x.shape[0]):
                    for k in range(z.shape[0]):
                        if x[i] > 0:
                            mask[i, j, k] = np.nan
                            continue
                        if self._sphere(x[i], y[j], z[k]) < (self.ex_r ** 2):
                            mask[i, j, k] = np.nan

            data = mask * unif_grid.data
            if not np.isnan(data).all():
                maxs.append(np.nanmax(data))
        return np.nanmax(maxs)

    def calculate_max_with_mask(self, var):
        maxs = []
        for _ref_level, _comp_index, unif_grid in var:
            if self.dim == 2:
                x, y = unif_grid.coordinates_from_grid()
            elif self.dim == 3:
                x, y, z = unif_grid.coordinates_from_grid()
            mask = np.ones(unif_grid.data.shape)
            for j in range(y.shape[0]):
                for i in range(x.shape[0]):
                    if x[i] > 0:
                        mask[i, j] = np.nan
                        continue

                    if self.dim == 2 and self._circle(x[i], y[j]) < (
                        self.ex_r ** 2
                    ):
                        mask[i, j] = np.nan
                    if self.dim == 3:
                        for k in range(z.shape[0]):
                            if self._sphere(x[i], y[j], z[k]) < (
                                self.ex_r ** 2
                            ):
                                mask[i, j, k] = np.nan

            data = mask * unif_grid.data
            if not np.isnan(data).all():
                maxs.append(np.nanmax(data))
        return np.nanmax(maxs)
