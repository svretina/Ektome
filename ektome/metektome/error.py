#!/usr/bin/env python3

import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
import bottleneck as bn

import ektome.metektome.simulation as sim
import ektome.metektome.cython_funcs as cf
import ektome.metektome.search3D as s

class Error:
    def __init__(self, vnl_sim_name,dim=2):
        self.sim_name1 = vnl_sim_name
        self.sim_name2 = f"excision{vnl_sim_name.split('vanilla')[1]}"
        self.dim = dim
        self.error_u = None
        self.error_psi = None
        try:
            self.vanilla = sim.Simulation(self.sim_name1,self.dim)
            self.excision = sim.Simulation(self.sim_name2,self.dim)
        except:
            print("File not found")
            raise FileNotFoundError("File not found")

        self.q = self.vanilla.mp/self.vanilla.mm
        # self._calculate_error_u()
        self._calculate_error_psi()
        self._calculate_error_psi_theoretical()

        self.ex_r = self.excision.ex_r
        self.par_b = self.vanilla.par_b
        self.psimax = self.calculate_max_with_mask3D(self.error_psi)
        self.psimaxt = self.calculate_max_with_mask3D(self.error_psi_t)
        # self.umax = self.calculate_max_with_mask(self.error_u)

    def _circle(self,x,y):
        return (x + self.vanilla.par_b) * (x + self.vanilla.par_b)\
            + y * y

    def _sphere(self,x,y,z):
        return (x + self.vanilla.par_b) * (x + self.vanilla.par_b) + y * y + z * z

    def _calculate_error_u(self):
        if ((self.vanilla.s1 == 0) &
            (self.vanilla.s2 == 0)):
            self.error_u = self.excision.u
        else:
            self.error_u = abs(self.vanilla.u - self.excision.u)/self.vanilla.u

    def _calculate_error_psi(self):
        self.error_psi = abs(self.vanilla.psi - self.excision.psi)/self.vanilla.psi

    def _calculate_error_psi_theoretical(self):
        temp = self.q / (4.0 * self.vanilla.par_b - self.q)
        self.error_psi_t = temp/self.vanilla.psi

    def _calculate_error_norm_with_mask(self):
        for _ref_level, _comp_index, unif_grid in self.error_psi:
            x, y = unif_grid.coordinates_from_grid()
            mask = np.ones(unif_grid.data.shape)
            for j in range(y.shape[0]):
                for i in range(x.shape[0]):
                    if x[i] > 0:
                        mask[i,j] = np.nan
                        continue
                    if self._circle(x[i],y[j]) < (self.ex_r**2):
                        mask[i,j] = np.nan
            data = mask * unif_grid.data
            data = data[~np.isnan(data)]
            data = data.reshape(-1)
            norm = np.linalg.norm(data)/np.sqrt(len(data))
        return norm


    def _sphere(self,x,y,z):
        return (x + self.par_b) * (x + self.par_b) + y * y + z * z

    def calculate_max_with_mask3D(self,var):
        maxs = []
        for _ref_level, _comp_index, unif_grid in var:
            x, y, z = unif_grid.coordinates_from_grid()
            dx = unif_grid.dx


            ybounds = [-self.ex_r - dx[1],
                       self.ex_r + dx[1]]
            zbounds = [-self.ex_r - dx[2],
                       self.ex_r + dx[2]]
            xbounds = [-self.par_b - self.ex_r - dx[0],
                       -self.par_b + self.ex_r + dx[0]]
            # mask = np.ones(unif_grid.data.shape)*np.nan
            # print(unif_grid.data.shape)
            mask = s.loop(self.par_b, self.ex_r,
                          x, y, z,
                          xbounds, ybounds, zbounds,
                          unif_grid.data.shape)
            # mask = np.asarray(mask)
            data = mask * unif_grid.data
            if not bn.allnan(data):
                maxs.append(bn.nanmax(data))
        return bn.nanmax(maxs)

    def calculate_max_with_mask(self, var):
        maxs = []
        for _ref_level, _comp_index, unif_grid in var:
            if self.vanilla.dim == 2:
                x, y = unif_grid.coordinates_from_grid()
                dx, dy = unif_grid.dx
                ybounds = [-self.ex_r - dy, self.ex_r + dy]
                xbounds = [-self.par_b - self.ex_r - dx,
                           -self.par_b + self.ex_r + dx]
                mask = np.ones(unif_grid.data.shape)*np.nan
                for j in range(y.shape[0]):
                    if y[j] < ybounds[0] or y[j] > ybounds[1]:
                        continue
                    for i in range(x.shape[0]):
                        if x[i] < xbounds[0] or x[i] > xbounds[1]:
                            continue
                        if self.dim == 2 and self._circle(x[i],y[j]) >= (self.ex_r**2):
                            mask[i,j] = 1

            elif self.vanilla.dim == 3:
                x, y, z = unif_grid.coordinates_from_grid()
                dx = unif_grid.dx

                ybounds = [-self.ex_r - dx[1], self.ex_r + dx[1]]
                zbounds = [-self.ex_r - dx[2], self.ex_r + dx[2]]
                xbounds = [-self.par_b - self.ex_r - dx[0],
                           -self.par_b + self.ex_r + dx[0]]
                mask = np.ones(unif_grid.data.shape)*np.nan
                for j in range(y.shape[0]):
                    if y[j] < ybounds[0] or y[j] > ybounds[1]:
                        continue
                    for i in range(x.shape[0]):
                        if x[i] < xbounds[0] or x[i] > xbounds[1]:
                            continue
                        for k in range(z.shape[0]):
                            if z[k] < zbounds[0] or z[k] > zbounds[1]:
                                continue
                            if self._sphere(x[i],y[j],z[k]) >= (self.ex_r**2):
                                mask[i,j,k] = 1


            data = mask * unif_grid.data
            if not bn.allnan(data):
                maxs.append(bn.nanmax(data))
        return bn.nanmax(maxs)



    def error_report(self):
        error_dict = {
            "q": self.q,
            "b": self.vanilla.par_b,
            "ex_r": self.ex_r,
            # 1st BH
            "s1x": self.vanilla.s1x,
            "s1y": self.vanilla.s1y,
            "s1z": self.vanilla.s1z,
            "s1": self.vanilla.s1,
            # 2nd BH
            "s2x": self.vanilla.s2x,
            "s2y": self.vanilla.s2y,
            "s2z": self.vanilla.s2z,
            "s2": self.vanilla.s2,
            "max_error_psi": self.psimax,
            "max_error_psi_theoretical": self.psimaxt
        }
        return pd.DataFrame.from_dict([error_dict])
