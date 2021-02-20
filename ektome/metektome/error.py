#!/usr/bin/env python3

import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
import ektome.metektome.simulation as sim


class Error:
    def __init__(self, vnl_sim_name):
        self.sim_name1 = vnl_sim_name
        self.sim_name2 = f"excision{vnl_sim_name.split('vanilla')[1]}"

        self.vanilla = sim.Simulation(self.sim_name1)
        self.excision = sim.Simulation(self.sim_name2)

        # self._calculate_error_u()
        self._calculate_error_psi()
        self._calculate_error_psi_theoretical()

        self.ex_r = self.excision.ex_r
        self.psimax = self.calculate_max_with_mask(self.error_psi)
        self.psimaxt = self.calculate_max_with_mask(self.error_psi_t)
        # self.umax = self.calculate_max_with_mask(self.error_u)

    def _circle(self,x,y):
        return (x + self.vanilla.par_b) * (x + self.vanilla.par_b) + y * y

    def _calculate_error_u(self):
        if ((self.vanilla.p1 == 0 ) &
            (self.vanilla.p2 == 0 ) &
            (self.vanilla.s1 == 0 ) &
            (self.vanilla.s2 == 0 )):
            self.error_u = self.excision.u
        else:
            self.error_u = abs(self.vanilla.u - self.excision.u)\
                /self.vanilla.u

    def _calculate_error_psi(self):
        self.error_psi = abs(self.vanilla.psi
                             - self.excision.psi)\
                             /self.vanilla.psi

    def _calculate_error_psi_theoretical(self):
        temp = self.vanilla.mp / (4.0 * self.vanilla.par_b
                                  - self.vanilla.mp )
        self.error_psi_t = temp/self.vanilla.psi

    def _calculate_error_norm_with_mask(self):
        for ref_level, comp_index, unif_grid in self.error_psi:
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
            # plt.hist(data, bins=len(data))
            # plt.savefig("hist.png")
        return norm


    def calculate_max_with_mask(self, var):
        maxs = []
        for ref_level, comp_index, unif_grid in var:
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
            if not np.isnan(data).all():
                maxs.append(np.nanmax(data))
        return np.nanmax(maxs)

    def error_report(self):
        error_dict = {
            "q": self.vanilla.mp,
            "b": self.vanilla.par_b,
            "ex_r": self.ex_r,
            # 1st BH
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
            "max_error_psi": self.psimax,
            "max_error_psi_theoretical": self.psimaxt
        }
        return pd.DataFrame.from_dict([error_dict])
