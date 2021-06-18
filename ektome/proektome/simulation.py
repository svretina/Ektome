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

"""This module provides :py:class:`~Simulation`.
"""
import numpy as np


class Simulation:
    def __init__(self,q=None,b=None,
                 px1=None,py1=None,pz1=None,
                 sx1=None,sy1=None,sz1=None,
                 px2=None,py2=None,pz2=None,
                 sx2=None,sy2=None,sz2=None,
                 exr=None):
        self.q = q
        self.b = b
        self.px1 = px1
        self.py1 = py1
        self.pz1 = pz1

        self.sx1 = sx1
        self.sy1 = sy1
        self.sz1 = sz1

        self.px2 = px2
        self.py2 = py2
        self.pz2 = pz2

        self.sx2 = sx2
        self.sy2 = sy2
        self.sz2 = sz2
        self.exr = exr
        self.vanilla_base_name , self.excision_base_name = self.create_base_names()

    @staticmethod
    def norm(x, y, z):
        return np.sqrt(x * x + y * y + z * z)

    def calculate_momentum(self):
        """Calculates the norm of the momentum vector for the simulation.
        :returns: The norm value
        :rtype: np.float
        """
        return self.norm(self.px1, self.py1, self.pz1), self.norm(self.px2, self.py2, self.pz2)


    def calculate_spin(self):
        """Calculates the norm of the effective spin for the simulation.
        :returns: The norm value
        :rtype: np.float
        """

        return self.norm(self.sx1, self.sy1, self.sz1), self.norm(self.sx2, self.sy2, self.sz2)

    def create_base_names(self):
        """Creates a name which serves as a base to name various
        output files and directories.

        :returns: Base names for vanilla and excision case.
        :rtype: str

        """
        tmp = []
        d = self.__dict__
        keys = ["q","b", "sx1","sy1","sz1","sx2","sy2","sz2"]
        for i in keys:
            if "q" in i:
                tmp.append(f"_{i}{int(d[i])}")
            elif "b" in i:
                tmp.append(f"_{i}{int(d[i])}")
            else:
                tmp.append(f"_{i}{d[i]}")

        suffix = "".join(tmp)
        vnl_name = f"vanilla{suffix}"
        exc_name = f"excision{suffix}"
        return vnl_name, exc_name
