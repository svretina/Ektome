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

"""This module provides helper functions to calculate keplerian
orbits.
The functions available are:
"""

import os
import sys
import math
import numpy as np
from PyAstronomy import pyasl
from jhuki.twopunctures import prepare_quasicircular_inspiral

import ektome.globals as glb


class Binary:
    def __init__(self, q):
        """Constructor for Binary class.
        :param m1: Mass of 1st binary component.
        :type m1: float
        :param m2: Mass of 2nd binary component.
        :type m2: float
        """
        self.m1 = 1
        self.m2 = q
        self.m_tot = self.m1 + self.m2
        self.mu = (self.m1 * self.m2) / self.m_tot
        self.m_chirp = (self.m1 * self.m2) ** (3.0 / 5.0) / (self.m_tot) ** (
            1.0 / 5.0
        )
        self.f_isco = 1.0 / (6 * np.sqrt(6) * 2.0 * np.pi) * (1.0 / self.m_tot)

    def semimajor(self, N_orb):
        """Calculates the semi major axis for a binary of
        masses $m_1$ and $m_2$ for $N_{orb}$ number of
        orbits outwards from the ISCO.

        :param N_orb: Number of orbits from ISCO and outwards
        :type N_orb: float
        :returns: The semi major axis for $N_{orb}$ number of
        orbits outwards of ISCO
        :rtype: float
        """
        denom = (2.0 * np.pi) ** (2.0 / 3.0)
        term = 2.0 * N_orb * 32 * np.pi ** (8.0 / 3.0) * self.m_chirp ** (
            5.0 / 3.0
        ) + self.f_isco ** (-5.0 / 3.0)
        return (self.m_tot ** (1.0 / 3.0) / denom) * term ** (2.0 / 5.0)

    @staticmethod
    def v_at(f, r):
        """Calculates the velocity at given frequency
        and distance.

        :param f: Frequency.
        :type f: float
        :param r: Distance.
        :type r: float
        :returns: The velocity
        :rtype: float
        """
        return 2 * np.pi * f * r

    def freq_at(self, semi_major):
        """Calculates the frequency at a given semi major axis.
        :param semi_major: The semi major axis to compute
        the frequency at
        :type semi_major: float
        :returns: The frequency.
        :rtype: float
        """
        return 1 / (2 * np.pi) * np.sqrt(self.m_tot / semi_major ** 3)

    def circular_binary(self, N_orb):
        """Calculates the velocity vectors for a binary.

        :param N_orb: The number of orbits from ISCO
        :type N_orb: int
        :returns:
        """
        semi_major = round(self.semimajor(N_orb))
        f = self.freq_at(semi_major)
        ke1 = pyasl.KeplerEllipse(-semi_major, 1 / f)
        ke2 = pyasl.KeplerEllipse(semi_major, 1 / f)

        vel1 = np.round(ke1.xyzVel(0.0), 2)
        vel2 = np.round(ke2.xyzVel(0.0), 2)
        return semi_major, vel1, vel2

    @staticmethod
    def quasicircular_inspiral(q, distance, spin_minus, spin_plus):
        twopunctures = prepare_quasicircular_inspiral(
            q,
            distance,
            1 + q,
            spin_plus,
            spin_minus,
        )
        p1 = np.asarray(twopunctures.momenta_minus, dtype=np.float64)
        p2 = np.asarray(twopunctures.momenta_plus, dtype=np.float64)

        return p1, p2
