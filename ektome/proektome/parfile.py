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


from string import Template
# from precactus import grid as pg
from jhuki import grid as pg
from jhuki.twopunctures import TwoPunctures
import numpy as np
import ektome.globals as glb


class Parameter_File:
    _lines = """
Cactus::cctk_run_title = "QC-0"
Cactus::cctk_full_warnings         = yes
Cactus::highlight_warning_messages = no
Cactus::terminate       = "time"
Cactus::cctk_final_time = 0.0

ActiveThorns = "IOUtil"

IO::out_dir = $sim_dir

ActiveThorns = "AEILocalInterp"

#ActiveThorns = "BLAS LAPACK"

ActiveThorns = "Fortran"

ActiveThorns = "GSL"

ActiveThorns = "GenericFD"

ActiveThorns = "HDF5"

ActiveThorns = "LocalInterp"

ActiveThorns = "LoopControl"

ActiveThorns = "Slab"

ActiveThorns = "TGRtensor"

ActiveThorns = "SummationByParts"

SummationByParts::order = 4

ActiveThorns = "InitBase"

ActiveThorns = "Carpet CarpetLib CarpetInterp CarpetReduce CarpetSlab"

Carpet::verbose           = no
Carpet::veryverbose       = no
Carpet::schedule_barriers = no
Carpet::storage_verbose   = no
CarpetLib::output_bboxes  = no

# Carpet::max_refinement_levels = 10

Carpet::use_buffer_zones = yes

Carpet::prolongation_order_space = 5
Carpet::prolongation_order_time  = 2

Carpet::convergence_level = 0

Carpet::init_fill_timelevels = yes
#Carpet::init_3_timelevels = yes

Carpet::poison_new_timelevels = yes
CarpetLib::poison_new_memory  = yes

Carpet::output_timers_every      = 5120
CarpetLib::print_timestats_every = 5120
CarpetLib::print_memstats_every  = 5120


ActiveThorns = "NaNChecker"
NaNChecker::check_every     = 1 # 512
NaNChecker::action_if_found = "just warn"
NaNChecker::check_vars      = "
        ADMBase::metric
        ADMBase::curv
        ADMBase::lapse
        ADMBase::shift
        ADMBase::dtlapse
        ADMBase::dtshift
"
ActiveThorns = "Boundary CartGrid3D CoordBase SymBase ReflectionSymmetry CarpetRegrid2"

# driver::ghost_size       = 3
# CoordBase::boundary_size_x_lower     = 3
# CoordBase::boundary_size_y_lower     = 3
# CoordBase::boundary_size_z_lower     = 3
# CoordBase::boundary_size_x_upper     = 3
# CoordBase::boundary_size_y_upper     = 3
# CoordBase::boundary_size_z_upper     = 3

ReflectionSymmetry::reflection_x = "no"
ReflectionSymmetry::reflection_y = "no"
# ReflectionSymmetry::reflection_z = "yes"

ReflectionSymmetry::reflection_z = "no"
ReflectionSymmetry::avoid_origin_x = "no"
ReflectionSymmetry::avoid_origin_y = "no"
ReflectionSymmetry::avoid_origin_z = "no"

#CoordBase::boundary_shiftout_x_lower = 1
CoordBase::boundary_shiftout_z_lower = 1


$grid_section

ActiveThorns = "MoL Time"

MoL::ODE_Method             = "RK4"
MoL::MoL_Intermediate_Steps = 4
MoL::MoL_Num_Scratch_Levels = 1


Time::dtfac = 0.25



ActiveThorns = "ADMBase ADMCoupling ADMMacros CoordGauge SpaceMask StaticConformal TmunuBase"

ADMMacros::spatial_order = 4

ActiveThorns = "TwoPunctures"

ADMBase::lapse_timelevels = 3
ADMBase::shift_timelevels = 3
ADMBase::metric_timelevels = 3

ADMBase::metric_type = "physical"
$initial_data_section

TwoPunctures::keep_u_around  = yes
TwoPunctures::excision_radius= $exr
TwoPunctures::TP_epsilon     = 1.0e-2
TwoPunctures::TP_Tiny        = 1.0e-2
TwoPunctures::verbose        = yes
TwoPunctures::grid_setup_method = evaluation
TwoPunctures::npoints_A      = 60
TwoPunctures::npoints_B      = 60
TwoPunctures::npoints_phi    = 32
TwoPunctures::newton_maxit   = 15
TwoPunctures::newton_tol     = 1e-18

ActiveThorns = "CarpetIOBasic"

IOBasic::outInfo_every      = 128
IOBasic::outInfo_reductions = "norm2"
IOBasic::outInfo_vars       = "
        Carpet::physical_time_per_hour
        # ML_ADMConstraints::H
        # SphericalSurface::sf_radius
        # QuasiLocalMeasures::qlm_spin[0]
"
ActiveThorns = "CarpetIOScalar"

IOScalar::one_file_per_group = yes

# IOScalar::outScalar_every = 128
# IOScalar::outScalar_vars  = "
#     # CarpetReduce::weight
#     # ADMBase::metric
#     # ADMBase::curv
#     # ADMBase::lapse
#     # ADMBase::shift
#     # ADMBase::dtlapse
#     # ADMBase::dtshift
#     # WEYLSCAL4::Psi4r
#     # WEYLSCAL4::Psi4i
#     # ML_ADMConstraints::ML_Ham
#     # ML_ADMConstraints::ML_mom
#     # SphericalSurface::sf_radius
#     # QuasiLocalMeasures::qlm_newman_penrose
#     # QuasiLocalMeasures::qlm_weyl_scalars
#     # QuasiLocalMeasures::qlm_ricci_scalars
#     # QuasiLocalMeasures::qlm_twometric
#     # QuasiLocalMeasures::qlm_killing_vector
#     # QuasiLocalMeasures::qlm_killed_twometric
#     # QuasiLocalMeasures::qlm_invariant_coordinates
#     # QuasiLocalMeasures::qlm_3determinant
# "
ActiveThorns = "CarpetIOASCII"

IOASCII::one_file_per_group = yes

IOASCII::output_symmetry_points = no
IOASCII::out3D_ghosts           = no

# IOASCII::out0D_every = 128
# IOASCII::out0D_vars  = "
#     # Carpet::timing
#     # CarpetReduce::weight
#     # ADMBase::metric
#     # ADMBase::curv
#     # ADMBase::lapse
#     # ADMBase::shift
#     # ADMBase::dtlapse
#     # ADMBase::dtshift
#     # WEYLSCAL4::Psi4r
#     # WEYLSCAL4::Psi4i
#     # ML_ADMConstraints::ML_Ham
#     # ML_ADMConstraints::ML_mom
#     # SphericalSurface::sf_active
#     # SphericalSurface::sf_valid
#     # SphericalSurface::sf_info
#     # SphericalSurface::sf_radius
#     # SphericalSurface::sf_origin
#     # SphericalSurface::sf_coordinate_descriptors
#     # QuasiLocalMeasures::qlm_state
#     # QuasiLocalMeasures::qlm_grid_int
#     # QuasiLocalMeasures::qlm_grid_real
#     # QuasiLocalMeasures::qlm_scalars
#     # QuasiLocalMeasures::qlm_multipole_moments
# "

# IOASCII::out1D_every = 128
# IOASCII::out1D_vars  = "
#     # CarpetReduce::weight
#     # ADMBase::metric
#     # ADMBase::curv
#     # ADMBase::lapse
#     # ADMBase::shift
#     # ADMBase::dtlapse
#     # ADMBase::dtshift
#     # WEYLSCAL4::Psi4r
#     # WEYLSCAL4::Psi4i
#     # ML_ADMConstraints::ML_Ham
#     # ML_ADMConstraints::ML_mom
#     # SphericalSurface::sf_radius
#     # QuasiLocalMeasures::qlm_shapes
#     # QuasiLocalMeasures::qlm_coordinates
#     # QuasiLocalMeasures::qlm_tetrad_l
#     # QuasiLocalMeasures::qlm_tetrad_n
#     # QuasiLocalMeasures::qlm_tetrad_m
#     # QuasiLocalMeasures::qlm_newman_penrose
#     # QuasiLocalMeasures::qlm_weyl_scalars
#     # QuasiLocalMeasures::qlm_ricci_scalars
#     # QuasiLocalMeasures::qlm_twometric
#     # QuasiLocalMeasures::qlm_killing_vector
#     # QuasiLocalMeasures::qlm_killed_twometric
#     # QuasiLocalMeasures::qlm_invariant_coordinates
#     # QuasiLocalMeasures::qlm_3determinant
#     # TwoPunctures::puncture_u
#     # Twopunctures::my_psi
# "

IOASCII::out2D_every = 128
IOASCII::out2D_vars  = "
        # SphericalSurface::sf_radius
	TwoPunctures::puncture_u
	Twopunctures::my_psi
"



Activethorns = "CarpetIOHDF5"

IOHDF5::out_every              = 1
IOHDF5::one_file_per_group     = yes
IOHDF5::output_symmetry_points = no
IOHDF5::out3D_ghosts           = no
IOHDF5::compression_level      = 1
IOHDF5::use_checksums          = yes
# IOHDF5::out_vars               = "
#         # CarpetReduce::weight
#         # ADMBase::metric
#         # ADMBase::curv
#         # ADMBase::lapse
#         # ADMBase::shift
#         # ADMBase::dtlapse
#         # ADMBase::dtshift
#         # WEYLSCAL4::Psi4r
#         # WEYLSCAL4::Psi4i
#         # ML_ADMConstraints::ML_Ham
#         # ML_ADMConstraints::ML_mom
# "

IOHDF5::out2D_every             = 1
IOHDF5::out2D_xy                = yes
IOHDF5::out2D_xz                = no
IOHDF5::out2D_yz                = no

IOHDF5::out2D_vars              = "
	TwoPunctures::puncture_u
	TwoPunctures::my_psi
"
IOHDF5::out3D_every             = 1
IOHDF5::out3D_vars              = "
	TwoPunctures::puncture_u
	TwoPunctures::my_psi
"

# IOHDF5::checkpoint                  = yes
# IO::checkpoint_dir                  = $sim_dir
# IO::checkpoint_ID                   = yes
# IO::checkpoint_every_walltime_hours = 6.0
# IO::checkpoint_on_terminate         = yes

# IO::recover     = "autoprobe"
# IO::recover_dir = $sim_dir
# ActiveThorns = "Formaline"
# ActiveThorns = "TimerReport"
# TimerReport::out_every                  = 5120
# TimerReport::out_filename               = "TimerReport"
# TimerReport::output_all_timers_together = yes
# TimerReport::output_all_timers_readable = yes
# TimerReport::n_top_timers               = 20
"""

    def __init__(self, simulation, excision, N=0):
        """Constructor of Parameter class
        :param simulation: Dictionary with simulation parameters
        :type simulation: dict
        :param base_name: String to name output files
        :type base_name: str
        :param N: Number of refinment radii to add from touching
        limit
        :type N: int
        """
        self.q = simulation.q
        self.N = N
        if excision:
            self.base_name = simulation.excision_base_name
        else:
            self.base_name = simulation.vanilla_base_name

        self.par_b = simulation.b
        # I just pass here the values
        # they should alreayd account for bound orbits
        self.px_minus = simulation.px1
        self.py_minus = simulation.py1
        self.pz_minus = simulation.pz1

        self.sx_minus = simulation.sx1
        self.sy_minus = simulation.sy1
        self.sz_minus = simulation.sz1

        self.px_plus = simulation.px2
        self.py_plus = simulation.py2
        self.pz_plus = simulation.pz2

        self.sx_plus = simulation.sx2
        self.sy_plus = simulation.sy2
        self.sz_plus = simulation.sz2

        self.m_plus = self.q
        self.m_minus = 1.0

        if "vanilla" in self.base_name:
            self.exr = 0.0
        else:
            self.exr = simulation.exr
        self.dx_plus = self.m_plus / 100
        self.dy_plus = self.m_plus / 100
        self.dz_plus = self.m_plus / 100
        self.dx_minus = self.m_minus / 100
        self.dy_minus = self.m_minus / 100
        self.dz_minus = self.m_minus / 100

        self.sim_dir = f"{glb.simulations_path}/{self.base_name}"

        self.EH_minus = glb.min_r
        self.EH_plus = self.m_plus / 2

        self._calculate_new_dx_minus()
        self._calculate_N_before_overlap()
        self._calculate_refinement_radii()
        self._calculate_outer_boundary()
        self._calculate_initial_data_section()
        self._calculate_grid_section()
        self._calculate_parfile_path()
        self._calculate_parfile_content()

    def _calculate_initial_data_section(self):
        spin_plus = (self.sx_plus, self.sy_plus, self.sz_plus)
        spin_minus = (self.sx_minus, self.sy_minus, self.sz_minus)
        momenta_plus = (self.px_plus, self.py_plus, self.pz_plus)
        momenta_minus = (self.px_minus, self.py_minus, self.pz_minus)
        distance = 2*self.par_b
        twopunctures = TwoPunctures(self.m_plus,
                                    self.m_minus,
                                    distance,
                                    momenta_plus=momenta_plus,
                                    momenta_minus=momenta_minus,
                                    chi_plus=spin_plus,
                                    chi_minus=spin_minus)
        self.twopunctures = twopunctures
        self.initial_data_section = twopunctures.parfile_code

    def _calculate_parfile_path(self):
        """Calculates the path of the parameter file."""
        self.parfile_path = f"{glb.parfiles_path}/{self.base_name}.par"

    def write_parfile(self):
        """Writes the parameter file."""
        with open(self.parfile_path, "w") as parfile:
            parfile.write(self.parfile_content)

    def _calculate_parfile_content(self):
        """Calculates the content of the parameter file"""
        self.parfile_content = Template(Parameter_File._lines).substitute(self.__dict__)

    def _calculate_grid_section(self):
        """Calculates the grid section, taking into account
        the AMR, using a PreCactus module.
        """
        center1 = pg.RefinementCenter(
            self.rr_minus,
            dx_fine=self.dx_minus,
            cfl_fine=0.5,
            center_num=1,
            position=(self.twopunctures.coord_x_minus, 0, 0),
        )

        # Same but with different center_num and position
        center2 = pg.RefinementCenter(
            self.rr_plus,
            dx_fine=self.dx_plus,
            cfl_fine=0.5,
            center_num=2,
            position=(self.twopunctures.coord_x_plus, 0, 0),
        )
        grid_not_synced = pg.Grid(
            (center1, center2), outer_boundary=self.outer_boundary
        )
        grid_synced = pg.set_dt_max_grid(grid_not_synced, dt_max=1)
        # print(grid_synced.parfile_code)
        self.grid_section = grid_synced.parfile_code

    def _calculate_outer_boundary(self):
        """Calculates the outer boundary to be 2 times
        the maximum refinement radii's.
        """
        out_bndr = max(max(self.rr_minus, self.rr_plus))
        self.outer_boundary = out_bndr
        while self.outer_boundary < self.par_b + out_bndr:
            self.outer_boundary = self.outer_boundary * 2

    def _calculate_refinement_radii(self):
        """Calculates the refinement radii as the event horizon
        of each black hole times powers of 2.
        """
        self.rr_plus = tuple(
            self.EH_plus * 2 ** level for level in range(self.N_plus + self.N)
        )
        self.rr_minus = tuple(
            self.EH_minus * 2 ** level for level in range(self.N_minus + self.N)
        )

    def _calculate_N_before_overlap(self, d=0):
        """Calculates the maximum number of refinement radii
        before the two radii overlap. This can be seen as a limit
        before they overlap. It's ok if they overlap.
        """
        n_plus = np.log2((2 * self.par_b - d) / (self.m_plus + self.m_minus * 2 ** self.factor)) + 1
        n_minus = n_plus + self.factor
        self.N_plus = int(round(n_plus))
        self.N_minus = int(round(n_minus))

    def _calculate_new_dx_minus(self):
        """Calculates the new resolution of the smaller black
        hole according to the resolution of the big black hole,
        so that its consistent with AMR.
        """
        self.factor = np.floor(np.log2(self.q))
        self.dx_minus = self.dx_plus / 2.0 ** self.factor
        self.dy_minus = self.dy_plus / 2.0 ** self.factor
        self.dz_minus = self.dz_plus / 2.0 ** self.factor
