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


"""This module provides helper functions to submit a simulation.
The functions available are:
- :py:func:`~.create_base_names`: Creates a name for various output files and directories.
- :py:func:`~.create_sim_dir`: Creates the output directory.
- :py:func:`~.submit_simulation`: Calls Condor to submit the simulation.
- :py:func:`~.submit`: Creates parameter and submition files,
as well as output directory and submits the simulation.

"""

import os
import time
import numpy as np
import htcondor as htc
import ektome.globals as glb
import ektome.proektome.parfile as p
import ektome.proektome.create_sub_file as csf


def create_base_names(simulation):
    """Creates a name which serves as a base to name various
    output files and directories.

    :param simulation: Dictionary containing simulation info
    :type simulation: dict
    :returns: Base names for vanilla and excision case.
    :rtype: str

    """
    tmp = []
    for i in simulation:
        if (simulation[i] - int(simulation[i])) == 0.0:
            tmp.append(f"_{i}{int(simulation[i])}")
        else:
            tmp.append(f"_{i}{simulation[i]}")

    suffix = "".join(tmp)
    vnl_name = "vanilla{suffix}"
    exc_name = "excision{suffix}"
    return vnl_name, exc_name


def create_sim_dir(base_name):
    """Creates the output directory for the simulation within the
    proj_path/simulations/ directory.

    :param base_name: Base name to name the output directory.
    :returns: Directory path
    :rtype: str
    """
    dir_name = f"{glb.simulations_path}/{base_name}"

    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    else:
        pass
    return dir_name


def submit_simulation(sub_file_dict):
    """Calls Condor to submit a submition file.

    :param sub_file_path: Path to the submition file
    :type sub_file_path: str
    """
    job = htc.Submit(sub_file_dict)
    schedd = htc.Schedd()
    with schedd.transaction() as txn:
        try:
            cluster_id = job.queue(txn)
            time.sleep(0.1)
            return cluster_id
        except RuntimeError:
            print("==============================")
            print("Simulation Submittion FAILED")
            print(sub_file_dict)
            print("==============================")
            return 0


def norm(x, y, z):
    return np.sqrt(x * x + y * y + z * z)


def calculate_momentum(simulation):
    """Calculates the norm of the momentum vector for the simulation.

    :param simulation: Dictionary containing simulation info
    :type simulation: dict
    :returns: The norm value
    :rtype: np.float
    """
    px1 = simulation["px1"]
    py1 = simulation["py1"]
    pz1 = simulation["pz1"]

    px2 = simulation["px2"]
    py2 = simulation["py2"]
    pz2 = simulation["pz2"]
    return norm(px1, py1, pz1), norm(px2, py2, pz2)


def calculate_spin(simulation):
    """Calculates the norm of the effective spin for the simulation.

    :param simulation: Dictionary containing simulation info
    :type simulation: dict
    :returns: The norm value
    :rtype: np.float
    """
    sx1 = simulation["sx1"]
    sy1 = simulation["sy1"]
    sz1 = simulation["sz1"]

    sx2 = simulation["sx2"]
    sy2 = simulation["sy2"]
    sz2 = simulation["sz2"]
    return norm(sx1, sy1, sz1), norm(sx2, sy2, sz2)


def write_submit_metadata(cluster_id, base_name):
    """Writes the cluster_id and the simulation name into a
    metadata file.

    :param cluster_id: The cluster_id from HTCondor.
    :type cluster_id: int
    :param base_name: The name of the simulation.
    :type base_name: str
    """
    with open(glb.metadata_path, "a") as metadata:
        metadata.write(f"{base_name},{cluster_id}\n")


def submit(simulation):
    """Main submit function. Creates parameter and submition files
    creates output directory and submits a simulation through Condor.

    :param simulation: Dictionary containing simulation info
    :type simulation: dict
    """
    s1, s2 = calculate_spin(simulation)
    p1, p2 = calculate_momentum(simulation)

    if (s1 > 1) or (s2 > 1):
        print("Spin > 1")
        return 1
    if (p1 > 1) or (p2 > 1):
        print("Momentum > 1")
        return 1

    vanilla_base_name, excision_base_name = create_base_names(simulation)

    # Create parameter files
    vanilla_parfile = p.Parameter_File(simulation, vanilla_base_name, 1)
    vanilla_parfile.write_parfile()

    excision_parfile = p.Parameter_File(simulation, excision_base_name, 1)
    excision_parfile.write_parfile()

    # Create submition files
    vnl_sub = csf.create_sub_dict(vanilla_base_name)
    exc_sub = csf.create_sub_dict(excision_base_name)

    # Create the simulation output directories
    vanilla_sim_dir = create_sim_dir(vanilla_base_name)
    excision_sim_dir = create_sim_dir(excision_base_name)

    # Submit the simulations
    vnl_id = submit_simulation(vnl_sub)
    exc_id = submit_simulation(exc_sub)

    write_submit_metadata(vnl_id, vanilla_base_name)
    write_submit_metadata(exc_id, excision_base_name)
    return 0
