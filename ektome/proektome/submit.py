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
from htcondor import HTCondorIOError as htc_error
import ektome.globals as glb
import ektome.proektome.parfile as p
import ektome.proektome.create_sub_file as csf


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
    :returns: Cluster ID number
    :rtype: int
    """
    job = htc.Submit(sub_file_dict)
    schedd = htc.Schedd()
    cluster_id = 0
    sucess = False
    while not sucess:
        try:
            with schedd.transaction() as txn:
                cluster_id = job.queue(txn)
                sucess = True
        except:
            sucess = False
            print("FAILED")
    return cluster_id


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

def get_info_from_folder_name(folder_name):
    # excision_q100_b639_sx1-0.1_sy1-0.9_sz1-0.1_sx20.1_sy20.1_sz20.9
    pieces = folder_name.split("_")
    d = {}
    d["q"] = int(pieces[1][1:])
    d["exr"] = int(d["q"]/2)
    d["sx1"] = float(pieces[3][3:])
    d["sy1"] = float(pieces[4][3:])
    d["sz1"] = float(pieces[5][3:])
    d["sx2"] = float(pieces[6][3:])
    d["sy2"] = float(pieces[7][3:])
    d["sz2"] = float(pieces[8][3:])
    return d


def check_output_dir(dir1):
    if not os.path.exists(dir1):
        create_sim_dir(dir1)

def check_existence_simulation(simulation_folder_name):
    """Checks if a simulation already exists and is finished.
    :param simulation_folder_name: Name of the simulation folder
    :type simulation_folder_name: str
    :returns: Boolean
    :rtype: bool
    """
    folder_path = f"{glb.simulations_path}/{simulation_folder_name}"
    file_path = f"{folder_path}/twopunctures.xyz.h5"
    if os.path.exists(file_path):
        return True
    else:
        return False


def submit(simulation):
    """Main submit function. Creates parameter and submition files
    creates output directory and submits a simulation through Condor.

    :param simulation: Dictionary containing simulation info
    :type simulation: dict
    """
    vanilla_base_name = simulation.vanilla_base_name
    excision_base_name = simulation.excision_base_name
    submit = True
    print("Submiting simulations for:")
    print(vanilla_base_name.split("vanilla_")[1])
    print(30*"==")

    if not check_existence_simulation(vanilla_base_name):
        vanilla_parfile = p.Parameter_File(simulation,
                                       excision=False, N=1)
        vanilla_parfile.write_parfile()
        vnl_sub = csf.create_sub_dict(vanilla_base_name)
        vanilla_sim_dir = create_sim_dir(vanilla_base_name)
        check_output_dir(vanilla_sim_dir)
        if submit:
            vnl_id = submit_simulation(vnl_sub)
            write_submit_metadata(vnl_id, vanilla_base_name)


    excision_parfile = p.Parameter_File(simulation,
                                        excision=True, N=1)

    excision_parfile.write_parfile()

    exc_sub = csf.create_sub_dict(excision_base_name)
    excision_sim_dir = create_sim_dir(excision_base_name)
    check_output_dir(excision_sim_dir)

    if submit:
        exc_id = submit_simulation(exc_sub)
        write_submit_metadata(exc_id, excision_base_name)

    return 0
