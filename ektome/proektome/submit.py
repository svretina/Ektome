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

"""This module provides helper functions to submit a simulation."""

import logging
import time
from typing import Any, Dict

try:
    import htcondor as htc
except ImportError:
    htc = None  # type: ignore

import ektome.globals as glb
import ektome.proektome.create_sub_file as csf
import ektome.proektome.parfile as p
from ektome.exceptions import JobSubmissionError

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def create_sim_dir(base_name: str) -> str:
    """Creates the output directory for the simulation.

    Args:
        base_name: Base name for the output directory.

    Returns:
        The path to the created directory as a string.
    """
    dir_path = glb.SIMULATIONS_PATH / base_name
    if not dir_path.exists():
        logger.info(f"Creating simulation directory: {dir_path}")
        dir_path.mkdir(parents=True, exist_ok=True)

    return str(dir_path)


def submit_simulation(sub_file_dict: Dict[str, Any], max_retries: int = 5) -> int:
    """Calls HTCondor to submit a simulation job.

    Args:
        sub_file_dict: Dictionary containing HTCondor submission parameters.
        max_retries: Maximum number of submission attempts.

    Returns:
        The Cluster ID assigned by HTCondor.

    Raises:
        JobSubmissionError: If the job fails to submit after max_retries.
    """
    if htc is None:
        logger.error("HTCondor module not found. Cannot submit job.")
        raise JobSubmissionError("HTCondor module not found.")

    job = htc.Submit(sub_file_dict)
    schedd = htc.Schedd()
    
    for attempt in range(1, max_retries + 1):
        try:
            with schedd.transaction() as txn:
                cluster_id: int = job.queue(txn)
                logger.info(f"Successfully submitted job. Cluster ID: {cluster_id}")
                return cluster_id
        except Exception as exc:
            logger.warning(f"Submission attempt {attempt}/{max_retries} failed: {exc}")
            if attempt == max_retries:
                msg = f"Failed to submit job after {max_retries} attempts."
                raise JobSubmissionError(msg) from exc
            time.sleep(2**attempt)  # Exponential backoff

    return 0  # Should not reach here


def write_submit_metadata(cluster_id: int, base_name: str) -> None:
    """Writes the cluster_id and the simulation name into a metadata file.

    Args:
        cluster_id: The cluster ID from HTCondor.
        base_name: The name of the simulation.
    """
    try:
        with glb.METADATA_PATH.open("a") as metadata:
            metadata.write(f"{base_name},{cluster_id}\n")
    except OSError as exc:
        logger.error(f"Failed to write metadata for {base_name}: {exc}")


def get_info_from_folder_name(folder_name: str) -> Dict[str, Any]:
    """Parses simulation information from the folder name.

    Args:
        folder_name: The name of the simulation folder.

    Returns:
        A dictionary containing parsed parameters.
    """
    pieces = folder_name.split("_")
    try:
        d = {
            "q": int(pieces[1][1:]),
            "sx1": float(pieces[3][3:]),
            "sy1": float(pieces[4][3:]),
            "sz1": float(pieces[5][3:]),
            "sx2": float(pieces[6][3:]),
            "sy2": float(pieces[7][3:]),
            "sz2": float(pieces[8][3:]),
        }
        d["exr"] = d["q"] / 2
        return d
    except (IndexError, ValueError) as exc:
        logger.error(f"Failed to parse folder name '{folder_name}': {exc}")
        return {}


def check_existence_simulation(simulation_folder_name: str) -> bool:
    """Checks if a simulation already exists and has produced output.

    Args:
        simulation_folder_name: Name of the simulation folder.

    Returns:
        True if the simulation exists and is considered finished.
    """
    file_path = glb.SIMULATIONS_PATH / simulation_folder_name / "twopunctures.xyz.h5"
    return file_path.exists()


def submit(simulation: Any) -> int:
    """Main submission entry point.

    Creates parameter and submission files, sets up directories, and 
    submits the simulation to HTCondor.

    Args:
        simulation: An object or dictionary containing simulation information.

    Returns:
        Exit code (0 for success).
    """
    v_base = simulation.vanilla_base_name
    e_base = simulation.excision_base_name

    logger.info(f"Preparing submission for: {v_base.replace('vanilla_', '')}")
    logger.info("=" * 30)

    if not check_existence_simulation(v_base):
        logger.info(f"Processing vanilla simulation: {v_base}")
        vanilla_parfile = p.Parameter_File(simulation, excision=False, N=1)
        vanilla_parfile.write_parfile()
        
        v_sub = csf.create_sub_dict(v_base)
        create_sim_dir(v_base)
        
        try:
            v_id = submit_simulation(v_sub)
            write_submit_metadata(v_id, v_base)
        except JobSubmissionError as exc:
            logger.error(f"Vanilla submission failed: {exc}")

    logger.info(f"Processing excision simulation: {e_base}")
    excision_parfile = p.Parameter_File(simulation, excision=True, N=1)
    excision_parfile.write_parfile()

    e_sub = csf.create_sub_dict(e_base)
    create_sim_dir(e_base)

    try:
        e_id = submit_simulation(e_sub)
        write_submit_metadata(e_id, e_base)
    except JobSubmissionError as exc:
        logger.error(f"Excision submission failed: {exc}")

    return 0


def main() -> None:
    """CLI entry point for submission."""
    # This would normally parse arguments and call submit()
    pass

if __name__ == "__main__":
    main()
