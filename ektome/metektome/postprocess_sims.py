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

"""This module provides functions to post process the simulation data."""

import logging
from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional

import pandas as pd

import ektome.globals as glb
import ektome.metektome.error as e

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def get_processed_sims(folders: List[str]) -> List[str]:
    """Identifies folders that have already been processed.

    Args:
        folders: List of folder names to check.

    Returns:
        List of folders containing 'psi.png'.
    """
    processed = []
    for folder in folders:
        if (glb.SIMULATIONS_PATH / folder / "psi.png").exists():
            processed.append(folder)
    return processed


def get_unfinished_sims(folders: List[str]) -> List[str]:
    """Identifies folders where the simulation is not yet finished.

    Args:
        folders: List of folder names to check.

    Returns:
        List of folders missing 'TwoPunctures.bbh'.
    """
    unfinished = []
    for folder in folders:
        if not (glb.SIMULATIONS_PATH / folder / "TwoPunctures.bbh").exists():
            unfinished.append(folder)
    return unfinished


def get_folders(recalc: bool = False) -> List[str]:
    """Gets the list of vanilla simulation folders that need processing.

    Args:
        recalc: If True, includes already processed folders.

    Returns:
        A list of folder names.
    """
    logger.info("Scanning simulations directory...")
    if not glb.SIMULATIONS_PATH.exists():
        logger.warning(f"Simulations path {glb.SIMULATIONS_PATH} does not exist.")
        return []

    sim_folders = [f.name for f in glb.SIMULATIONS_PATH.iterdir() if f.is_dir()]
    vanilla_folders = [x for x in sim_folders if x.startswith("vanilla")]

    unfinished = get_unfinished_sims(vanilla_folders)
    
    if recalc:
        # Exclude only unfinished
        folders = [f for f in vanilla_folders if f not in unfinished]
    else:
        # Exclude processed AND unfinished
        processed = get_processed_sims(vanilla_folders)
        exclude = set(processed) | set(unfinished)
        folders = [f for f in vanilla_folders if f not in exclude]

    logger.info(f"Found {len(folders)} folders to process.")
    return folders


def get_error_dataframe(folder: str) -> pd.DataFrame:
    """Calculates error information for a single simulation folder.

    Args:
        folder: The folder name.

    Returns:
        A DataFrame containing the error report.
    """
    logger.debug(f"Processing folder: {folder}")
    try:
        err = e.Error(folder, dim=3)
        return err.error_report()
    except Exception as exc:
        logger.error(f"Failed to process {folder}: {exc}")
        return pd.DataFrame()


def run_postprocessing(
    save: bool = True, recalc: bool = False, n_cpus: Optional[int] = None
) -> pd.DataFrame:
    """Runs the post-processing loop over folders, optionally in parallel.

    Args:
        save: Whether to save the results to CSV.
        recalc: Whether to re-process folders.
        n_cpus: Number of CPUs to use for parallel processing.

    Returns:
        A concatenated DataFrame of all results.
    """
    folders = get_folders(recalc)
    if not folders:
        return pd.DataFrame()

    logger.info(f"Starting post-processing with {n_cpus or 'default'} workers...")
    
    results = []
    with ProcessPoolExecutor(max_workers=n_cpus) as executor:
        results = list(executor.map(get_error_dataframe, folders))

    if not results:
        return pd.DataFrame()

    error_info = pd.concat(results, ignore_index=True)
    
    if save:
        output_path = glb.RESULTS_PATH / "error_data3D.csv"
        glb.RESULTS_PATH.mkdir(parents=True, exist_ok=True)
        
        mode = "w" if recalc else "a"
        header = recalc or not output_path.exists()
        
        error_info.to_csv(output_path, mode=mode, header=header, index=False)
        logger.info(f"Results saved to {output_path}")

    return error_info


def main() -> None:
    """CLI entry point for post-processing."""
    run_postprocessing()


if __name__ == "__main__":
    main()
