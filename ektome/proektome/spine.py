"""This module provides helper functions to manage submission of simulations.
The functions available are:
- :py:func:`~.read_config`: Reads a config file and outputs a dictionary with the values of the config file.
- :py:func:`~.create_dirs`: Checks if output directories exists and if not, creates them.
- :py:func:`~.get_ini_file`: Searches the directory for a .ini config file and returns its path.

"""

import configparser as cfg
import itertools
import logging
from typing import Any, Dict, List, Optional

import numpy as np

import ektome.globals as glb
import ektome.proektome.binary as bnr
import ektome.proektome.simulation as simulation
import ektome.proektome.submit as sb
from ektome.exceptions import ConfigurationError

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def read_config(path: Optional[glb.Path] = None) -> Dict[str, Any]:
    """Reads a config file and returns a dictionary with its contents.

    Args:
        path: Path to the config file. If None, uses default.

    Returns:
        Dictionary with the values of the config file.

    Raises:
        ConfigurationError: If the config file cannot be read or is missing.
    """
    config_path = path or glb.CONFIG_PATH
    if not config_path.exists():
        ini_file = get_ini_file()
        if ini_file:
            config_path = glb.PROJECT_PATH / ini_file
        else:
            raise ConfigurationError(f"Configuration file not found: {config_path}")

    logger.info(f"Reading configuration from: {config_path}")
    parser = cfg.ConfigParser()
    try:
        files_read = parser.read(config_path)
        if not files_read:
            raise ConfigurationError(f"Failed to read config file: {config_path}")
    except Exception as exc:
        raise ConfigurationError(f"Error parsing config file: {exc}")

    config: Dict[str, Any] = {}
    for section in parser.sections():
        config[section] = dict(parser[section])
    return config


def create_dirs() -> None:
    """Ensures all required output directories exist."""
    logger.info("Initializing project directories...")
    for directory in glb.REQUIRED_DIRECTORIES:
        if not directory.exists():
            logger.debug(f"Creating directory: {directory}")
            directory.mkdir(parents=True, exist_ok=True)


def get_ini_file() -> Optional[str]:
    """Searches for a .ini config file and returns its name.

    Returns:
        Name of the .ini config file or None.
    """
    for item in glb.PROJECT_PATH.iterdir():
        if item.is_file() and item.suffix == ".ini":
            return item.name
    return None


def clear_submit_metadata() -> None:
    """Removes the submission metadata file if it exists."""
    if glb.METADATA_PATH.exists():
        logger.info(f"Clearing metadata file: {glb.METADATA_PATH}")
        glb.METADATA_PATH.unlink()


def my_arange(start: float, step: float, end: float) -> np.ndarray:
    """Custom arange for consistent floating point steps.

    Args:
        start: Start value.
        step: Step size.
        end: End value.

    Returns:
        Numpy array of values.
    """
    num = int(round((end - start) / step)) + 1
    tmp = np.linspace(start, end, num)
    return np.round(tmp, 3)


def create_config_arrays(cfg_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Creates numpy arrays spaced according to the steps provided in the config file.

    Args:
        cfg_dict: Dictionary containing the config information.

    Returns:
        Dictionary with the numpy arrays.
    """
    config_array: Dict[str, Any] = {}
    for section, content in cfg_dict.items():
        try:
            if "array" not in content:
                config_array[section] = np.nan
                continue

            temp_array = content["array"].split(",")
            values: List[np.ndarray] = []
            for item in temp_array:
                if ":" in item:
                    start, step, end = [float(j) for j in item.split(":")]
                    values.append(my_arange(start, step, end))
                else:
                    values.append(np.array([float(item)]))
            config_array[section] = np.sort(np.concatenate(values))
        except (ValueError, KeyError) as exc:
            logger.warning(f"Failed to parse section {section}: {exc}")
            config_array[section] = np.nan
    return config_array


def create_simulation_dict_and_submit(cfg_arr: Dict[str, Any]) -> None:
    """Creates simulate and submits simulations based on config arrays.

    Args:
        cfg_arr: Dictionary with the config info.
    """
    # Required keys check
    required_keys = [
        "minus_spin_x", "minus_spin_y", "minus_spin_z",
        "plus_spin_x", "plus_spin_y", "plus_spin_z",
        "mass_ratio", "number_of_orbits"
    ]
    for key in required_keys:
        if key not in cfg_arr or (isinstance(cfg_arr[key], float) and np.isnan(cfg_arr[key])):
            logger.error(f"Missing or invalid configuration for: {key}")
            return

    spin1 = list(itertools.product(
        cfg_arr["minus_spin_x"],
        cfg_arr["minus_spin_y"],
        cfg_arr["minus_spin_z"],
    ))

    spin2 = list(itertools.product(
        cfg_arr["plus_spin_x"],
        cfg_arr["plus_spin_y"],
        cfg_arr["plus_spin_z"],
    ))

    logger.info(f"Starting submission loop for {len(cfg_arr['mass_ratio'])} mass ratios...")
    for q in cfg_arr["mass_ratio"]:
        exr = q / 2.0
        binary = bnr.Binary(q)
        for s1 in spin1:
            if np.linalg.norm(s1) >= 1:
                logger.warning(f"S1 magnitude {np.linalg.norm(s1)} >= 1, skipping.")
                continue
            for s2 in spin2:
                if np.linalg.norm(s2) >= 1:
                    logger.warning(f"S2 magnitude {np.linalg.norm(s2)} >= 1, skipping.")
                    continue
                for n_orb in cfg_arr["number_of_orbits"]:
                    b = binary.semimajor(n_orb)
                    p1, p2 = binary.quasicircular_inspiral(q, 2 * b, s1, s2)
                    sim = simulation.Simulation(
                        q=q, b=b,
                        px1=p1[0], py1=p1[1], pz1=p1[2],
                        sx1=s1[0], sy1=s1[1], sz1=s1[2],
                        px2=p2[0], py2=p2[1], pz2=p2[2],
                        sx2=s2[0], sy2=s2[1], sz2=s2[2],
                        exr=exr,
                    )
                    sb.submit(sim)


def main() -> None:
    """Main execution point."""
    try:
        create_dirs()
        config_dict = read_config()
        config_arr = create_config_arrays(config_dict)
        create_simulation_dict_and_submit(config_arr)
    except Exception as exc:
        logger.exception(f"Fatal error during execution: {exc}")


if __name__ == "__main__":
    main()
