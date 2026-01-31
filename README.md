# 🌌 Ektome: Modern Simulation Management for HTCondor

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)

**Ektome** (from Greek *εκτομή*, "excision") is a high-performance Python toolkit designed to streamline the lifecycle of Intermediate Mass Ratio Inspiral (IMRI) simulations. It provides a robust, type-safe interface for preparing, submitting, and post-processing EinsteinToolkit simulations on HTCondor clusters.

---

## 🚀 Key Features

- **⚡ Automated Submission**: Seamlessly handle HTCondor transactions and parallel job queuing.
- **🔪 Excision Logic**: Intelligent preparation of parameter files for black hole excision techniques.
- **📊 Advanced Post-Processing**: Parallel processing of simulation results with masked error analysis.
- **🛠️ Modern DX**: Built with `uv`, `ruff`, and `mypy` for a top-tier development experience.
- **🛡️ Robustness**: Exponential backoff for cluster transactions and custom exception handling.

---

## 📦 Installation

Experience the speed of `uv` to install Ektome and its dependencies:

```bash
# Clone the repository
git clone https://github.com/svretina/Ektome.git
cd Ektome

# Install dependencies and setup environment
uv sync
```

---

## 🚦 Quick Start

### Submitting a Simulation

```python
from ektome.proektome.submit import submit
from my_sim_config import SimulationConfig

# Initialize your simulation parameters
sim = SimulationConfig(q=100, b=6.4, excision=True)

# Generate parfiles and submit to HTCondor
submit(sim)
```

### Post-Processing Results

```python
from ektome.metektome.postprocess_sims import run_postprocessing

# Process all new simulations and save results
results = run_postprocessing(n_cpus=8)

print(results.head())
```

---

## 🛠️ Development

We use modern Python tooling to ensure code quality:

```bash
# Run linting
uv run ruff check .

# Run type checks
uv run mypy .

# Run tests
uv run pytest
```

---

## 📜 License

This project is licensed under the **GNU GPLv3**. See the [LICENSE](LICENSE) file for details.
