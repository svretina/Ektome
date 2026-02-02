import sys
from unittest.mock import MagicMock

# Define all external modules that should be mocked
EXTERNAL_MODULES = [
    "htcondor",
    "jhuki",
    "jhuki.twopunctures",
    "PyAstronomy",
    "PyAstronomy.pyasl",
    "kuibit",
    "pandas",
    "matplotlib",
    "matplotlib.pyplot",
]

for mod in EXTERNAL_MODULES:
    sys.modules[mod] = MagicMock()
