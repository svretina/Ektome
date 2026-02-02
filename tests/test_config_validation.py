import pytest
import numpy as np
from pathlib import Path
from ektome.proektome.spine import create_config_arrays

def test_create_config_arrays_parsing():
    """Test that config strings are correctly converted to numpy arrays."""
    cfg = {
        "mass_ratio": {"array": "1.0, 2.0, 4.0"},
        "plus_spin_z": {"array": "0.1:0.1:0.3"}
    }
    
    arrays = create_config_arrays(cfg)
    
    assert len(arrays["mass_ratio"]) == 3
    assert 1.0 in arrays["mass_ratio"]
    assert 4.0 in arrays["mass_ratio"]
    
    # Test arange logic
    assert len(arrays["plus_spin_z"]) == 3 # 0.1, 0.2, 0.3
    assert arrays["plus_spin_z"][1] == 0.2

def test_create_config_arrays_missing_key():
    """Test handling of missing 'array' key in section."""
    cfg = {
        "invalid_section": {"not_an_array": "foo"}
    }
    import numpy as np
    arrays = create_config_arrays(cfg)
    assert np.isnan(arrays["invalid_section"])

def test_create_config_arrays_malformed():
    """Test handling of malformed array strings."""
    cfg = {
        "bad_section": {"array": "1.0:foo:2.0"}
    }
    import numpy as np
    arrays = create_config_arrays(cfg)
    assert np.isnan(arrays["bad_section"])
