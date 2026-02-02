import pytest
import sys
from unittest.mock import MagicMock, patch
from pathlib import Path

# External dependencies are mocked in conftest.py
import htcondor as htc
mock_htc = htc

from ektome.proektome.submit import submit_simulation
from ektome.exceptions import JobSubmissionError

@pytest.fixture(autouse=True)
def reset_mock():
    mock_htc.reset_mock()

def test_submit_simulation_success():
    """Verify that submit_simulation calls HTCondor correctly on success."""
    # Resetting specifically for this test
    mock_htc.Submit.reset_mock()
    mock_htc.Schedd.reset_mock()
    
    mock_job = MagicMock()
    mock_schedd = MagicMock()
    
    mock_htc.Submit.return_value = mock_job
    mock_htc.Schedd.return_value = mock_schedd
    
    # Mock transaction context manager
    mock_txn = MagicMock()
    mock_schedd.transaction.return_value.__enter__.return_value = mock_txn
    mock_job.queue.return_value = 12345
    
    cluster_id = submit_simulation({"executable": "/bin/ls"})
    
    assert cluster_id == 12345
    mock_htc.Submit.assert_called_once()
    mock_job.queue.assert_called_once_with(mock_txn)

def test_submit_simulation_retry_logic():
    """Verify that submit_simulation retries on failure."""
    # Reset state
    mock_htc.Schedd.reset_mock()
    
    mock_job = MagicMock()
    mock_schedd = MagicMock()
    
    mock_htc.Submit.return_value = mock_job
    mock_htc.Schedd.return_value = mock_schedd
    
    # Simple retry case: fail once, then succeed
    mock_cm = MagicMock()
    mock_schedd.transaction.side_effect = [
        Exception("Fail"),
        mock_cm
    ]
    mock_cm.__enter__.return_value = MagicMock()
    mock_job.queue.return_value = 54321
    
    with patch("time.sleep"):
        cluster_id = submit_simulation({"executable": "/bin/ls"}, max_retries=5)
        assert cluster_id == 54321
        assert mock_schedd.transaction.call_count == 2

def test_submit_simulation_failure():
    """Verify that JobSubmissionError is raised after max retries."""
    mock_schedd = MagicMock()
    mock_htc.Schedd.return_value = mock_schedd
    mock_schedd.transaction.side_effect = Exception("Permanent failure")
    
    with patch("time.sleep"):
        with pytest.raises(JobSubmissionError):
            submit_simulation({"executable": "/bin/ls"}, max_retries=2)
