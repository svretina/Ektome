import pytest
import sys
from unittest.mock import MagicMock, patch

# Mock specialized dependencies before importing the module that uses them
sys.modules["htcondor"] = MagicMock()
sys.modules["jhuki"] = MagicMock()
sys.modules["jhuki.twopunctures"] = MagicMock()
sys.modules["jhuki.grid"] = MagicMock()
sys.modules["kuibit"] = MagicMock()
sys.modules["kuibit.simdir"] = MagicMock()

from ektome.proektome.submit import submit_simulation, htc
from ektome.exceptions import JobSubmissionError

def test_submit_simulation_success():
    mock_htc_submit = MagicMock()
    mock_schedd = MagicMock()
    
    with patch("htcondor.Submit", return_value=mock_htc_submit), \
         patch("htcondor.Schedd", return_value=mock_schedd), \
         patch("ektome.proektome.submit.htc") as patched_htc:
        
        patched_htc.Submit.return_value = mock_htc_submit
        patched_htc.Schedd.return_value = mock_schedd
        
        mock_txn = MagicMock()
        mock_schedd.transaction.return_value.__enter__.return_value = mock_txn
        mock_htc_submit.queue.return_value = 12345
        
        cluster_id = submit_simulation({"executable": "test"})
        assert cluster_id == 12345

def test_submit_simulation_retry_success():
    mock_htc_submit = MagicMock()
    mock_schedd = MagicMock()
    
    with patch("htcondor.Submit", return_value=mock_htc_submit), \
         patch("htcondor.Schedd", return_value=mock_schedd), \
         patch("ektome.proektome.submit.htc") as patched_htc, \
         patch("time.sleep"):
        
        patched_htc.Submit.return_value = mock_htc_submit
        patched_htc.Schedd.return_value = mock_schedd
        
        mock_txn = MagicMock()
        mock_schedd.transaction.return_value.__enter__.return_value = mock_txn
        
        # Fail once, then succeed
        mock_htc_submit.queue.side_effect = [Exception("Cluster full"), 67890]
        
        cluster_id = submit_simulation({"executable": "test"}, max_retries=2)
        assert cluster_id == 67890

def test_submit_simulation_failure():
    mock_htc_submit = MagicMock()
    mock_schedd = MagicMock()
    
    with patch("htcondor.Submit", return_value=mock_htc_submit), \
         patch("htcondor.Schedd", return_value=mock_schedd), \
         patch("ektome.proektome.submit.htc") as patched_htc, \
         patch("time.sleep"):
        
        patched_htc.Submit.return_value = mock_htc_submit
        patched_htc.Schedd.return_value = mock_schedd
        
        mock_txn = MagicMock()
        mock_schedd.transaction.return_value.__enter__.return_value = mock_txn
        mock_htc_submit.queue.side_effect = Exception("Permanent error")
        
        with pytest.raises(JobSubmissionError):
            submit_simulation({"executable": "test"}, max_retries=3)

def test_submit_no_htcondor_module():
    with patch("ektome.proektome.submit.htc", None):
        with pytest.raises(JobSubmissionError, match="HTCondor module not found"):
            submit_simulation({"executable": "test"})
