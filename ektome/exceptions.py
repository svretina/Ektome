"""Custom exception classes for the Ektome package."""

class EktomeError(Exception):
    """Base class for all Ektome exceptions."""
    pass

class ClusterError(EktomeError):
    """Base class for cluster-related failures."""
    pass

class JobSubmissionError(ClusterError):
    """Raised when an HTCondor job submission fails."""
    pass

class TransactionError(ClusterError):
    """Raised when an HTCondor transaction fails."""
    pass

class PostProcessingError(EktomeError):
    """Raised when post-processing of simulation data fails."""
    pass

class SimulationNotFoundError(PostProcessingError):
    """Raised when a specified simulation directory or file is not found."""
    pass

class ConfigurationError(EktomeError):
    """Raised when there is an error in the configuration file or its parsing."""
    pass
