# src/bio_validator/exceptions.py

class BioValidatorError(Exception):
    """Base exception for all errors thrown by the bio_validator package."""
    pass

class FastaValidationError(BioValidatorError):
    """Raised when a file breaks structural FASTA conventions."""
    pass

class MarkerMismatchError(BioValidatorError):
    """Raised when the biological contents do not match the expected marker (e.g., not COI)."""
    pass

