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

class OtuMappingError(BioValidatorError):
    """Raised when and OTU table structure is malformed or mismatches its accompanying FASTA file."""
    pass