# tests/test_validator.py
import pytest
import io
from pathlib import Path
from unittest.mock import patch
from bio_validator.validator import FastaValidator
from bio_validator.exceptions import FastaValidationError, MarkerMismatchError

# =====================================================================
# 1. STRUCTURAL INTEGRITY TESTS (Using pytest's built-in tmp_path)
# =====================================================================

def test_valid_fasta_passes_structural_check(tmp_path):
    """Verifies that a perfectly structured FASTA file passes without errors."""
    # Create a temporary file safely using the tmp_path fixture
    fasta_file = tmp_path / "good.fasta"
    fasta_file.write_text(">seq1\nACTGATCGATCG\n>seq2\nATCGATCGATCG\n")
    
    validator = FastaValidator(file_path=fasta_file)
    
    # Execution should run cleanly to completion (returns None)
    assert validator.check_structural_integrity() is None


def test_fastq_input_raises_fasta_validation_error(tmp_path):
    """Verifies that providing a FASTQ file format triggers our custom structural exception."""
    fastq_file = tmp_path / "bad.fastq"
    fastq_file.write_text("@ERR12345\nACTGATCGATCG\n+\nIIIIIIIIIIII\n")
    
    validator = FastaValidator(file_path=fastq_file)
    
    # Assert that this block explicitly raises a FastaValidationError
    with pytest.raises(FastaValidationError) as exc_info:
        validator.check_structural_integrity()
        
    assert "appears to be a FASTQ file" in str(exc_info.value)


def test_invalid_characters_raise_error(tmp_path):
    """Verifies that forbidden non-alphabetic/wildcard characters are safely rejected."""
    corrupt_file = tmp_path / "corrupt.fasta"
    corrupt_file.write_text(">seq1\nACTG*TCG123\n")  # Contains illegal characters '*' and numbers
    
    validator = FastaValidator(file_path=corrupt_file)
    
    with pytest.raises(FastaValidationError):
        validator.check_structural_integrity()


# =====================================================================
# 2. REMOTE BIOLOGICAL BLAST MOCKING TESTS
# =====================================================================

# Mock Data: A sample fake XML string that mirrors what NCBI returns for a good COI alignment
MOCK_NCBI_XML_PASS = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_id>gnl|BL_ORD_ID|0</Hit_id>
          <Hit_def>Homo sapiens cytochrome c oxidase subunit 1 (COI) mRNA</Hit_def>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

@patch("bio_validator.validator.NCBIWWW.qblast")
def test_blast_check_passes_with_coi_hits(mock_qblast, tmp_path):
    """Intercepts the internet call and feeds our validator mock COI database responses."""
    # 1. Setup a valid sample file so the subsampler functions correctly
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(">seq1\nACTGATCGATCG\n")
    
    # 2. Configure our internet mock to instantly return our fake passing XML string
    mock_handle = io.StringIO(MOCK_NCBI_XML_PASS)
    mock_qblast.return_value = mock_handle
    
    validator = FastaValidator(file_path=fasta_file, sample_size=1)
    
    # 3. Execute. The method will run cleanly without hitting the real internet!
    assert validator.run_blast_check() is None
    
    # Double check that our internal machinery actually tried to trigger the API wrapper
    mock_qblast.assert_called_once()

