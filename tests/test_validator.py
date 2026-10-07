# tests/test_validator.py
import pytest
import io
from pathlib import Path
from unittest.mock import patch
from bio_validator.validator import FastaValidator
from bio_validator.exceptions import FastaValidationError

# =====================================================================
# 1. STRUCTURAL INTEGRITY TESTS
# =====================================================================

def test_valid_fasta_passes_structural_check(tmp_path):
    """Verifies that a perfectly structured FASTA file passes without errors."""
    fasta_file = tmp_path / "good.fasta"
    fasta_file.write_text(">seq1\nACTGATCGATCG\n>seq2\nATCGATCGATCG\n")
    
    validator = FastaValidator(file_path=fasta_file)
    # Smooth execution to completion ensures a clean pass
    assert validator.check_structural_integrity() is None


def test_fastq_input_raises_fasta_validation_error(tmp_path):
    """Verifies that providing a FASTQ file format triggers our custom structural exception."""
    fastq_file = tmp_path / "bad.fastq"
    fastq_file.write_text("@ERR12345\nACTGATCGATCG\n+\nIIIIIIIIIIII\n")
    
    validator = FastaValidator(file_path=fastq_file)
    
    with pytest.raises(FastaValidationError) as exc_info:
        validator.check_structural_integrity()
        
    assert "appears to be a FASTQ file" in str(exc_info.value)


def test_invalid_characters_raise_error(tmp_path):
    """Verifies that forbidden non-alphabetic/wildcard characters are safely rejected."""
    corrupt_file = tmp_path / "corrupt.fasta"
    corrupt_file.write_text(">seq1\nACTG*TCG123\n")
    
    validator = FastaValidator(file_path=corrupt_file)
    
    with pytest.raises(FastaValidationError):
        validator.check_structural_integrity()


# =====================================================================
# 2. REMOTE BIOLOGICAL BLAST MOCKING TESTS
# =====================================================================

@patch("bio_validator.validator.NCBIXML.parse")
@patch("bio_validator.validator.NCBIWWW.qblast")
def test_blast_check_passes_with_coi_hits(mock_qblast, mock_parse, tmp_path):
    """Intercepts both the network and parsing layer to cleanly verify validation flow."""
    # 1. Setup a valid sample file so the local subsampler runs cleanly
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(">seq1\nACTGATCGATCG\n")
    
    # 2. FIX: Wrap the string in io.StringIO so it has a .read() method!
    mock_qblast.return_value = io.StringIO("fake_xml_data")
    
    # 3. Create mock record objects that perfectly emulate Biopython's structure
    from unittest.mock import MagicMock
    mock_record = MagicMock()
    mock_alignment = MagicMock()
    
    # Set the top hit title directly to match our COI keyword logic pattern
    mock_alignment.title = "Homo sapiens cytochrome c oxidase subunit 1 (COI) mRNA"
    mock_record.alignments = [mock_alignment]
    
    # Make the mock parser return our fake records list
    mock_parse.return_value = [mock_record]
    
    validator = FastaValidator(file_path=fasta_file, sample_size=1)
    
    # 4. Execute the pipeline block
    assert validator.run_blast_check() is None
    
    # Ensure our internal modules were appropriately called by the manager method
    mock_qblast.assert_called_once()
    mock_parse.assert_called_once()
