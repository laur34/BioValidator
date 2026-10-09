# BioValidator

A production-grade, memory-efficient CLI tool to validate FASTA file structures, verify if it's COI by performing an NCBI test BLAST of a subsample of sequences, and check consistency of fasta headers with an OTU table.

## Installation

After cloning or downloading the repository, it is simple to install.

1. Navigate to its root directory:
```bash
cd /path/to/Biovalidator
```

2. Install using `pip`:
```bash
pip install -e .
```

The optional "-e" switch is for **editable development mode** (allowing you to tweak code logic dynamically without re-installing).

## Pipeline Features & Architecture

The validation pipeline performs three sequential validation phases using a memory-optimized streaming approach:

1. **Structural Verification:** Streams the input file line-by-line (safely validating files >50GB without system memory spikes). It enforces strict FASTA header alignment rules, strips quality indicators to catch unintended FASTQ submissions, and checks characters against clean IUPAC code parameters (wildcards are forbidden).
2. **Biological Marker Alignment:** Pulls a subsample (default size of 10) of sequence strings spaced deterministically across your entire dataset. It groups them into one batched remote NCBI BLAST network request payload to protect against server rate limits, and uses keyword threshold filtering to verify if the file likely contains COI marker fragments.
3. **OTU Cross-Referencing (Optional):** Performs an O(N) symmetric set comparison checking your FASTA headers against an accompanying OTU table to guarantee that sequence IDs between files align.


## Usage Examples

Once installed, the global application tool is accessible natively from any terminal window using the `bio-validate` command name wrapper.

### 1. Basic FASTA Validation
Validates structural formatting and runs a 10-sequence remote BLAST check:
```bash
bio-validate /path/to/my_sequences.fasta
```

### 2. Adjusting the Biological Sampling Quota
Increases or decreases the number of internal cross-section reads dispatched to NCBI:
```bash
bio-validate --sample-size 5 /path/to/my_sequences.fasta
```

### 3. Integrated FASTA & OTU Matrix Cross-Validation
Validates the structural and biological markers, then ensures the identifiers match your OTU table perfectly. By default, the matrix delimiter matches tab-separated formats (`.tsv`):
```bash
bio-validate --otu-table /path/to/otu_table.tsv /path/to/my_sequences.fasta
```

## Running Automated Unit Tests 🧪

Automated testing layouts are built using the `pytest` testing ecosystem framework.
Pytest will automatically run all the tests in the test folder. test_validator.py uses Mocking to circumvent the need to connect to NCBI's liver servers for the BLAST.
To install and use pytest:

```bash
# Install pytest framework dependencies
pip install pytest

# Run the test execution runner suite
pytest
```
