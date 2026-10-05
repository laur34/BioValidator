# BioValidator
This is a command-line program to validate whether a file is a valid fasta file, and if it is, whether its sequences are COI, by performing a small test BLAST.

## Installation
To install it (UNIX/LINUX), simply cd into the BioValidator directory once you have downloaded it, and run  ```pip install -e . ```

## Usage
```bio-validate [-h] [--sample-size SAMPLE_SIZE] fasta_file ```

Default sample size is 10. It is the number of sequences to sample from the fasta file, for a test BLAST.
You must supply a fasta file, as the input.