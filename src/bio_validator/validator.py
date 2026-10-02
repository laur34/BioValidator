from pathlib import Path
import random, re
from bio_validator.exceptions import FastaValidationError
import urllib.error
from Bio import Blast  # Requires Biopython
from Bio.Blast import NCBIWWW, NCBIXML
from bio_validator.exceptions import FastaValidationError, MarkerMismatchError
import io  # Add this near your other imports at the top

class FastaValidator:
    """Encapsulates the business logic for validating FASTA file integrity and contents."""
    
    def __init__(self, file_path: str | Path, sample_size: int = 10):
        self.file_path = Path(file_path)
        self.sample_size = sample_size

    def validate_file_exists(self) -> None:
        """Defensive check before any parsing happens."""
        if not self.file_path.exists():
            raise FileNotFoundError(f"Target file missing: {self.file_path}")
        if not self.file_path.is_file():
            raise IsADirectoryError(f"Expected a file, but found a directory: {self.file_path}")

    def check_structural_integrity(self) -> None:
        """Streams the file to verify format. Raises FastaValidationError if corrupt.
        
        Returns:
            None. Execution completes smoothly if valid.
        """
        self.validate_file_exists()

        valid_sequence_pattern = re.compile(r"^[A-Z\-\s]+$", re.IGNORECASE)
        has_seen_header = False
        lines_checked = 0
        max_lines_to_sample = 50000  

        with open(self.file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                cleaned_line = line.strip()
                if not cleaned_line:
                    continue

                lines_checked += 1

                if not has_seen_header:
                    if cleaned_line.startswith(">"):
                        has_seen_header = True
                        continue
                    elif cleaned_line.startswith("@"):
                        raise FastaValidationError(
                            f"File starts with '@'. This appears to be a FASTQ file, not a FASTA."
                        )
                    else:
                        raise FastaValidationError(
                            f"First non-empty line does not start with '>'. Found: '{cleaned_line[:20]}...'"
                        )

                if cleaned_line.startswith(">"):
                    continue
                
                if not valid_sequence_pattern.match(cleaned_line):
                    raise FastaValidationError(
                        f"Invalid characters detected in sequence on line {lines_checked}: '{cleaned_line[:20]}...'"
                    )

                if lines_checked >= max_lines_to_sample:
                    break

        if not has_seen_header:
            raise FastaValidationError("The file appears to be completely empty or lacks valid headers.")
      
    def _count_sequences(self) -> int:
        """Pass 1: Fast stream to count total sequence headers."""
        self.validate_file_exists()
        count = 0
        with open(self.file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
        return count

    def subsample_sequences(self) -> list[dict[str, str]]:
        """Pass 2: Pulls exactly N evenly spaced sequences without loading the file into RAM.
        
        Returns:
            A list of dictionaries, where each dict has 'header' and 'sequence' keys.
        """
        total_records = self._count_sequences()
        
        # Defensive edge case: file has fewer sequences than our desired sample size
        if total_records <= self.sample_size:
            # If it's small, just grab everything
            return self._extract_all_sequences()

        # Calculate step size to space out our sample evenly across the file
        stride = total_records // self.sample_size
        
        # Determine the exact indices we want to save (e.g., [0, stride, stride*2, ...])
        target_indices = {i * stride for i in range(self.sample_size)}

        sampled_records = []
        current_record_index = -1
        current_header = ""
        current_seq_fragments = []

        with open(self.file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                cleaned_line = line.strip()
                if not cleaned_line:
                    continue

                if cleaned_line.startswith(">"):
                    # We hit a new header! Process the previous record first if we were saving it
                    if current_record_index in target_indices and current_header:
                        sampled_records.append({
                            "header": current_header,
                            "sequence": "".join(current_seq_fragments)
                        })
                        current_seq_fragments = []

                    current_record_index += 1
                    current_header = cleaned_line

                    # Quick optimization short-circuit: stop reading if we have filled our sample quota
                    if len(sampled_records) == self.sample_size:
                        break
                        
                else:
                    # It's a sequence line. Append it only if this record index is one we want
                    if current_record_index in target_indices:
                        current_seq_fragments.append(cleaned_line)

        return sampled_records

    def _extract_all_sequences(self) -> list[dict[str, str]]:
        """Fallback method to grab everything if the file is smaller than sample_size."""
        records = []
        current_header = ""
        current_seq_fragments = []

        with open(self.file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                cleaned_line = line.strip()
                if not cleaned_line:
                    continue
                if cleaned_line.startswith(">"):
                    if current_header:
                        records.append({"header": current_header, "sequence": "".join(current_seq_fragments)})
                        current_seq_fragments = []
                    current_header = cleaned_line
                else:
                    current_header and current_seq_fragments.append(cleaned_line)
                    
            if current_header:
                records.append({"header": current_header, "sequence": "".join(current_seq_fragments)})
                
        return records
           

    def run_blast_check(self) -> None:
        """Samples sequences, runs a remote BLAST query, and verifies they are COI.
        
        Raises:
            MarkerMismatchError: If the sequences don't align with Cytochrome C Oxidase I.
            RuntimeError: If a network or API communication error occurs.
        """
        # 1. Grab our 10 representative sequences from disk
        sampled_records = self.subsample_sequences()
        if not sampled_records:
            raise FastaValidationError("Cannot run BLAST validation; no sequences found in file.")

        # 2. Convert our list of dicts back into a single multi-FASTA string for the API call
        # This keeps us down to exactly ONE network request instead of ten separate hits.
        multi_fasta_query = ""
        for record in sampled_records:
            multi_fasta_query += f"{record['header']}\n{record['sequence']}\n"

        print(f"📡 Sending {len(sampled_records)} sequences to NCBI BLAST (this can take a minute)...")

        try:
            # 1. Set the global module identification properties. This satisfies NCBI guidelines.
            NCBIWWW.email = "lv70xo@gmail.com"
            # 2. Execute the remote BLAST search safely using ONLY strict API parameters
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=multi_fasta_query,
            )
            
            # Read the XML response object from the network stream
            blast_results_raw = result_handle.read()
            result_handle.close()

        except urllib.error.URLError as e:
            raise RuntimeError(f"Network connection to NCBI BLAST failed: {e}")

            # Any other API or stream parsing failure
            raise RuntimeError(f"An error occurred while communicating with NCBI: {e}")

        # 4. Parse and evaluate the biological content
        # Compile a regex to look for variations of COI / Cytochrome Oxidase Subunit 1
        coi_keywords = re.compile(r"(coi|cox1|cytochrome\s+c\s+oxidase\s+subunit\s+1)", re.IGNORECASE)
        
        # FIX: Wrap the raw string directly in io.StringIO
        blast_records = NCBIXML.parse(io.StringIO(blast_results_raw))
        
        confirmed_coi_count = 0
        total_queries_evaluated = 0

        for blast_record in blast_records:
            total_queries_evaluated += 1
            
            # Defensive check: Did this sequence hit anything at all?
            if not blast_record.alignments:
                continue  # Dark matter sequence, or sequence too short to hit anything

            # Snag the definition line of the #1 absolute top hit for this sequence
            top_hit_title = blast_record.alignments[0].title # Added indexing protection [0]
            
            if coi_keywords.search(top_hit_title):
                confirmed_coi_count += 1

        # 5. Make an architectural decision based on the sample cross-section
        # If less than 70% of hits match COI, reject the file as containing non-COI sequences.
        if total_queries_evaluated == 0 or (confirmed_coi_count / total_queries_evaluated) < 0.7:
            raise MarkerMismatchError(
                f"Biological verification failed. Only {confirmed_coi_count}/{total_queries_evaluated} "
                f"sampled sequences matched the COI marker database entries."
            )

        # Smooth exit implies success!
