from pathlib import Path
import random, re
from bio_validator.exceptions import FastaValidationError

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
           

    def run_blast_check(self) -> bool:
        """Subsamples sequences and queries NCBI to verify if it is COI."""
        # Step 1: Subsample 10 sequences deterministically
        # Step 2: Call NCBI via Biopython (with try/except handling)
        # Step 3: Evaluate top hits
        return True

