# src/bio_validator/cli.py
import sys
import argparse
from bio_validator.validator import FastaValidator
from bio_validator.exceptions import FastaValidationError, MarkerMismatchError

def main() -> None:
    parser = argparse.ArgumentParser(description="Production-grade FASTA validator.")
    parser.add_argument("fasta_file", type=str, help="Path to input FASTA file.")
    parser.add_argument("--sample-size", type=int, default=10)
    args = parser.parse_args()

    print(f"🧬 Initializing validation for: {args.fasta_file}")
    validator = FastaValidator(file_path=args.fasta_file, sample_size=args.sample_size)
    
    try:
        # Run validations — they will either succeed or throw an explicit exception
        validator.check_structural_integrity()
        print("✅ Structural integrity check passed.")
        
        validator.run_blast_check()
        print("🎉 Validation Complete: File is structurally sound and verified as COI!")
        sys.exit(0)

    # 🛡️ Defensive Catch Blocks
    except FileNotFoundError as e:
        print(f"❌ Input Error: {e}")
        sys.exit(2)
    except FastaValidationError as e:
        print(f"❌ Data Integrity Failure: {e}")
        sys.exit(1)
    except MarkerMismatchError as e:
        print(f"❌ Biological Verification Failure: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"💥 Unexpected system failure occurred: {e}")
        sys.exit(3)
