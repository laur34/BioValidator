# src/bio_validator/cli.py
import sys
import argparse
from bio_validator.validator import FastaValidator
from bio_validator.exceptions import FastaValidationError, MarkerMismatchError, OtuMappingError

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Production-grade FASTA structural, marker, and OTU table validator."
    )
    
    # Required positional argument
    parser.add_argument(
        "fasta_file", 
        type=str, 
        help="Path to the input FASTA file to validate."
    )
    
    # Optional sampling argument
    parser.add_argument(
        "--sample-size", 
        type=int, 
        default=10,
        help="Number of sequences to subsample for BLAST validation (default: 10)."
    )
    
    # New Optional OTU table validation argument
    parser.add_argument(
        "--otu-table", 
        type=str, 
        default=None,
        help="Optional path to a companion OTU table (.tsv or .csv) to cross-validate IDs."
    )

    args = parser.parse_args()

    print(f"🧬 Initializing validation pipeline for: {args.fasta_file}")
    validator = FastaValidator(file_path=args.fasta_file, sample_size=args.sample_size)
    
    try:
        # Step 1: Run the structural integrity validation
        validator.check_structural_integrity()
        print("✅ Structural integrity check passed.")
        
        # Step 2: Run the biological taxonomy validation
        validator.run_blast_check()
        print("✅ Biological marker verification passed (Verified as COI).")
        
        # Step 3: Run the optional OTU table cross-reference if the user provided one
        if args.otu_table:
            print(f"📊 Companion OTU table provided. Cross-validating: {args.otu_table}")
            validator.validate_otu_table(otu_table_path=args.otu_table)
            print("✅ OTU table mapping validation passed.")
            
        print("🎉 Pipeline Complete: All submitted datasets are fully validated and aligned!")
        sys.exit(0)

    # 🛡️ Clean, Granular Defensive Catch Blocks
    except FileNotFoundError as e:
        print(f"❌ Input Error: {e}")
        sys.exit(2)
    except FastaValidationError as e:
        print(f"❌ Data Integrity Failure: {e}")
        sys.exit(1)
    except MarkerMismatchError as e:
        print(f"❌ Biological Verification Failure: {e}")
        sys.exit(1)
    except OtuMappingError as e:
        print(f"❌ OTU Mapping Failure: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"💥 Unexpected system failure occurred: {e}")
        sys.exit(3)

if __name__ == "__main__":
    main()
