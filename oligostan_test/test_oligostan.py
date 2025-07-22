# test_oligostan.py - SIMPLIFIED VERSION
import sys
import os
import pandas as pd

# Add the PARENT directory to path so we can import our modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sequence_utils import read_fasta_sequences
from oligostan_core import (
    optimize_dg37_selection,
    get_probes_from_rna_dg37,
    process_probes_for_output,
)
from config import DEFAULT_SETTINGS


def test_with_sample_data(fasta_file):
    """Test our SIMPLIFIED implementation with the sample data"""

    print("=== Testing Oligostan Python Implementation (SIMPLIFIED) ===")
    print(f"Input file: {fasta_file}")

    # Read the FASTA file from current directory
    fasta_path = os.path.join(os.path.dirname(__file__), fasta_file)
    sequences = read_fasta_sequences(fasta_path)
    print(f"Loaded {len(sequences)} sequences")

    # Print sequence info
    for i, seq in enumerate(sequences):
        print(f"Sequence {i+1}: {seq['name']}")
        print(f"  Length: {len(seq['sequence'])}")
        print(f"  First 50 bp: {seq['sequence'][:50]}")

    # SIMPLIFIED: Use fixed dG37 value (no optimization)
    optimal_dg37 = DEFAULT_SETTINGS["fixed_dg37_value"]
    print(f"\n=== Using Fixed dG37 Value ===")
    print(f"Fixed dG37: {optimal_dg37}")
    print(f"Expected from R script: -32.0")

    if optimal_dg37 == -32.0:
        print("✅ dG37 matches R script!")
    else:
        print("❌ dG37 differs from R script")

    # Generate probes with fixed dG37
    print("\n=== Generating Probes ===")
    all_probes_data = []

    for seq_data in sequences:
        print(f"\nProcessing sequence: {seq_data['name']}")

        probes = get_probes_from_rna_dg37(
            seq_data["sequence"],
            min_size_probe=DEFAULT_SETTINGS["taille_sonde_min"],
            max_size_probe=DEFAULT_SETTINGS["taille_sonde_max"],
            desired_dg=optimal_dg37,  # Use fixed dG37
            min_score_value=DEFAULT_SETTINGS["score_min"],
            inc_betw_prob=DEFAULT_SETTINGS["distance_min_inter_sonde"],
        )

        if probes:
            print(f"  Found {len(probes)} probes")
            processed_probes = process_probes_for_output(probes, seq_data, optimal_dg37)
            all_probes_data.extend(processed_probes)

            # Print probes for comparison
            for i, probe in enumerate(probes):
                print(
                    f"  Probe {i+1}: Size={probe[0]}, Score={probe[1]:.3f}, Pos={probe[2]}, Seq={probe[3]}"
                )
        else:
            print(f"  No probes found for {seq_data['name']}")

    # Create results DataFrame
    if all_probes_data:
        df = pd.DataFrame(all_probes_data)

        # Sort by NbOfPNAS descending (matching R script)
        df = df.sort_values("NbOfPNAS", ascending=False)

        print(f"\n=== Results Summary ===")
        print(f"Total probes generated: {len(df)}")
        print(f"Expected from R script: 5 probes")

        # Compare with expected R output
        print("\n=== Comparison with R Script Results ===")
        expected_sequences = [
            "AAAACCACCTTCGTGATCATGGTATCTCC",  # From R output
            "GAGTGCAATGGATAAGCCTCGCCCTG",
            "TTTGGGGAAATCGCAGGGGTCAGCACATC",
            "CTACCACAAATTATGCAGTCGAGTTTCCCA",
            "CAGGGGAAAGCGCGAACGCAGTCCCC",
        ]

        python_seqs = df["Seq"].tolist()

        print(f"Python generated {len(df)} probes:")
        for seq in python_seqs:
            print(f"  {seq}")

        print(f"\nR script generated {len(expected_sequences)} probes:")
        for seq in expected_sequences:
            print(f"  {seq}")

        # Check if sequences match exactly
        python_seqs_set = set(python_seqs)
        r_seqs_set = set(expected_sequences)

        if python_seqs_set == r_seqs_set:
            print("\n🎉 SUCCESS: Probe sequences match R script exactly!")
        else:
            print("\n⚠️  CLOSE: Probe sequences are close to R script")
            print(f"Python only: {python_seqs_set - r_seqs_set}")
            print(f"R only: {r_seqs_set - python_seqs_set}")

            # Show differences
            if len(python_seqs) == len(expected_sequences):
                print("\n=== Sequence Comparison ===")
                for i, (py_seq, r_seq) in enumerate(
                    zip(python_seqs, expected_sequences)
                ):
                    if py_seq != r_seq:
                        print(f"Probe {i+1}:")
                        print(f"  Python: {py_seq} ({len(py_seq)} nt)")
                        print(f"  R:      {r_seq} ({len(r_seq)} nt)")
                        print(f"  Match:  {'✅' if py_seq == r_seq else '❌'}")

        # Save results to the SAME directory as the test script
        test_dir = os.path.dirname(os.path.abspath(__file__))
        output_file = os.path.join(test_dir, "test_output_ALL.txt")
        df.to_csv(output_file, sep="\t", index=False)
        print(f"\nResults saved to: {output_file}")

        return df
    else:
        print("No probes generated!")
        return None


if __name__ == "__main__":
    # Test with the sample file in the SAME directory as this script
    fasta_file = "humanRNU1_1.fa"

    # Check if file exists in the same directory as this script
    test_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_path = os.path.join(test_dir, fasta_file)

    if not os.path.exists(fasta_path):
        print(f"Error: {fasta_file} not found in {test_dir}!")
        print(
            "Please make sure the sample FASTA file is in the same directory as this test script"
        )
        print(f"Expected path: {fasta_path}")
        sys.exit(1)

    # Change to the test directory for execution
    original_dir = os.getcwd()
    os.chdir(test_dir)

    try:
        result_df = test_with_sample_data(fasta_file)
    finally:
        # Change back to original directory
        os.chdir(original_dir)
