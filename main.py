# main.py
import tkinter as tk
from tkinter import filedialog, messagebox
from rich.progress import track
import os
import pandas as pd
from pathlib import Path

from sequence_utils import read_fasta_sequences, create_output_directory
from oligostan_core import optimize_dg37_selection, get_probes_from_rna_dg37
from filters import is_ok_4_pnas_filter, is_ok_4_gc_filter, dustmasker_filter
from config import DEFAULT_SETTINGS, FLAP_SEQUENCES


def select_fasta_files():
    """GUI file selection"""
    root = tk.Tk()
    root.withdraw()

    files = filedialog.askopenfilenames(
        title="Select FASTA files for probe design",
        filetypes=[("FASTA files", "*.fa *.fasta *.fas"), ("All files", "*.*")],
    )

    root.destroy()
    return files


def process_single_file(file_path):
    """Process a single FASTA file"""
    try:
        # Read sequences
        sequences = read_fasta_sequences(file_path)

        # Create output directory
        output_dir = create_output_directory(file_path)
        base_name = os.path.splitext(os.path.basename(file_path))[0]

        # Generate dG37 range
        dg37_range = [
            DEFAULT_SETTINGS["the_dg37_min"] + i * DEFAULT_SETTINGS["the_dg37_step"]
            for i in range(
                int(
                    (
                        DEFAULT_SETTINGS["the_dg37_max"]
                        - DEFAULT_SETTINGS["the_dg37_min"]
                    )
                    / DEFAULT_SETTINGS["the_dg37_step"]
                )
                + 1
            )
        ]

        # Optimize dG37 selection
        optimal_dg37 = optimize_dg37_selection(
            sequences, dg37_range, **DEFAULT_SETTINGS
        )

        # Generate probes with optimal dG37
        all_probes_data = []
        for seq_data in sequences:
            probes = get_probes_from_rna_dg37(
                seq_data["sequence"],
                min_size_probe=DEFAULT_SETTINGS["taille_sonde_min"],
                max_size_probe=DEFAULT_SETTINGS["taille_sonde_max"],
                desired_dg=optimal_dg37,
                min_score_value=DEFAULT_SETTINGS["score_min"],
                inc_betw_prob=DEFAULT_SETTINGS["distance_min_inter_sonde"],
            )

            if probes:
                all_probes_data.extend(
                    process_probes_for_output(probes, seq_data, optimal_dg37)
                )

        # Generate output files
        generate_output_files(all_probes_data, output_dir, base_name)

        return True

    except Exception as e:
        raise Exception(f"Error processing {file_path}: {str(e)}")


def process_probes_for_output(probes, seq_data, dg37_value):
    """Process probes and add all required information"""
    processed_probes = []

    for probe in probes:
        probe_size, score, position, sequence = probe

        # Calculate all required metrics
        probe_info = {
            "dGOpt": dg37_value,
            "ProbesNames": seq_data["id"],
            "theStartPos": len(seq_data["sequence"]) - position + 1,
            "theEndPos": len(seq_data["sequence"]) - position + 1 + probe_size,
            "ProbeSize": probe_size,
            "Seq": sequence,
            "dGScore": score,
            "dG37": dg37_value,  # This would be recalculated in full implementation
            "GCpc": (sequence.count("G") + sequence.count("C")) / len(sequence) * 100,
            "GCFilter": is_ok_4_gc_filter(
                sequence, DEFAULT_SETTINGS["min_gc"], DEFAULT_SETTINGS["max_gc"]
            ),
            "aCompFilter": is_ok_4_pnas_filter(sequence, [1]),
            "aStackFilter": is_ok_4_pnas_filter(sequence, [2]),
            "cCompFilter": is_ok_4_pnas_filter(sequence, [3]),
            "cStackFilter": is_ok_4_pnas_filter(sequence, [4]),
            "cSpecStackFilter": is_ok_4_pnas_filter(sequence, [5]),
            "PNASFilter": is_ok_4_pnas_filter(
                sequence, DEFAULT_SETTINGS["pnas_filter_option"]
            ),
            "HybFlpX": sequence + FLAP_SEQUENCES["X"],
            "HybFlpY": sequence + FLAP_SEQUENCES["Y"],
            "HybFlpZ": sequence + FLAP_SEQUENCES["Z"],
        }

        # Calculate PNAS sum
        probe_info["NbOfPNAS"] = sum(
            [
                probe_info["aCompFilter"],
                probe_info["aStackFilter"],
                probe_info["cCompFilter"],
                probe_info["cStackFilter"],
                probe_info["cSpecStackFilter"],
            ]
        )

        processed_probes.append(probe_info)

    return processed_probes


def generate_output_files(probes_data, output_dir, base_name):
    """Generate output CSV files with exact R script structure"""
    if not probes_data:
        # Write empty file if no probes found
        with open(os.path.join(output_dir, f"Probes_{base_name}_FILT.txt"), "w") as f:
            f.write("No probes found after filtering. Change filtering parameters.\n")
        return

    # Create DataFrame
    df = pd.DataFrame(probes_data)

    # Sort by PNAS compliance (descending)
    df = df.sort_values("NbOfPNAS", ascending=False)

    # Save raw results (ALL)
    raw_filename = os.path.join(output_dir, f"Probes_{base_name}_ALL.txt")
    df.to_csv(raw_filename, sep="\t", index=False)

    # Filter for final results
    filtered_df = df[(df["GCFilter"] == True) & (df["PNASFilter"] == True)]

    # Save filtered results (FILT)
    filt_filename = os.path.join(output_dir, f"Probes_{base_name}_FILT.txt")
    filtered_df.to_csv(filt_filename, sep="\t", index=False)


def main():
    """Main application entry point"""
    print("Oligostan Python - smiFISH Probe Design Tool")
    print("=" * 50)

    # Select FASTA files
    files = select_fasta_files()

    if not files:
        print("No files selected. Exiting.")
        return

    print(f"Selected {len(files)} files for processing")

    # Process files with progress tracking
    success_count = 0
    for file_path in track(files, description="Processing files..."):
        try:
            process_single_file(file_path)
            success_count += 1
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue

    print(f"\nBatch processing completed!")
    print(f"Successfully processed {success_count}/{len(files)} files")


if __name__ == "__main__":
    main()
