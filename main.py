# main.py - UPDATED with dustmasker functionality restored
import tkinter as tk
from tkinter import filedialog, messagebox
from rich.progress import track
import os
import pandas as pd
from pathlib import Path

from sequence_utils import read_fasta_sequences, create_output_directory
from oligostan_core import (
    optimize_dg37_selection,
    get_probes_from_rna_dg37,
    process_probes_for_output,
)
from filters import (
    is_ok_4_pnas_filter,
    is_ok_4_gc_filter,
    dustmasker_filter,
)  # FIXED: Added back dustmasker_filter
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

        # Use fixed dG37 value (simplified approach)
        optimal_dg37 = DEFAULT_SETTINGS["fixed_dg37_value"]

        # Generate probes with fixed dG37
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
                # UPDATED: Pass dustmasker parameters
                all_probes_data.extend(
                    process_probes_for_output(
                        probes, seq_data, optimal_dg37, **DEFAULT_SETTINGS
                    )
                )

        # Generate output files
        generate_output_files(all_probes_data, output_dir, base_name)

        return True

    except Exception as e:
        raise Exception(f"Error processing {file_path}: {str(e)}")


def generate_output_files(probes_data, output_dir, file_base_name):
    """Generate CSV files with exact column structure as R script"""
    if not probes_data:
        # Write empty file if no probes found
        with open(
            os.path.join(output_dir, f"Probes_{file_base_name}_FILT.txt"), "w"
        ) as f:
            f.write("No probes found after filtering. Change filtering parameters.\n")
        return

    # Create DataFrame
    df = pd.DataFrame(probes_data)

    # Sort by PNAS compliance (descending)
    df = df.sort_values("NbOfPNAS", ascending=False)

    # Save raw results (ALL)
    raw_filename = os.path.join(output_dir, f"Probes_{file_base_name}_ALL.txt")
    df.to_csv(raw_filename, sep="\t", index=False)

    # Filter for final results - UPDATED: Include dustmasker in filter logic
    use_dustmasker = DEFAULT_SETTINGS.get("use_dustmasker", False)

    if use_dustmasker:
        # Include dustmasker in filter criteria
        filtered_df = df[
            (df["GCFilter"] == 1)
            & (df["PNASFilter"] == 1)
            & (df["MaskedFilter"] == 1)  # Include dustmasker filter
        ]
    else:
        # Original filter criteria (dustmasker disabled)
        filtered_df = df[(df["GCFilter"] == 1) & (df["PNASFilter"] == 1)]

    # Save filtered results (FILT)
    filt_filename = os.path.join(output_dir, f"Probes_{file_base_name}_FILT.txt")
    filtered_df.to_csv(filt_filename, sep="\t", index=False)


def main():
    """Main application entry point"""
    print("Oligostan Python - smiFISH Probe Design Tool")
    print("=" * 50)

    # Show dustmasker status
    if DEFAULT_SETTINGS.get("use_dustmasker", False):
        print("🔍 dustmasker filter: ENABLED")
    else:
        print("🔍 dustmasker filter: DISABLED (default, matching R script)")

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
            print(f"✅ Successfully processed: {os.path.basename(file_path)}")
        except Exception as e:
            print(f"❌ Error processing {os.path.basename(file_path)}: {e}")
            continue

    print(f"\nBatch processing completed!")
    print(f"Successfully processed {success_count}/{len(files)} files")


if __name__ == "__main__":
    main()
