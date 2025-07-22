# sequence_utils.py
from Bio import SeqIO
from Bio.Seq import Seq
import os


def read_fasta_sequences(file_path):
    """Read FASTA sequences and return as list with reverse complement"""
    sequences = []

    # Extract base filename without extension (e.g., "humanRNU1_1" from "humanRNU1_1.fa")
    base_filename = os.path.splitext(os.path.basename(file_path))[0]

    for record in SeqIO.parse(file_path, "fasta"):
        # Use base filename instead of sequence header (FIXED!)
        seq_id = base_filename

        # Reverse complement to work from probe perspective (matching R script)
        rev_comp_seq = str(record.seq.reverse_complement())

        sequences.append(
            {
                "id": seq_id,
                "name": base_filename,  # Use filename base
                "sequence": rev_comp_seq,
            }
        )

    return sequences


def create_output_directory(input_file_path):
    """Create output directory in same location as input file"""
    input_dir = os.path.dirname(input_file_path)
    base_name = os.path.splitext(os.path.basename(input_file_path))[0]
    output_dir = os.path.join(input_dir, f"Probes_{base_name}")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return output_dir
