# oligostan_core.py - UPDATED with dustmasker integration
import pandas as pd
import numpy as np
from thermodynamics import dg_calc_rna_37, dg37_score_calc
from filters import (
    is_ok_4_pnas_filter,
    is_ok_4_gc_filter,
    dustmasker_filter,
    is_it_ok_4_a_comp,
    is_it_ok_4_a_stack,
    is_it_ok_4_c_comp,
    is_it_ok_4_c_stack,
    is_it_ok_4_c_spec_stack,
)
from config import DEFAULT_SETTINGS, FLAP_SEQUENCES


def which_max_r(x):
    """Exact translation of R's WhichMax function"""
    max_val = np.max(x)
    max_indices = np.where(x == max_val)[0]

    if len(max_indices) >= 2:
        return [0, max_val]  # Return 0 for ties (R behavior)
    else:
        return [max_indices[0] + 1, max_val]  # Return 1-based index


def get_probes_from_rna_dg37(
    seq,
    min_size_probe=26,
    max_size_probe=32,
    desired_dg=-32,
    min_score_value=0.9,
    inc_betw_prob=2,
):
    """Exact translation of getProbesFromRNAdG37 from R"""
    if isinstance(seq, list):
        seq = "".join(seq).upper()

    diff_size = max_size_probe - min_size_probe

    # R: dGCalc.RNA.37(Seq, ProbeLength = MaxSizeProbe) -> TheTmsTmp
    the_tms_tmp = dg_calc_rna_37(seq, probe_length=max_size_probe)
    nb_of_probes = len(the_tms_tmp)

    # Build matrix like R: start with max size, then add smaller sizes
    if diff_size > 0:
        # R: for (i in seq(DiffSize - 1, 0, -1))
        # R builds columns from right to left: max_size, max_size-1, ..., min_size
        all_columns = [the_tms_tmp]  # Start with max size

        for i in range(diff_size - 1, -1, -1):
            probe_length = min_size_probe + i
            dg_values = dg_calc_rna_37(seq, probe_length=probe_length)
            # Truncate to match shortest length
            min_len = min(len(dg_values), nb_of_probes)
            all_columns.insert(0, dg_values[:min_len])  # Insert at beginning

        # Build matrix: rows = positions, columns = probe sizes (min to max)
        max_len = max(len(col) for col in all_columns)
        the_tms_matrix = np.full((max_len, diff_size + 1), np.nan)

        for col_idx, column in enumerate(all_columns):
            the_tms_matrix[: len(column), col_idx] = column

    else:
        # Only one size
        the_tms_matrix = np.array(the_tms_tmp).reshape(-1, 1)

    # R: dG37ScoreCalc(TheTmsTmp, Desireddg) -> TmScores
    tm_scores = np.zeros_like(the_tms_matrix)
    for col in range(the_tms_matrix.shape[1]):
        valid_mask = ~np.isnan(the_tms_matrix[:, col])
        if np.any(valid_mask):
            tm_scores[valid_mask, col] = dg37_score_calc(
                the_tms_matrix[valid_mask, col], desired_dg
            )

    # R: t(apply(TmScores, 1, WhichMax)) -> BestScores
    best_scores = []
    for row in range(tm_scores.shape[0]):
        row_data = tm_scores[row, :]
        if not np.all(np.isnan(row_data)) and not np.all(row_data == 0):
            valid_data = row_data[~np.isnan(row_data)]
            if len(valid_data) > 0:
                result = which_max_r(valid_data)
                best_scores.append(result)
            else:
                best_scores.append([0, 0])
        else:
            best_scores.append([0, 0])

    best_scores = np.array(best_scores)

    # R: BestScores[, 1] + (MinSizeProbe - 1) -> BestScores[, 1]
    # Convert column index to actual probe size
    for i in range(len(best_scores)):
        if best_scores[i, 0] != 0:  # Not a tie
            # R adds (MinSizeProbe - 1) to the column index
            best_scores[i, 0] = best_scores[i, 0] + (min_size_probe - 1)
        else:  # Tie case - use minimum size
            best_scores[i, 0] = min_size_probe

    # R: cbind(BestScores, seq(1:length(BestScores[, 1]))) -> BestScores
    positions = np.arange(1, len(best_scores) + 1).reshape(-1, 1)
    best_scores_with_pos = np.column_stack([best_scores, positions])

    # R: BestScores[BestScores[, 2] >= MinScoreValue, ] -> ValidedScores
    validated_scores = best_scores_with_pos[
        best_scores_with_pos[:, 1] >= min_score_value
    ]

    if len(validated_scores) == 0:
        return None

    # Apply spacing constraint exactly like R
    the_probes = []

    # R: ValidedScores[order(ValidedScores[, 3]), ] -> ValidedScores
    if len(validated_scores.shape) == 1:
        validated_scores = validated_scores.reshape(1, -1)

    validated_scores = validated_scores[np.argsort(validated_scores[:, 2])]

    pointeur = 0  # R starts with 0
    while pointeur < len(seq):
        # R: ValidedScores[ValidedScores[, 3] >= Pointeur, ]
        valid_tmp = validated_scores[validated_scores[:, 2] >= pointeur]

        if len(valid_tmp) > 0:
            if len(valid_tmp.shape) == 1:
                valid_tmp = valid_tmp.reshape(1, -1)

            # Take first valid probe
            probe_size = int(valid_tmp[0, 0])
            score = valid_tmp[0, 1]
            position = int(valid_tmp[0, 2])

            # R: substr(Seq, start = ValiTmp[1, 3], stop = (ValiTmp[1, 3] + ValiTmp[1, 1] - 1))
            probe_seq = seq[position - 1 : position - 1 + probe_size].upper()

            the_probes.append([probe_size, score, position, probe_seq])

            # R: Pointeur <- (ValiTmp[1, 3] + ValiTmp[1, 1] + IncBetwProb)
            pointeur = position + probe_size + inc_betw_prob
        else:
            break

    return the_probes


def optimize_dg37_selection(sequences, dg37_range=None, **params):
    """SIMPLIFIED: Always return fixed dG37 value (no optimization)"""
    return params.get("fixed_dg37_value", -32.0)


def process_probes_for_output(probes, seq_data, dg37_value, **params):
    """Process probes exactly like R script - UPDATED with optional dustmasker"""
    processed_probes = []

    # RESTORED: Apply dustmasker filter if enabled
    use_dustmasker = params.get("use_dustmasker", False)
    max_masked_percent = params.get("max_masked_percent", 0.1)

    if use_dustmasker and probes:
        # Extract sequences for dustmasker
        probe_sequences = [probe[3] for probe in probes]
        dustmasker_results, masked_percentages = dustmasker_filter(
            probe_sequences, max_masked_percent
        )
    else:
        # Default: pass all probes (MaskedFilter <- FALSE behavior)
        dustmasker_results = [True] * len(probes) if probes else []
        masked_percentages = [0.0] * len(probes) if probes else []

    for i, probe in enumerate(probes):
        probe_size, score, position, sequence = probe

        # Ensure sequence is uppercase
        sequence = sequence.upper()

        # R position calculation:
        # (seqlength - ProbeList[[probeListNb]][i, 3] + 1) -> EndPosTmp
        # (EndPosTmp - ProbeList[[probeListNb]][i, 1]) -> StartPosTmp
        seq_length = len(seq_data["sequence"])
        the_end_pos = seq_length - position + 1
        the_start_pos = the_end_pos - probe_size  # FIXED: Removed +1 to match R exactly

        # Recalculate actual dG37 for this probe
        actual_dg37 = dg_calc_rna_37(sequence, probe_length=len(sequence))[0]

        # Calculate GC percentage
        gc_count = sequence.count("G") + sequence.count("C")
        gc_percentage = gc_count / len(sequence)

        # All filters
        gc_filter_pass = (
            1
            if is_ok_4_gc_filter(
                sequence, DEFAULT_SETTINGS["min_gc"], DEFAULT_SETTINGS["max_gc"]
            )
            else 0
        )
        a_comp_pass = 1 if is_it_ok_4_a_comp(sequence) else 0
        a_stack_pass = 1 if is_it_ok_4_a_stack(sequence) else 0
        c_comp_pass = 1 if is_it_ok_4_c_comp(sequence) else 0
        c_stack_pass = 1 if is_it_ok_4_c_stack(sequence) else 0
        c_spec_pass = 1 if is_it_ok_4_c_spec_stack(sequence) else 0
        pnas_filter_pass = (
            1
            if is_ok_4_pnas_filter(sequence, DEFAULT_SETTINGS["pnas_filter_option"])
            else 0
        )

        # RESTORED: dustmasker filter results
        dustmasker_pass = (
            1 if (i < len(dustmasker_results) and dustmasker_results[i]) else 0
        )
        repeat_masker_pc = masked_percentages[i] if i < len(masked_percentages) else 0.0

        # PNAS sum
        nb_of_pnas = (
            a_comp_pass + a_stack_pass + c_comp_pass + c_stack_pass + c_spec_pass
        )

        # Format exactly like R output
        probe_info = {
            "dGOpt": dg37_value,  # NUMBER, not string
            "ProbesNames": seq_data["name"],  # Now uses filename base
            "theStartPos": the_start_pos,  # FIXED: Now matches R positions
            "theEndPos": the_end_pos,
            "ProbeSize": probe_size,
            "Seq": sequence,
            "dGScore": score,
            "dG37": actual_dg37,
            "GCpc": gc_percentage,
            "GCFilter": gc_filter_pass,
            "aCompFilter": a_comp_pass,
            "aStackFilter": a_stack_pass,
            "cCompFilter": c_comp_pass,
            "cStackFilter": c_stack_pass,
            "cSpecStackFilter": c_spec_pass,
            "NbOfPNAS": nb_of_pnas,
            "PNASFilter": pnas_filter_pass,
            "MaskedFilter": dustmasker_pass,  # RESTORED: dustmasker filter result
            "RepeatMaskerPC": repeat_masker_pc,  # RESTORED: masked percentage
            "InsideUTR": 0,
            "HybFlpX": sequence + FLAP_SEQUENCES["X"],
            "HybFlpY": sequence + FLAP_SEQUENCES["Y"],
            "HybFlpZ": sequence + FLAP_SEQUENCES["Z"],
        }

        processed_probes.append(probe_info)

    return processed_probes
