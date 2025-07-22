# thermodynamics.py
import pandas as pd
import numpy as np
from config import DG37_VALUES


def convert_rna_seq_2_delta_g_at_37(rna_seq):
    """Exact translation of ConvertRNASeq2DeltaGat37"""
    rna_seq = rna_seq.upper()
    nb_base = len(rna_seq)

    dg = []
    dim = []

    for i in range(nb_base - 1):
        dimer = rna_seq[i : i + 2]
        dim.append(dimer)
        dg.append(DG37_VALUES.get(dimer, 0))

    return pd.DataFrame({"dim": dim, "dG": dg})


def dg_calc_rna_37(rna_seq, probe_length=31, salt_conc=0.115):
    """Exact translation of dGCalc.RNA.37"""
    rna_seq_conv = convert_rna_seq_2_delta_g_at_37(rna_seq)

    # Rolling sum for probe_length-1 consecutive dimers
    rolling_sum = []
    dg_values = rna_seq_conv["dG"].values

    for i in range(len(dg_values) - probe_length + 2):
        if i + probe_length - 1 <= len(dg_values):
            window_sum = np.sum(dg_values[i : i + probe_length - 1])
            rolling_sum.append(window_sum)

    # Apply salt correction: dG - ((log(SaltConc) * -0.175) - 0.2)
    all_dg = [dg - ((np.log(salt_conc) * -0.175) - 0.2) for dg in rolling_sum]

    return all_dg


def dg37_score_calc(the_dg37, desired_dg=-33):
    """Exact translation of dG37ScoreCalc"""
    if isinstance(the_dg37, (list, np.ndarray)):
        return [(-0.1 * abs(dg - desired_dg)) + 1 for dg in the_dg37]
    else:
        return (-0.1 * abs(the_dg37 - desired_dg)) + 1
