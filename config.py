# config.py
# Exact parameter values from R script
DEFAULT_SETTINGS = {
    "score_min": 0.9,
    "taille_sonde_max": 32,
    "taille_sonde_min": 26,
    "distance_min_inter_sonde": 2,
    "min_gc": 0.4,
    "max_gc": 0.6,
    "max_masked_percent": 0.1,
    "min_probe_per_transcript": 0,
    "pnas_filter_option": [1, 2, 4],
    "salt_conc": 0.115,
    # SIMPLIFIED: Use fixed dG37 value
    "fixed_dg37_value": -32.0,  # Always use -32 like R script behavior
}

# FLAP sequences - exact from R script
FLAP_SEQUENCES = {
    "X": "CCTCCTAAGTTTCGAGCTGGACTCAGTG",
    "Y": "TTACACTCGGACCTCGTCGACATGCATT",
    "Z": "CCAGCTTCTAGCATCCATGCCCTATAAG",
}

# Thermodynamic parameters
DG37_VALUES = {
    "AA": -0.2,
    "AC": -1.5,
    "AG": -0.9,
    "AT": -1.0,
    "CA": -1.0,
    "CC": -2.2,
    "CG": -1.2,
    "CT": -1.4,
    "GA": -0.8,
    "GC": -2.4,
    "GG": -1.5,
    "GT": -1.0,
    "TA": -0.3,
    "TC": -1.4,
    "TG": -1.0,
    "TT": -0.4,
}
