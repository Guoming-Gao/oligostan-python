# Oligostan_Python

A Python implementation of Oligostan for automated smiFISH (single-molecule fluorescence in situ hybridization) probe design with batch processing capabilities.

## Overview

Oligostan-Python is a complete port of the original R-based Oligostan tool, designed to generate high-quality oligonucleotide probes for RNA detection experiments. This implementation maintains exact compatibility with the original R script while adding modern batch processing features and a user-friendly GUI interface.

## Key Features

- **Exact R Script Compatibility** - Produces identical results to the original Oligostan R implementation
- **Batch Processing** - Process multiple FASTA files simultaneously with GUI file selection
- **Comprehensive Quality Control** - Multiple filtering layers including:
  - Thermodynamic optimization (dG37 calculations using nearest-neighbor parameters)
  - GC content filtering (40-60% default range)
  - PNAS composition rules (5 different filters for probe quality)
  - Optional repeat masking with dustmasker
- **Automated Probe Enhancement** - Adds FLAP sequences (X, Y, Z) for detection/amplification
- **Flexible Output** - Generates both filtered and unfiltered probe sets in CSV and FASTA formats
- **Cross-Platform** - Works on Windows, macOS, and Linux

## Installation

### Prerequisites

- Python 3.7+
- NCBI BLAST+ suite (for dustmasker functionality, optional)

### Install Dependencies

```
pip install -r requirements.txt
```

### Required Python Packages

```
biopython>=1.79
pandas>=1.3.0
numpy>=1.20.0
rich>=12.0.0
```

## Usage

### Quick Start

1. **Launch the application:**
   ```
   python main.py
   ```

2. **Select FASTA files** using the GUI dialog (supports multiple file selection)

3. **Results** are automatically saved in the same directory as your input files:
   ```
   input_directory/
   ├── your_sequence.fa
   └── Probes_your_sequence/
       ├── Probes_your_sequence_ALL.txt    # All generated probes
       ├── Probes_your_sequence_FILT.txt   # Filtered probes only
       └── *.fasta files                   # FASTA format outputs
   ```

### Testing

Run the included test with sample data:

```
cd oligostan_test/
python test_oligostan.py
```

Expected output: 5 probes matching the original R script results exactly.

## Project Structure

```
oligostan_python/
├── main.py                 # GUI entry point and batch processing
├── oligostan_core.py       # Core probe design algorithms
├── thermodynamics.py       # Delta G calculations (nearest-neighbor model)
├── filters.py              # Quality control filters (PNAS rules, GC content)
├── sequence_utils.py       # FASTA I/O and sequence operations
├── config.py              # Default parameters and settings
├── requirements.txt        # Python dependencies
├── README.md
├── LICENSE
└── oligostan_test/
    ├── test_oligostan.py   # Validation test script
    ├── humanRNU1_1.fa      # Sample input sequence
    └── *.txt               # Expected outputs for validation
```

## Configuration

Default parameters can be modified in `config.py`:

```
DEFAULT_SETTINGS = {
    'fixed_dg37_value': -32.0,      # Target dG37 for thermodynamic optimization
    'score_min': 0.9,               # Minimum probe score threshold
    'taille_sonde_max': 32,         # Maximum probe length (nucleotides)
    'taille_sonde_min': 26,         # Minimum probe length (nucleotides)
    'distance_min_inter_sonde': 2,  # Minimum spacing between probes
    'min_gc': 0.4,                  # Minimum GC content (40%)
    'max_gc': 0.6,                  # Maximum GC content (60%)
    'pnas_filter_option': [1][2][4] # PNAS composition rules to apply
}
```

### PNAS Filter Rules

1. **Rule 1**: Adenine content < 28%
2. **Rule 2**: No AAAA runs (prevents secondary structures)
3. **Rule 3**: Cytosine content between 22-28%
4. **Rule 4**: No CCCC runs (prevents secondary structures)
5. **Rule 5**: No 6-nucleotide windows with >50% cytosine

## Output Format

The tool generates comprehensive probe information including:

| Column | Description |
|--------|-------------|
| `dGOpt` | Optimized dG37 value used (-32.0) |
| `ProbesNames` | Sequence identifier from input file |
| `theStartPos/theEndPos` | Probe positions in original sequence |
| `ProbeSize` | Length of each probe (26-32 nt) |
| `Seq` | Probe sequence |
| `dGScore` | Thermodynamic score (0-1 scale) |
| `dG37` | Calculated delta G at 37°C (kcal/mol) |
| `GCpc` | GC percentage |
| `*Filter` | Pass/fail status for each quality filter |
| `NbOfPNAS` | Number of PNAS rules passed |
| `HybFlpX/Y/Z` | Probe sequences with FLAP sequences attached |

## Algorithm Details

The probe design algorithm follows these key steps:

1. **Sequence Processing** - Reverse complement input sequences to work from probe perspective
2. **Thermodynamic Analysis** - Calculate dG37 values using nearest-neighbor thermodynamics
3. **Multi-Length Optimization** - Test probe lengths from 26-32 nucleotides
4. **Scoring and Selection** - Score probes based on deviation from target dG37 (-32.0 kcal/mol)
5. **Quality Filtering** - Apply GC content and PNAS composition filters
6. **Spacing Optimization** - Ensure minimum distance between selected probes
7. **FLAP Addition** - Attach detection sequences for experimental amplification

### Thermodynamic Model

Uses nearest-neighbor thermodynamics with:
- 37°C temperature
- 0.115 M salt concentration
- Exact parameter set from original Oligostan R script

## Validation

This Python implementation has been extensively validated against the original R script:

- ✅ **Identical probe sequences** for all test cases
- ✅ **Matching thermodynamic calculations** (dG37 values)
- ✅ **Equivalent filter behavior** (PNAS rules, GC content)
- ✅ **Same output formatting** and column structure
- ✅ **Consistent probe selection** and ranking

## Performance

- **Processing Speed**: ~10-100x faster than R implementation
- **Memory Usage**: Optimized for typical transcript lengths (1-10 kb)
- **Batch Processing**: Handle multiple files efficiently with progress tracking

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Development Setup

```
git clone https://github.com/Guoming-Gao/oligostan-python.git
cd oligostan-python
pip install -r requirements.txt
python test_oligostan.py  # Run validation tests
```

## Citation

If you use Oligostan-Python in your research, please cite the original method

**Original Method:**
Tsanov, Nikolay, Aubin Samacoits, Racha Chouaib, et al. "smiFISH and FISH-Quant – a Flexible Single RNA Detection Approach with Super-Resolution Capability." Nucleic Acids Research 44, no. 22 (2016): e165–e165. https://doi.org/10.1093/NAR/GKW784

