# Pan-Methylomes Framework

This script (`pan_methylomes_framework.py`) processes bacterial DNA methylation data to build a comprehensive pan-methylome profile. It identifies target motifs across multiple strains, calculates methylation frequencies for both sense and antisense strands, and outputs summary tables for downstream comparative analyses.

## Features

- Search for user-specified DNA motifs (default: `GATC`) in aligned gene sequences.
- Extract methylation information from Tombo output.
- Calculate methylation frequency, coverage, and overlap across strains.
- Output separate results for sense and antisense strands.

## Requirements

- Python 3.7+
- Dependencies:
  - `pandas`
  - `numpy`
  - `scipy`
  - `biopython`

Install required Python packages:
```bash
pip install pandas numpy scipy biopython
```

Additionally, `mafft` is recommended for sequence alignment.

### Required Arguments

| Argument         | Description |
|------------------|-------------|
| `--index`        | Path to strain sequence index file |
| `--gene_index`   | Path to gene index file |
| `--methy`        | Path to Tombo methylation results |
| `--seq_dir`      | Path to sequence files |
| `--motif`        | DNA motif to search (default: `GATC`) |
| `--gene_info`    | Path to gene information files |
| `--out`          | Output directory |
| `--gene_suffix`  | Sequence file suffix (default: `fa`) |


## Input Format

### `--index` (Strain Index CSV)
CSV file containing strain IDs to be analyzed.

### `--gene_index` (Gene Index CSV)
CSV file listing gene IDs to be analyzed.

### `--methy` (Methylation Results Directory)
Folder containing Tombo methylation output files (`combine_freq_cover.tsv`) for each gene. The information in each column is chromosome number, methylation start position, methylation end position, number of sites detected, methylation ratio, and methylation direction.

### `--gene_info` (Gene Information Directory)
Folder containing per-strain gene information CSV files with genomic coordinates and overlap data. The information in each column is chromosome number, gene start position, gene end position, gene direction, gene overlap, and gene name.


## Output
- `pan_methy_sense.tab` – Methylation statistics for the sense strand
- `pan_methy_antisense.tab` – Methylation statistics for the antisense strand


### Output format
| Gene_ID | Site | Motif count |Methylation count| Strain count|Overlapping info|Length|Methylated_Strain_List |
|---------|------|-----|---|--|-|---|----------------------|
| 00001   | 125  | 80|60|84|0 |1000| strainA,strainC,...     |

