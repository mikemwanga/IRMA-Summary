# IRMA-Summary

Workflow to summarise IRMA assembly output into tabular summary files and publication-ready plots.

## What This Tool Does

`IRMA-Summary` scans a parent directory that contains multiple IRMA sample folders (for example `barcode01`, `barcode02`, ...), extracts assembly/read-depth information, calculates summary metrics, then writes:

- TSV tables to: `IRMA_Summary_<timestamp>/summary_data/`
- PDF plots to: `IRMA_Summary_<timestamp>/figures/`
- Per-sample depth PDFs to: `IRMA_Summary_<timestamp>/figures/individual_sample_depth/`

## Expected Input Structure

Run the script on a directory that contains IRMA sample folders. Each sample folder should include:

- `tables/*coverage.txt`
- `tables/READ_COUNTS.txt`
- `intermediate/0-ITERATIVE-REFERENCES/R0-*.ref`
- optional subtype FASTA files (`*_HA_*.fasta`, `*_NA_*.fasta`)

## How To Run

```bash
python3 main.py <PATH_TO_PARENT_DIR>
```

Example:

```bash
python3 main.py dataset
```

## Core Calculations

### 1. Average read depth per gene per sample
Input:
- `tables/*coverage.txt`

Calculation:
- `coverage = mean(Consensus_Count)` for each gene segment in each sample

Used in:
- `depth_coverage_data.tsv`
- `depth_coverage_data_log_10.tsv`
- heatmap figures

### 2. Absolute assembled length per gene
Input:
- per-sample coverage tables (`*coverage.txt`)

Calculation:
- `seq_len = number of rows in each gene coverage table` (assembled positions)

Used in:
- `gene_absolute_length_coverage_data.tsv`

### 3. Percentage gene coverage relative to IRMA reference
Input:
- assembled lengths (above)
- `intermediate/0-ITERATIVE-REFERENCES/R0-*.ref`

Calculation:
- `gene_percent_coverage = (assembled_length / reference_length) * 100`

Used in:
- `gene-percentage-length-coverage-data.tsv`
- `gene-proportion-length-coverage-data.tsv`
- gene-length heatmap figures

### 4. Assembled reads percentage
Input:
- `tables/READ_COUNTS.txt` (`initial`, `match` rows)

Calculations:
- `total_reads = initial`
- `total_assembled = match`
- `prop_assembled = (match / initial) * 100`
- `notassembled = initial - match`
- `prop_notassembled = 100 - prop_assembled`

Used in:
- `assembled-reads-data.tsv`
- `read_distribution-Read Proportion(%).pdf`
- `read_distribution-Read Frequency.pdf`
- `assembled_reads_Subset-*.pdf`

IRMA provenance reference:
- Values are sourced from IRMA `READ_COUNTS.txt` records:
  - `1-initial` -> `total_reads`
  - `3-match` -> `total_assembled`
- `prop_assembled` is computed as `(3-match / 1-initial) * 100`.

### 5. Log-transformed depth for line plots
Input:
- merged per-position depth tables (`<sample>-read-coverage-table.tsv`)

Calculation:
- `depth_log = log10(Consensus_Count + 1)`

Used in:
- `coverage-lines-log-counts.pdf`

### 6. Subtype inference
Input:
- sample FASTA filenames (HA/NA patterns)

Logic:
- collect HA/NA subtype tokens from filenames
- fallback to `Hx`/`Nx` when missing

Used in:
- `subtypes.tsv`
- `influenza_subtypes.tsv`
- sample naming in summary matrices

### 7. Gene count table
Input:
- `tables/READ_COUNTS.txt`

Logic:
- remove non-gene records (`initial`, `failQC`, `passQC`, `chimeric`, `nomatch`, `match`, `altmatch`)
- keep per-sample gene read counts

Used in:
- `gene_count.tsv`

## Output Files: Summary Tables (`summary_data/`)

### Global summary tables

- `assembled-reads-data.tsv`
  - assembled/non-assembled read metrics per sample
- `depth_coverage_data.tsv`
  - average read depth per gene per sample
- `depth_coverage_data_log_10.tsv`
  - log10-transformed depth matrix
- `gene_absolute_length_coverage_data.tsv`
  - absolute assembled gene lengths
- `gene-percentage-length-coverage-data.tsv`
  - % coverage by gene relative to reference
- `gene-proportion-length-coverage-data.tsv`
  - same proportional coverage matrix (kept for compatibility)
- `gene_count.tsv`
  - per-gene read counts per sample
- `subtypes.tsv`
  - sample ID and inferred subtype
- `influenza_subtypes.tsv`
  - subtype table in compatibility format

### Intermediate per-sample position tables (for plotting)

- `<sample>-read-coverage-table.tsv`
  - position-level depth merged with reference coordinates
  - used to generate line plots

## Output Files: Figures (`figures/`)

### Heatmaps

- `coverage-depth-log-counts-Subset_*.pdf`
  - heatmap of log10 depth values
- `coverage_length_Subset-*.pdf`
  - heatmap of percentage gene coverage

### Read distribution

- `read_distribution-Read Proportion(%).pdf`
  - stacked bars of assembled% vs not assembled%
- `read_distribution-Read Frequency.pdf`
  - stacked bars of assembled vs non-assembled read counts
- `assembled_reads_Subset-*.pdf`
  - bars of assembled percentage with total read annotations

### Line plots

- `coverage-lines-log-counts.pdf`
  - faceted multi-sample line plot (log depth) by gene

### Individual sample depth curves

- `individual_sample_depth/<sample>_depth.pdf`
  - concatenated gene-position coverage curve for each sample

## Notes

- The script creates a new timestamped output directory each run.
- Existing temporary sample `tmp/` folders are recreated per run.
- If static image export fails in restricted environments (Kaleido/Chromium sandbox issue), TSV tables may still be produced successfully.
