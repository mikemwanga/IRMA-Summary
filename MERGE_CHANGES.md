# IRMA-Summary Merge Changes

This merge combines stable behavior from root scripts with selected upgrades from `irma_dev/scripts`.

## Files Updated

- `config.py`
- `figures.py`
- `summarise_IRMA.py`

## Where Changes Happened

1. `config.py`
- Added `plotly.io as pio` and enabled `pio.kaleido.scope.mathjax = None`.
- Added shared visualization constants (`textsize`, `bgcolor`, `grid_color`, etc.).
- Added `read_colormap` for new assembled vs non-assembled distribution plots.

2. `figures.py`
- Added `figure_aesthetics()` for consistent styling.
- Added `plot_read_distribution()` for proportion and absolute read distribution charts.
- Added `plot_all_line_plots()` for faceted all-sample log coverage plot.
- Added `export_data(data, name, index=True)` to avoid accidental index/schema regressions.
- Fixed output filenames to consistent semantics:
  - log heatmaps: `coverage-depth-log-counts-*.pdf`
  - length heatmaps: `coverage_length_*.pdf`
- Kept legacy `plot_assembled_reads()` output (`assembled_reads_*.pdf`) and made it tolerant of either `assembled` or `prop_assembled` input column.

3. `summarise_IRMA.py`
- Merged robustness checks from `irma_dev`:
  - skip missing per-gene coverage/reference temp files.
  - subtype fallback to `Hx/Nx` and global fallback to `HxNx`.
- Added safe empty-input handling for coverage-table concat (`process_data()`).
- Fixed `irma_dev` regression where percentage filename could receive absolute lengths:
  - `gene-percentage-length-coverage-data.tsv` now exports proportional percentage matrix.
- Added compatibility outputs:
  - `subtypes.tsv` and `influenza_subtypes.tsv`.
- Added new plots:
  - `read_distribution-Read Proportion(%).pdf`
  - `read_distribution-Read Frequency.pdf`
  - `coverage-lines-log-counts.pdf`
- Preserved legacy outputs and workflow order where possible.

## Notes

- `irma_dev` directory itself was not rewritten in this merge. Clean merged pipeline lives in root scripts.
- Any existing untracked/generated files in repository were left untouched.
