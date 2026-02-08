# `rt_pct_strain_source_2time_point.R` Detailed README

This README documents the full behavior of:

- `Jean/pipeline_v2/rt_pct_strain_source_2time_point.R`

It is intended for day-to-day use, debugging, and onboarding new users.

---

## 1. What This Script Does

This script runs a qPCR analysis pipeline for **two time points only**.

At a high level, it:

1. Reads one or more qPCR CSV files.
2. Cleans and standardizes sample/gene fields.
3. Extracts time labels from sample names.
4. Validates that data is compatible with 2-time-point comparison.
5. Builds per-replicate Ct tables (`control_gene` + target gene).
6. Computes delta and fold-change summaries.
7. Exports Excel workbooks.
8. Produces figure panels (single and combined).

---

## 2. Input Requirements

### 2.1 File Format

Input files are expected to be qPCR CSV exports that contain, at minimum, columns equivalent to:

- `Well`
- `Sample name`
- `Sample type`
- `Dye`
- `Gene`
- `Ct`

The parser is robust to vendor metadata rows above the table and header irregularities (e.g. extra commas, hidden chars, blank columns).

### 2.2 Time Label Requirement (Strict)

This script is for **2-time-point analysis**. It requires explicit time tokens in sample names, such as:

- `b73_0h`, `b73_1h`
- `opr21h-2`
- `S_1H`
- `1hr`

If no time token is detected in a file, the script stops with an error.

### 2.3 Number of Time Points (Strict)

After parsing:

- If no time points are detected -> error.
- If more than 2 unique time points are detected -> error.

---

## 3. CLI Usage

```bash
Rscript Jean/pipeline_v2/rt_pct_strain_source_2time_point.R [options]
```

### 3.1 Core Options

- `--raw-dir PATH`
  Input directory for raw CSV files.

- `--files CSV`
  Comma-separated file list (absolute or relative). Overrides `--pattern`.

- `--pattern REGEX`
  Regex for file selection in `--raw-dir` (used only if `--files` is not provided).

- `--parent PATH`
  Parent output directory.

- `--out-name NAME`
  Output folder name under `--parent`.

- `--out-dir PATH`
  Full output path. Overrides `--parent` + `--out-name`.

### 3.2 Analysis Control Options

- `--control-gene GENE`
  Control gene for delta calculation.
  If omitted, auto-detect order is:
  1. `5.8S`
  2. `TUB`

- `--control-time T`
  Control time point used in 2-time comparisons.
  Default: `0h`.

- `--sample-name-mode M`
  Time extraction mode:
  - `auto` (default): flexible extraction (embedded formats allowed)
  - `strain_time`: strict suffix extraction from `strain_time` format

- `--exclude-samples CSV`
  Comma-separated sample names to remove before analysis.
  Match is done on normalized sample names (lowercase, spaces removed).

- `--max-sample N`
  Maximum retained replicate count after outlier dropping.
  Default: `3`.

### 3.3 Calibration Options

- `--perform-calibration true|false`
  Enable or disable cross-file calibration.

- `--calibration-sample S`
  Sample used for calibration.

- `--calibration-gene G`
  Gene used for calibration (default = control gene).

- `--calibration-time T`
  Time used for calibration (default = `0h`).

If calibration is enabled and sample is not provided, script tries to auto-find a common `*standard*` sample across files.

### 3.4 Plot Ordering Option

- `--sample-order CSV`
  Optional preferred order for combined plots.
  Strains not listed are still included.

### 3.5 Help

- `--help`
  Prints option summary and exits.

---

## 4. Sample Name Parsing Rules

The script normalizes sample names and extracts time labels.

### 4.1 Normalization

- Lowercase conversion.
- Removes trailing replicate suffixes such as `-2` / `_2`.
- Removes spaces.
- Converts `opr67` -> `opr78`.

### 4.2 Time Token Examples Recognized in `auto`

- `0h`, `1h`, `24h`
- `1hr`
- `1h-2`
- Embedded like `b731h-2`
- Underscore style like `lox10_1H`

### 4.3 `strain_time` Mode

In `strain_time`, time is extracted only from suffix pattern like:

- `strain_1h`
- `strain_1h-2`

If your names are not in this structure, use `auto`.

---

## 5. Validation Rules and Failure Conditions

The script intentionally fails fast in these cases:

1. File not found.
2. Required qPCR columns cannot be mapped.
3. No time token in any sample name within a file.
4. More than two unique time points after parsing.
5. Requested `--control-time` not present in detected time points.
6. `--control-gene` specified but not found in data.
7. Calibration requested but reference values missing.

---

## 6. Outputs

Output directory:

- `--out-dir`
- or `<parent>/<out-name>`

### 6.1 Excel Files

- `all_raw_results.xlsx`
- `all_raw_time_results.xlsx`
- `all_time_results.xlsx` (if sheets exist)

### 6.2 Figures

Saved under:

- `<out_dir>/figure`

Generated plot types include:

- Per-result bar plots (`*_result.png`)
- Combined-by-miR plots (`miR_*_combined.png`)
- Combined-by-strain plots (`strain_*_combined.png`)

Combined figures are saved with white background and square padded layout.

---

## 7. Common Messages (Not Errors)

When R loads packages, you may see masking messages like:

- `between, first, last`
- `filter, lag`

These are normal namespace messages, not pipeline failures.

---

## 8. Troubleshooting

### 8.1 "No time point token found"

- Check if sample names truly include time markers.
- Try `--sample-name-mode auto` first.
- If using strict `strain_time`, ensure format is `strain_0h`, `strain_1h`, etc.

### 8.2 "supports exactly two time points"

- This script is only for two-time-point experiments.
- Remove extra time points or use another pipeline for multi-time analysis.

### 8.3 "Control time not present"

- Set `--control-time` to one of the detected time labels.

### 8.4 Missing plots

- If all rows were filtered (e.g. Ct threshold / exclusions), plot lists can be empty.
- Check Excel outputs first to confirm upstream data availability.

### 8.5 Calibration issues

- Provide explicit `--calibration-sample`, `--calibration-gene`, `--calibration-time`.
- Confirm those combinations actually exist in each file.

---

## 9. Recommended Command Templates

### 9.1 Single-file, auto mode

```bash
Rscript Jean/pipeline_v2/rt_pct_strain_source_2time_point.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/your_file.csv" \
  --perform-calibration false \
  --out-name run_single
```

### 9.2 Two-file, explicit control gene and exclusions

```bash
Rscript Jean/pipeline_v2/rt_pct_strain_source_2time_point.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/f1.csv,/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/f2.csv" \
  --control-gene 5.8S \
  --control-time 0h \
  --exclude-samples "opr78,lox10" \
  --perform-calibration false \
  --out-name run_pair
```

### 9.3 Strict `strain_time` mode

```bash
Rscript Jean/pipeline_v2/rt_pct_strain_source_2time_point.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/f.csv" \
  --sample-name-mode strain_time \
  --control-time 0h \
  --out-name run_strict
```

### 9.4 With calibration

```bash
Rscript Jean/pipeline_v2/rt_pct_strain_source_2time_point.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/f1.csv,/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/f2.csv" \
  --perform-calibration true \
  --calibration-sample b73standard \
  --calibration-gene 5.8S \
  --calibration-time 0h \
  --out-name run_calibrated
```

---

## 10. Notes for Future Users

- Keep sample naming consistent within a project.
- Prefer explicit and stable time labels (`0h`, `1h`, `24h`).
- Use `--exclude-samples` for known problematic strains instead of manually editing source files.
- For >2 time points, do not force this script; use a multi-time-point pipeline instead.

