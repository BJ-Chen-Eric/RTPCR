# RT-qPCR Unified Pipeline (v2)

This README covers both:

- `work_pp/rt_pct_strain_source_unified.R`
- `pipeline_Rpack` (R package version)

The script supports both:

- 2-time-point analysis
- multi-time-point analysis

## R Package (pipelineRpack)

Main entrypoints:

- `pipelineRpack::run_rt_pct_strain_source_unified()`
- `pipelineRpack::rt_pct_strain_source_unified_cli()`

Install locally from package directory:

```bash
cd pipeline_v2/pipeline_Rpack
R CMD INSTALL .
```

Install from GitHub (subdirectory package):

```r
install.packages("remotes")
remotes::install_github("BJ-Chen-Eric/RTPCR", subdir = "pipeline_v2/pipeline_Rpack")
```

Private repo:

```r
Sys.setenv(GITHUB_PAT = "<your_token>")
remotes::install_github("BJ-Chen-Eric/RTPCR", subdir = "pipeline_v2/pipeline_Rpack")
```

Run package function:

```r
pipelineRpack::run_rt_pct_strain_source_unified(
  raw_dir = "RTpcr/raw",
  pattern = "20260203b73lox.csv",
  out_name = "output_test",
  analysis_mode = "auto",
  base_dir = getwd()
)
```

Run package CLI wrapper:

```bash
Rscript pipeline_v2/pipeline_Rpack/inst/scripts/rt_pct_strain_source_unified.R --help
```

## Get The Repository

Clone with SSH:

```bash
git clone git@github.com:BJ-Chen-Eric/RTPCR.git
cd RTPCR
```

From here, you can run the pipeline commands shown below.

## Run

From repository root:

```bash
Rscript pipeline_v2/work_pp/rt_pct_strain_source_unified.R --help
```

Or from `pipeline_v2/`:

```bash
Rscript work_pp/rt_pct_strain_source_unified.R --help
```

## Required R Packages

Install these packages before running:

```r
install.packages(c(
  "dplyr", "purrr", "fs", "seqinr", "stringr", "data.table", "ggplot2",
  "tidyr", "tibble", "openxlsx", "cowplot", "scales"
))
```

Notes:

- The pipeline checks required packages at runtime.
- Auto-install is disabled by default.
- To allow auto-install fallback, set `RT_AUTO_INSTALL_PKGS=true`.

## Input Data

Expected raw CSV format:

- qPCR export with a header row starting with `Well`
- columns containing sample name, sample type, gene, and Ct

Sample naming:

- Default parsing mode: `auto`
- Optional strict mode: `strain_time` (samples like `strain_0h`, `strain_1h`)

## Basic Usage

Single file:

```bash
Rscript pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
  --files "RTpcr/raw/example.csv" \
  --out-name run_demo
```

Multiple files:

```bash
Rscript pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
  --files "RTpcr/raw/file1.csv,RTpcr/raw/file2.csv" \
  --out-name run_merge
```

## Multi-Time Analysis

Use `analysis-mode multi` with `multi-compare-style all_time`:

```bash
Rscript pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
  --files "RTpcr/raw/file1.csv,RTpcr/raw/file2.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-name run_multi
```

## Important Flags

- `--raw-dir PATH`: directory containing raw CSV files
- `--files CSV`: comma-separated list of file paths (relative paths supported)
- `--out-name NAME`: output folder name under output parent
- `--control-gene GENE`: control gene (auto-detect: `5.8S`, then `TUB`)
- `--control-time T`: control time (default `0h`)
- `--analysis-mode auto|two|multi`
- `--multi-compare-style pairwise|all_time` (used in `multi`)
- `--sample-name-mode auto|strain_time`
- `--exclude-samples "sample_a,sample_b"`
- `--max-sample N`: outlier-trim sample cap (default `3`)
- `--perform-calibration true|false`
- `--calibration-sample`, `--calibration-gene`, `--calibration-time`
- `--sample-order "b73,lox10,opr78"`

## Default Paths

If not provided by flags:

- raw input defaults to `RTpcr/raw/` (resolved relative to script location)
- output parent defaults to `RTpcr/output/`

Result directory:

- `RTpcr/output/<out-name>/`

## Outputs

Main outputs:

- `all_raw_results.xlsx`
- `all_raw_time_results.xlsx`
- `all_time_results.xlsx` (when generated)
- `figure/*.png`

## Common Troubleshooting

- Error: missing control gene for some sample/time.
  Cause: target genes exist but no matching control-gene rows in the same sample/time.
  Fix: include the control-gene file in `--files`, or set the correct `--control-gene`.

- Error: no time token detected.
  Cause: sample names do not include recognizable time labels (for example `0h`, `1h`).
  Fix: rename samples or use `--sample-name-mode strain_time` if your naming matches that format.
