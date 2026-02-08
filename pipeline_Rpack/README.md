# pipelineRpack

RT-qPCR unified workflow packaged as an R package.

## Main entrypoints

- `pipelineRpack::run_rt_pct_strain_source_unified()`
- `pipelineRpack::rt_pct_strain_source_unified_cli()`

Core logic was migrated from:

- `work_pp/rt_pct_strain_source_unified.R`
- `work_pp/functions_rt.R`

## Install locally

From the package directory:

```bash
R CMD INSTALL .
```

Or in R:

```r
install.packages(".", repos = NULL, type = "source")
```

## Install from GitHub

Public repo:

```r
install.packages("remotes")
remotes::install_github("BJ-Chen-Eric/RTPCR", subdir = "<path/to/package/in-repo>")
```

Private repo:

```r
Sys.setenv(GITHUB_PAT = "<your_token>")
remotes::install_github("BJ-Chen-Eric/RTPCR", subdir = "<path/to/package/in-repo>")
```

If the package is at the repository root, omit `subdir`.

## Run as package function

```r
pipelineRpack::run_rt_pct_strain_source_unified(
  raw_dir = "RTpcr/raw",
  pattern = "20260203b73lox.csv",
  out_name = "output_test",
  analysis_mode = "auto",
  base_dir = getwd()
)
```

You can still pass `args = list(...)` if you prefer. Underscore names are supported (for example `raw_dir`, `analysis_mode`) and are mapped internally to CLI-style keys.

## Run with CLI wrapper script

```bash
Rscript inst/scripts/rt_pct_strain_source_unified.R --help
Rscript inst/scripts/rt_pct_strain_source_unified.R --raw-dir RTpcr/raw --pattern '20260203b73lox.csv'
```

## Notes

- Original script versions are retained under `work_pp/` for reference.
- Package build excludes `work_pp/` and `source/` via `.Rbuildignore`.
