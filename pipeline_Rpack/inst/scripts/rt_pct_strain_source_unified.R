#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (!requireNamespace("pipelineRpack", quietly = TRUE)) {
  stop("Package 'pipelineRpack' is not installed. Install it first with: R CMD INSTALL .")
}

pipelineRpack::rt_pct_strain_source_unified_cli(args = args, base_dir = getwd())
