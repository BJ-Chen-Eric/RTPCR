#!/usr/bin/env Rscript

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, cmd_args)
  if (length(idx) > 0) {
    return(normalizePath(sub(file_arg, "", cmd_args[idx][1]), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

print_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript rt_pct_strain_source_unified.R [options]\n\n",
    "Options:\n",
    "  --raw-dir PATH          Input raw CSV directory\n",
    "  --files CSV             Comma-separated file names/paths (overrides --pattern)\n",
    "  --pattern REGEX         Regex used to select raw files\n",
    "  --parent PATH           Parent output directory\n",
    "  --out-name NAME         Output folder name under --parent\n",
    "  --out-dir PATH          Full output directory (overrides --parent/--out-name)\n",
    "  --control-gene GENE     Control gene name (auto: 5.8S then TUB)\n",
    "  --perform-calibration   true/false, apply cross-file Ct calibration (default: false)\n",
    "  --calibration-sample S  Sample_name used for calibration (auto: '*standard*' common across files)\n",
    "  --calibration-gene G    Gene used for calibration (default: control gene)\n",
    "  --calibration-time T    Time used for calibration (default: 0h)\n",
    "  --control-time T        Control time point for 2-time comparison (default: 0h)\n",
    "  --analysis-mode M       'auto', 'two', or 'multi' (default: auto)\n",
    "  --multi-compare-style S 'pairwise' or 'all_time' for multi mode (default: pairwise)\n",
    "  --sample-name-mode M    'auto' or 'strain_time' (default: auto)\n",
    "  --exclude-samples CSV   Comma-separated sample names to exclude from analysis\n",
    "  --max-sample N          max_sample for drop_outliers_by_group (default: 3)\n",
    "  --sample-order CSV      Optional order for combined plots; unspecified strains are still included\n",
    "  --help                  Show this help\n",
    sep = ""
  )
}

parse_cli_args <- function(args) {
  opts <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (startsWith(arg, "--")) {
      kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
      key <- kv[[1]]
      if (length(kv) > 1) {
        val <- paste(kv[-1], collapse = "=")
      } else {
        if (i == length(args) || startsWith(args[[i + 1]], "--")) {
          val <- TRUE
        } else {
          i <- i + 1
          val <- args[[i]]
        }
      }
      opts[[key]] <- val
    }
    i <- i + 1
  }
  opts
}

script_path <- get_script_path()
script_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)

resolve_first_existing_dir <- function(paths) {
  for (p in paths) {
    if (dir.exists(p)) return(normalizePath(p, mustWork = FALSE))
  }
  normalizePath(paths[[1]], mustWork = FALSE)
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (isTRUE(x)) return(TRUE)
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Invalid boolean value: ", x)
}

parse_int <- function(x, default) {
  if (is.null(x)) return(default)
  n <- suppressWarnings(as.integer(x))
  if (is.na(n)) stop("Invalid integer value: ", x)
  n
}

normalize_sample_name <- function(x) {
  tolower(trimws(as.character(x)))
}

normalize_time_label <- function(x) {
  tolower(trimws(as.character(x)))
}

normalize_header_name <- function(x) {
  x <- tolower(trimws(x))
  gsub("[^a-z0-9]", "", x)
}

match_first_col <- function(norm_names, candidates, label, required = TRUE) {
  idx <- which(norm_names %in% candidates)
  if (length(idx) == 0) {
    # Fallback: tolerate hidden header artifacts by substring matching.
    idx <- which(vapply(norm_names, function(nm) any(grepl(paste(candidates, collapse = "|"), nm)), logical(1)))
  }
  if (length(idx) == 0) {
    if (required) stop("Cannot find required column: ", label)
    return(NA_integer_)
  }
  idx[1]
}

read_qpcr_rows <- function(file_path) {
  raw_lines <- readLines(file_path, warn = FALSE)
  raw_lines <- sub("^\ufeff", "", raw_lines)

  header_idx <- grep("^Well(,|$)", raw_lines)
  if (length(header_idx) == 0) {
    stop("Cannot find qPCR data header that starts with 'Well' in file: ", file_path)
  }

  data_lines <- raw_lines[header_idx[1]:length(raw_lines)]
  data_lines <- gsub("OPR6,7", "OPR67", data_lines, fixed = TRUE)
  data_lines <- gsub("OPR6/7", "OPR67", data_lines, fixed = TRUE)
  data_lines <- gsub("LOX4,5", "LOX45", data_lines, fixed = TRUE)
  data_lines <- gsub("LOX4/5", "LOX45", data_lines, fixed = TRUE)
  data_lines <- gsub("hr$", "h", data_lines)
  data_lines <- gsub("H,", "h,", data_lines, fixed = TRUE)

  con <- textConnection(data_lines)
  on.exit(close(con), add = TRUE)
  df <- utils::read.csv(con, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  keep <- rowSums(!is.na(df) & trimws(as.character(df)) != "") > 0
  df <- df[keep, , drop = FALSE]

  norm_names <- normalize_header_name(colnames(df))
  idx_well <- match_first_col(norm_names, c("well"), "Well", required = FALSE)
  idx_sample <- match_first_col(norm_names, c("samplename"), "Sample name", required = FALSE)
  idx_type <- match_first_col(norm_names, c("sampletype"), "Sample type", required = FALSE)
  idx_dye <- match_first_col(norm_names, c("dye"), "Dye", required = FALSE)
  idx_gene <- match_first_col(norm_names, c("gene"), "Gene", required = FALSE)
  idx_ct <- match_first_col(norm_names, c("ct"), "Ct", required = FALSE)
  idx_mean <- match_first_col(norm_names, c("meanct", "mean"), "Mean Ct", required = FALSE)

  # Fallback for vendor CSVs with malformed/empty header names.
  if (is.na(idx_well)) {
    plate_like <- vapply(df, function(col) {
      vals <- trimws(as.character(col))
      vals <- vals[vals != "" & !is.na(vals)]
      if (length(vals) == 0) return(FALSE)
      mean(grepl("^[A-H][0-9]{1,2}$", vals)) > 0.3
    }, logical(1))
    if (any(plate_like)) {
      idx_well <- which(plate_like)[1]
    } else {
      idx_well <- 1
    }
  }

  pick_or_offset <- function(idx, offset) {
    if (!is.na(idx)) return(idx)
    cand <- idx_well + offset
    if (cand <= ncol(df)) return(cand)
    NA_integer_
  }

  idx_sample <- pick_or_offset(idx_sample, 1)
  idx_type <- pick_or_offset(idx_type, 2)
  idx_dye <- pick_or_offset(idx_dye, 3)
  idx_gene <- pick_or_offset(idx_gene, 4)
  idx_ct <- pick_or_offset(idx_ct, 5)
  if (is.na(idx_mean)) idx_mean <- pick_or_offset(idx_mean, 6)

  if (any(is.na(c(idx_well, idx_sample, idx_type, idx_dye, idx_gene, idx_ct)))) {
    stop(
      "Cannot map required qPCR columns in file: ", file_path,
      ". Parsed headers: ", paste(colnames(df), collapse = " | ")
    )
  }

  v2 <- rep("", nrow(df))
  if ("v2" %in% norm_names) {
    idx_v2 <- match_first_col(norm_names, c("v2"), "V2", required = FALSE)
    if (!is.na(idx_v2)) v2 <- df[[idx_v2]]
  }

  out <- data.frame(
    well = as.character(df[[idx_well]]),
    V2 = as.character(v2),
    Sample_name = as.character(df[[idx_sample]]),
    Sample_type = as.character(df[[idx_type]]),
    sg = as.character(df[[idx_dye]]),
    Gene = as.character(df[[idx_gene]]),
    Ct = as.character(df[[idx_ct]]),
    mean = if (!is.na(idx_mean)) as.character(df[[idx_mean]]) else NA_character_,
    stringsAsFactors = FALSE
  )

  out
}

safe_t_test_p <- function(df, group_col = "p") {
  grp <- unique(df[[group_col]])
  if (length(grp) < 2) return(NA_real_)
  pval <- tryCatch(
    t.test(fold_change ~ p, data = df, var.equal = FALSE, alternative = "two.sided")$p.value,
    error = function(e) NA_real_
  )
  pval
}

p_to_stars <- function(pval) {
  if (is.na(pval)) return("ns")
  dplyr::case_when(
    pval <= 0.001 ~ "***",
    pval <= 0.01  ~ "**",
    pval <= 0.05  ~ "*",
    TRUE          ~ "ns"
  )
}

write_workbook_if_nonempty <- function(data_list, path) {
  if (length(data_list) == 0) {
    message("[warn] Skip workbook with no sheets: ", path)
    return(invisible(NULL))
  }
  write_workbook_from_list(data_list, path)
  message("[ok] Wrote workbook: ", path)
}

save_plot_if_nonempty <- function(plot_obj, file_path, width, height, dpi = 300) {
  if (is.null(plot_obj)) return(invisible(NULL))
  ggsave(filename = file_path, plot = plot_obj, width = width, height = height, dpi = dpi, bg = "white")
}

safe_file_component <- function(x) {
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "unnamed", x)
}

blank_white_plot <- function() {
  ggplot() +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

make_square_plot_grid <- function(plots) {
  plots <- plots[!vapply(plots, is.null, logical(1))]
  if (length(plots) == 0) return(NULL)
  side <- ceiling(sqrt(length(plots)))
  target_n <- side * side
  if (length(plots) < target_n) {
    pads <- replicate(target_n - length(plots), blank_white_plot(), simplify = FALSE)
    plots <- c(plots, pads)
  }
  cowplot::plot_grid(plotlist = plots, ncol = side)
}

make_vertical_plot_grid <- function(plots) {
  plots <- plots[!vapply(plots, is.null, logical(1))]
  if (length(plots) == 0) return(NULL)
  cowplot::plot_grid(plotlist = plots, ncol = 1)
}

panel_name_from_result_key <- function(keys, mode = "auto", perspective = "miR") {
  meta <- header_cleaning(keys, "_")
  if (nrow(meta) == 0) return(character(0))
  if (ncol(meta) < 3) {
    return(if (perspective == "miR") meta$V1 else meta$V2)
  }

  if (mode == "multi") {
    if (perspective == "miR") {
      return(paste(meta$V1, meta$V3, sep = "_"))
    }
    return(paste(meta$V2, meta$V3, sep = "_"))
  }

  if (perspective == "miR") return(meta$V1)
  meta$V2
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    print_usage()
    return(invisible(NULL))
  }

  source(file.path(script_dir, "functions_rt.R"))
  ensure_packages(c("tidyr", "tibble", "openxlsx", "cowplot", "scales"))
  suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(data.table)
    library(ggplot2)
    library(tidyr)
    library(tibble)
    library(openxlsx)
    library(cowplot)
  })

  default_parent <- resolve_first_existing_dir(c(
    file.path(script_dir, "..", "..", "RTpcr", "output"),
    file.path(script_dir, "..", "RTpcr", "output")
  ))
  default_raw <- resolve_first_existing_dir(c(
    file.path(script_dir, "..", "..", "RTpcr", "raw"),
    file.path(script_dir, "..", "RTpcr", "raw")
  ))

  parent <- if (!is.null(args$parent)) args$parent else default_parent
  raw_dir <- if (!is.null(args$`raw-dir`)) args$`raw-dir` else default_raw
  files_arg <- if (!is.null(args$files)) args$files else NULL
  pattern <- if (!is.null(args$pattern)) args$pattern else "20260203b73lox.csv"
  out_name <- if (!is.null(args$`out-name`)) args$`out-name` else paste0("output_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  out_dir <- if (!is.null(args$`out-dir`)) args$`out-dir` else file.path(parent, out_name)

  control_gene_arg <- args$`control-gene`
  control_gene <- if (!is.null(control_gene_arg)) toupper(control_gene_arg) else NULL
  perform_calibration <- parse_bool(args$`perform-calibration`, default = FALSE)
  calibration_sample_user <- args$`calibration-sample`
  calibration_gene_user <- args$`calibration-gene`
  calibration_time_user <- args$`calibration-time`
  calibration_sample <- if (!is.null(calibration_sample_user)) normalize_sample_name(calibration_sample_user) else NULL
  calibration_gene <- if (!is.null(calibration_gene_user)) toupper(calibration_gene_user) else NULL
  calibration_time <- if (!is.null(calibration_time_user)) calibration_time_user else "0h"
  control_time <- if (!is.null(args$`control-time`)) args$`control-time` else "0h"
  control_time <- normalize_time_label(control_time)
  analysis_mode <- if (!is.null(args$`analysis-mode`)) tolower(trimws(args$`analysis-mode`)) else "auto"
  if (!(analysis_mode %in% c("auto", "two", "multi"))) {
    stop("--analysis-mode must be one of: auto, two, multi")
  }
  multi_compare_style <- if (!is.null(args$`multi-compare-style`)) tolower(trimws(args$`multi-compare-style`)) else "pairwise"
  if (!(multi_compare_style %in% c("pairwise", "all_time"))) {
    stop("--multi-compare-style must be one of: pairwise, all_time")
  }
  sample_name_mode <- if (!is.null(args$`sample-name-mode`)) tolower(trimws(args$`sample-name-mode`)) else "auto"
  if (!(sample_name_mode %in% c("auto", "strain_time"))) {
    stop("--sample-name-mode must be one of: auto, strain_time")
  }
  exclude_samples <- character(0)
  if (!is.null(args$`exclude-samples`)) {
    exclude_samples <- strsplit(args$`exclude-samples`, ",", fixed = TRUE)[[1]]
    exclude_samples <- tolower(trimws(exclude_samples))
    exclude_samples <- exclude_samples[nzchar(exclude_samples)]
  }
  max_sample <- parse_int(args$`max-sample`, default = 3)
  if (max_sample < 1) stop("--max-sample must be >= 1")

  sample_order <- character(0)
  if (!is.null(args$`sample-order`)) {
    sample_order <- strsplit(args$`sample-order`, ",", fixed = TRUE)[[1]]
    sample_order <- tolower(trimws(sample_order))
  }

  if (!dir.exists(raw_dir)) stop("raw_dir does not exist: ", raw_dir)

  selected_file_paths <- character(0)
  if (!is.null(files_arg)) {
    files_input <- strsplit(files_arg, ",", fixed = TRUE)[[1]]
    files_input <- trimws(files_input)
    files_input <- files_input[nzchar(files_input)]
    if (length(files_input) == 0) stop("--files was provided but empty")

    selected_file_paths <- vapply(files_input, function(x) {
      if (grepl("^/", x)) {
        return(normalizePath(x, mustWork = FALSE))
      }
      normalizePath(file.path(raw_dir, x), mustWork = FALSE)
    }, character(1))
    missing_files <- selected_file_paths[!file.exists(selected_file_paths)]
    if (length(missing_files) > 0) {
      stop("The following files from --files do not exist: ", paste(missing_files, collapse = ", "))
    }
  } else {
    files <- list.files(raw_dir, pattern = pattern)
    if (length(files) == 0) {
      stop("No files matched pattern '", pattern, "' in ", raw_dir)
    }
    selected_file_paths <- normalizePath(file.path(raw_dir, files), mustWork = FALSE)
  }
  selected_file_paths <- sort(selected_file_paths)

  figure_dir <- file.path(out_dir, "figure")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  message("[config] raw_dir: ", raw_dir)
  message("[config] files: ", paste(basename(selected_file_paths), collapse = ", "))
  if (is.null(files_arg)) message("[config] pattern: ", pattern)
  message("[config] out_dir: ", out_dir)
  message("[config] perform_calibration: ", perform_calibration)
  message("[config] control_time: ", control_time)
  message("[config] analysis_mode: ", analysis_mode)
  message("[config] multi_compare_style: ", multi_compare_style)
  message("[config] sample_name_mode: ", sample_name_mode)
  message("[config] exclude_samples: ", ifelse(length(exclude_samples) == 0, "<none>", paste(exclude_samples, collapse = ", ")))
  message("[config] max_sample: ", max_sample)

  filter_csv <- list()
  file_labels <- make.unique(tools::file_path_sans_ext(basename(selected_file_paths)))
  for (idx in seq_along(selected_file_paths)) {
    fp <- selected_file_paths[[idx]]
    file_label <- file_labels[[idx]]
    time_extract_pattern <- if (sample_name_mode == "strain_time") {
      "(?<=_)[0-9]+\\s*h(?:r)?(?:\\-2)?$"
    } else {
      "[0-9]+\\s*h(?:r)?(?:\\-2)?"
    }
    sample_remove_pattern <- if (sample_name_mode == "strain_time") {
      "_[0-9]+\\s*h(?:r)?(?:\\-2)?$"
    } else {
      "[_-]?[0-9]+\\s*h(?:r)?(?:\\-2)?"
    }

    A <- read_qpcr_rows(fp)

    A <- A %>%
      filter(Sample_name != "", !Sample_type %in% c("Empty", "NTC")) %>%
      mutate(
        sample_raw = tolower(trimws(Sample_name)),
        Time_raw = str_extract(sample_raw, pattern = time_extract_pattern),
        has_time_token = !is.na(Time_raw),
        Time_raw = str_replace_all(Time_raw, pattern = "\\s+", replacement = ""),
        Time_raw = str_replace(Time_raw, pattern = "hr", replacement = "h"),
        Time_raw = str_remove(Time_raw, pattern = "\\-2$"),
        Sample_name = str_remove(sample_raw, pattern = sample_remove_pattern),
        Sample_name = str_remove(Sample_name, pattern = "[-_]2$"),
        Sample_name = tolower(Sample_name),
        Sample_name = ifelse(Sample_name == "opr67", "opr78", Sample_name),
        Sample_name = str_remove(Sample_name, pattern = " "),
        Gene = toupper(Gene),
        Ct = ifelse(Ct == "No Ct", 100, Ct) %>% as.numeric(),
        Time = ifelse(is.na(Time_raw), control_time, normalize_time_label(Time_raw)),
        p = paste(Sample_name, Time, sep = "_"),
        p_filter = paste(Sample_name, Time, Gene, sep = "_")
      ) %>%
      select(Sample_name, has_time_token, Time_raw, Time, Ct, Gene, p, p_filter)

    if (!any(A$has_time_token, na.rm = TRUE)) {
      bad_examples <- A %>%
        distinct(Sample_name) %>%
        head(8) %>%
        pull(Sample_name) %>%
        paste(collapse = ", ")
      stop(
        "No time point token like '<number>h' was found in file: ",
        basename(fp),
        ". This 2-time-point pipeline requires explicit time labels in sample names. ",
        "If you already use 'strain_time' format, run with --sample-name-mode strain_time. ",
        "Example parsed sample names: ", bad_examples
      )
    }

    if (length(exclude_samples) > 0) {
      before_n <- nrow(A)
      A <- A %>% filter(!(Sample_name %in% exclude_samples))
      removed_n <- before_n - nrow(A)
      if (removed_n > 0) {
        message("[exclude] ", file_label, ": removed ", removed_n, " rows by --exclude-samples")
      }
    }

    A <- A %>% select(-has_time_token, -Time_raw)

    A <- A %>%
      group_by(p_filter) %>%
      mutate(rep = row_number()) %>%
      ungroup()

    filter_csv[[file_label]] <- A %>% as.data.frame()
  }

  if (is.null(control_gene)) {
    detected_control <- NULL
    for (cg in c("5.8S", "TUB")) {
      found <- any(vapply(filter_csv, function(x) any(x$Gene == cg, na.rm = TRUE), logical(1)))
      if (found) {
        detected_control <- cg
        break
      }
    }
    if (is.null(detected_control)) {
      stop(
        "Control gene was not specified and auto-detection failed: neither 5.8S nor TUB was found. ",
        "Use --control-gene explicitly."
      )
    }
    control_gene <- detected_control
    message("[auto] control_gene detected: ", control_gene)
  } else {
    found_control <- any(vapply(filter_csv, function(x) any(x$Gene == control_gene, na.rm = TRUE), logical(1)))
    if (!found_control) {
      stop("Specified --control-gene '", control_gene, "' was not found in selected files.")
    }
  }

  if (is.null(calibration_gene)) {
    calibration_gene <- control_gene
  }

  if (perform_calibration && is.null(calibration_sample)) {
    standard_by_file <- lapply(filter_csv, function(x) {
      x %>%
        filter(str_detect(Sample_name, regex("standard", ignore_case = TRUE))) %>%
        pull(Sample_name) %>%
        unique() %>%
        as.character()
    })

    has_standard_in_all <- all(vapply(standard_by_file, function(v) length(v) > 0, logical(1)))
    if (!has_standard_in_all) {
      stop(
        "Calibration requested but no '*standard*' sample was found in at least one input file. ",
        "Provide --calibration-sample explicitly."
      )
    }

    common_standard <- Reduce(intersect, standard_by_file)
    if (length(common_standard) == 0) {
      stop(
        "Calibration requested but there is no common '*standard*' sample across all files. ",
        "Provide --calibration-sample explicitly."
      )
    }
    calibration_sample <- sort(common_standard)[1]
    message("[auto] calibration_sample detected: ", calibration_sample)
  }

  message("[config] control_gene: ", control_gene)
  message("[config] calibration_sample: ", ifelse(is.null(calibration_sample), "<none>", calibration_sample))
  message("[config] calibration_gene: ", calibration_gene)
  message("[config] calibration_time: ", calibration_time)

  if (perform_calibration && length(filter_csv) > 1) {
    ref_name <- names(filter_csv)[1]
    ref_ct <- filter_csv[[ref_name]] %>%
      filter(
        Sample_name == calibration_sample,
        Time == calibration_time,
        Gene == calibration_gene,
        Ct != 100
      ) %>%
      pull(Ct) %>%
      mean(na.rm = TRUE)

    if (!is.finite(ref_ct)) {
      stop(
        "Calibration failed: no valid reference Ct in file '", ref_name,
        "' for sample=", calibration_sample, ", gene=", calibration_gene, ", time=", calibration_time
      )
    }

    for (k in 2:length(filter_csv)) {
      nm <- names(filter_csv)[k]
      cur_ct <- filter_csv[[nm]] %>%
        filter(
          Sample_name == calibration_sample,
          Time == calibration_time,
          Gene == calibration_gene,
          Ct != 100
        ) %>%
        pull(Ct) %>%
        mean(na.rm = TRUE)

      if (!is.finite(cur_ct)) {
        stop(
          "Calibration failed: no valid Ct in file '", nm,
          "' for sample=", calibration_sample, ", gene=", calibration_gene, ", time=", calibration_time
        )
      }

      offset <- cur_ct - ref_ct
      filter_csv[[nm]]$Ct <- filter_csv[[nm]]$Ct - offset
      message("[calibration] ", nm, " offset=", round(offset, 4), " (aligned to ", ref_name, ")")
    }
  } else if (perform_calibration && length(filter_csv) <= 1) {
    message("[calibration] skipped (need at least 2 files)")
  }

  cname <- Reduce(intersect, lapply(filter_csv, colnames))
  for (i in names(filter_csv)) {
    filter_csv[[i]] <- filter_csv[[i]][, cname]
  }

  all <- filter_csv %>%
    list2df() %>%
    tibble::rownames_to_column("data_resoure") %>%
    mutate(
      data_resoure = str_remove(data_resoure, pattern = "\\..*") %>% str_remove(pattern = "\\.[a-z]+$"),
      p_filter = paste(p_filter, data_resoure, rep, sep = "_")
    )

  # Fail early when a sample-time has target genes but no control gene.
  pair_check <- all %>%
    group_by(p) %>%
    summarise(
      has_control = any(Gene == control_gene, na.rm = TRUE),
      has_target = any(Gene != control_gene, na.rm = TRUE),
      .groups = "drop"
    )
  missing_control_pairs <- pair_check %>% filter(has_target, !has_control)
  if (nrow(missing_control_pairs) > 0) {
    bad_pairs <- head(missing_control_pairs$p, 12)
    stop(
      "Some sample/time groups have target genes but no control gene '", control_gene, "'. ",
      "Cannot compute deltaCt for these groups. Missing control for: ",
      paste(bad_pairs, collapse = ", "),
      if (nrow(missing_control_pairs) > 12) " ..." else "",
      ". Include matching control-gene rows (often from an additional file) or choose a different --control-gene."
    )
  }

  sample <- unique(all$Sample_name)
  miR_df <- all %>% filter(Gene != control_gene)
  miR <- unique(miR_df$Gene)
  miR <- miR[!miR %in% c(toupper(control_gene), "5.8S", "U6", "")]
  sample <- sample[!sample %in% c("")]
  sample <- sample[!str_detect(sample, regex("standard$", ignore_case = TRUE))]
  time_points <- sort(unique(all$Time))

  if (length(time_points) == 0) {
    stop("No time points were detected after parsing input files.")
  }
  if (!(control_time %in% time_points)) {
    stop(
      "Control time '", control_time, "' is not present in detected time points: ",
      paste(time_points, collapse = ", ")
    )
  }
  if (analysis_mode == "auto") {
    analysis_mode <- ifelse(length(time_points) == 2, "two", "multi")
    message("[auto] analysis_mode resolved to: ", analysis_mode)
  }
  if (analysis_mode != "multi" && multi_compare_style != "pairwise") {
    message("[warn] --multi-compare-style is only used in multi mode; falling back to pairwise.")
    multi_compare_style <- "pairwise"
  }
  if (analysis_mode == "two" && length(time_points) != 2) {
    stop(
      "analysis-mode=two requires exactly 2 time points. Detected ",
      length(time_points), ": ", paste(time_points, collapse = ", ")
    )
  }
  if (analysis_mode == "multi" && length(time_points) < 2) {
    stop("analysis-mode=multi requires at least 2 time points.")
  }
  compare_times <- setdiff(time_points, control_time)
  if (length(compare_times) == 0) {
    stop("No non-control time points detected. control_time=", control_time)
  }

  result_raw <- list()
  for (t in time_points) {
    for (s1 in sample) {
      for (g in miR) {
        to_calc <- all %>%
          filter(p %in% c(paste0(s1, "_", t)), Gene %in% c(g, control_gene)) %>%
          select(p, Ct, Gene, rep, data_resoure)

        wide_reps <- make_wide_reps(to_calc, control_gene, g, include_data_resoure = TRUE)
        if (nrow(wide_reps) == 0) next

        wide_reps$pg <- paste0(wide_reps$rep, wide_reps[, 4], wide_reps[, 5])
        wide_reps <- wide_reps %>% distinct(pg, .keep_all = TRUE) %>% select(-pg)
        wide_reps$hr <- ifelse(is.na(wide_reps$hr), 0, wide_reps$hr)

        colnames(wide_reps)[c(4, 5)] <- c("con", "gene")
        wide_reps <- wide_reps %>% mutate(hr = as.numeric(hr))

        wide_reps <- wide_reps %>%
          mutate(
            species = str_extract(p, pattern = ".*_"),
            delta = gene - con
          ) %>%
          select(-species)

        result_raw[[paste(s1, g, t, sep = "_")]] <- as.data.frame(wide_reps)
      }
    }
  }

  if (length(result_raw) == 0) {
    stop("No valid result rows generated. Check control gene and input data.")
  }

  all_time <- result_raw %>% list2df()
  meta <- rownames(all_time) %>%
    header_cleaning("_") %>%
    mutate(miR = V2) %>%
    select(-c(V1, V2, V3))

  all_time <- cbind(meta, all_time)
  all_time$strain <- all_time$p %>% str_extract(pattern = ".*_") %>% str_remove(pattern = "_")

  species <- unique(all_time$strain)
  miR <- unique(all_time$miR)

  result_time_list <- list()
  result_raw_time <- list()

  for (s1 in species) {
    for (g in miR) {
      target_times <- if (analysis_mode == "multi" && multi_compare_style == "all_time") {
        time_points
      } else {
        compare_times
      }

      for (t in target_times) {
        t_num <- suppressWarnings(as.numeric(str_extract(as.character(t), "[0-9]+")))
        c_num <- suppressWarnings(as.numeric(str_extract(control_time, "[0-9]+")))
        if (is.na(t_num) || is.na(c_num)) next

        hr_keep <- if (analysis_mode == "multi" && multi_compare_style == "all_time") {
          suppressWarnings(as.numeric(str_extract(time_points, "[0-9]+")))
        } else {
          c(c_num, t_num)
        }

        wide_reps_dt <- all_time %>%
          filter(strain == s1, miR %in% g, hr %in% hr_keep) %>%
          mutate(p_filter = paste(p, data_resoure, rep, sep = "_")) %>%
          select(p, data_resoure, dup_resource, note, rep, con, gene, hr)

        wide_reps_dt <- wide_reps_dt %>%
          mutate(
            species = str_extract(p, pattern = ".*_"),
            delta = gene - con
          ) %>%
          select(-species)

        if (analysis_mode == "multi" && multi_compare_style == "all_time") {
          result_raw_time[[paste(s1, g, "alltime", sep = "_")]] <- as.data.frame(wide_reps_dt)
        } else if (analysis_mode == "multi") {
          result_raw_time[[paste(s1, g, t, "hr", sep = "_")]] <- as.data.frame(wide_reps_dt)
        } else {
          result_raw_time[[paste(s1, g, sep = "_")]] <- as.data.frame(wide_reps_dt)
        }

        wide_reps_dt <- wide_reps_dt %>% filter(con < 40 & gene < 40)
        if (length(unique(wide_reps_dt$p)) < 2) next

        wide_reps_dt_c <- as.data.table(wide_reps_dt[wide_reps_dt$hr == c_num, ])
        if (nrow(wide_reps_dt_c) == 0) next
        if (length(unique(wide_reps_dt_c$data_resoure)) != 1) {
          wide_reps_dt_c <- wide_reps_dt_c %>% filter(data_resoure != "20260108tcp43")
        }
        if (nrow(wide_reps_dt_c) == 0) next

        wide_reps_dt_c$dis <- wide_reps_dt_c$delta - median(wide_reps_dt_c$delta)
        wide_reps_dt_c <- wide_reps_dt_c %>%
          arrange(abs(dis)) %>%
          mutate(
            center = median(delta),
            mad = mad(delta, constant = 1),
            z_score = (delta - center) / mad,
            n = sum(!is.na(delta))
          ) %>%
          arrange(abs(z_score)) %>%
          distinct(rep, .keep_all = TRUE) %>%
          select(-c(center, mad, z_score, n, dis))

        wide_reps_dt_s <- as.data.table(wide_reps_dt[wide_reps_dt$hr != c_num & wide_reps_dt$hr %in% hr_keep, ])
        if (nrow(wide_reps_dt_s) == 0) next
        if (length(unique(wide_reps_dt_s$data_resoure)) != 1) {
          wide_reps_dt_s <- wide_reps_dt_s %>% filter(data_resoure != "20260108tcp43")
        }
        if (nrow(wide_reps_dt_s) == 0) next

        wide_reps_dt_s$dis <- wide_reps_dt_s$delta - median(wide_reps_dt_s$delta)
        wide_reps_dt_s <- wide_reps_dt_s %>%
          arrange(abs(dis)) %>%
          mutate(
            center = median(delta),
            mad = mad(delta, constant = 1),
            z_score = abs((delta - center) / mad),
            n = sum(!is.na(delta))
          ) %>%
          arrange(abs(z_score)) %>%
          distinct(if (analysis_mode == "multi" && multi_compare_style == "all_time") paste(rep, hr) else rep, .keep_all = TRUE) %>%
          select(-c(center, mad, z_score, n, dis))

        key <- paste(s1, g, sep = "_")
        if (key %in% c("xxx")) {
          wide_reps_dt_c <- wide_reps_dt_c %>% filter(rep %in% c(1:3))
          wide_reps_dt_s <- wide_reps_dt_s %>% filter(rep %in% c(1:3))
        } else {
          wide_reps_dt_c <- drop_outliers_by_group(wide_reps_dt_c, max_sample = max_sample)
          if (analysis_mode == "multi" && multi_compare_style == "all_time") {
            wide_reps_dt_s <- wide_reps_dt_s %>%
              group_by(hr) %>%
              group_split() %>%
              lapply(function(df_hr) drop_outliers_by_group(as.data.frame(df_hr), max_sample = max_sample)) %>%
              list2df()
          } else {
            wide_reps_dt_s <- drop_outliers_by_group(wide_reps_dt_s, max_sample = max_sample)
          }
        }

        required_merge_cols <- c("p", "data_resoure", "dup_resource", "note", "rep", "con", "gene", "hr", "delta")
        for (nm in required_merge_cols) {
          if (!(nm %in% colnames(wide_reps_dt_c))) wide_reps_dt_c[[nm]] <- NA
          if (!(nm %in% colnames(wide_reps_dt_s))) wide_reps_dt_s[[nm]] <- NA
        }
        wide_reps <- dplyr::bind_rows(
          wide_reps_dt_c[, required_merge_cols, drop = FALSE],
          wide_reps_dt_s[, required_merge_cols, drop = FALSE]
        ) %>% as.data.frame()
        if (nrow(wide_reps) == 0) next

        wide_reps <- wide_reps %>%
          mutate(
            species = str_extract(p, pattern = ".*_"),
            delta = gene - con,
            pre_fold_change = 2^(-delta),
            control = mean(pre_fold_change[hr == c_num], na.rm = TRUE),
            fold_change = 2^(-delta) / control
          ) %>%
          select(-species)

        result <- wide_reps %>%
          group_by(p) %>%
          mutate(mean_fc = if_else(row_number() == 1, mean(fold_change, na.rm = TRUE), NA_real_)) %>%
          ungroup() %>%
          as.data.frame()

        result <- replace(result, is.na(result), "")
        if (analysis_mode == "multi" && multi_compare_style == "all_time") {
          result_time_list[[paste(s1, g, "alltime", sep = "_")]] <- result
        } else if (analysis_mode == "multi") {
          result_time_list[[paste(s1, g, t, "hr", sep = "_")]] <- result
        } else {
          result_time_list[[paste(s1, g, sep = "_")]] <- result
        }

        if (analysis_mode == "multi" && multi_compare_style == "all_time") {
          break
        }
      }
    }
  }

  write_workbook_if_nonempty(result_raw, file.path(out_dir, "all_raw_results.xlsx"))
  write_workbook_if_nonempty(result_raw_time, file.path(out_dir, "all_raw_time_results.xlsx"))
  write_workbook_if_nonempty(result_time_list, file.path(out_dir, "all_time_results.xlsx"))

  if (length(result_time_list) == 0) {
    message("[warn] No time-comparison result sheets were generated; using fallback plotting from raw delta values.")

    fallback_df <- all_time %>%
      filter(con < 40, gene < 40) %>%
      mutate(
        fold_change = 2^(-delta),
        hr = ifelse(is.na(hr), 0, hr),
        hr = as.character(hr)
      )

    if (nrow(fallback_df) > 0) {
      fallback_sum <- fallback_df %>%
        group_by(strain, miR, hr) %>%
        summarise(
          mean_fc = mean(fold_change, na.rm = TRUE),
          sd_fc = sd(fold_change, na.rm = TRUE),
          .groups = "drop"
        )

      fallback_strain_plots <- list()
      for (st in unique(fallback_sum$strain)) {
        sub <- fallback_sum %>% filter(strain == st)
        if (nrow(sub) == 0) next
        sub <- sub %>%
          mutate(
            panel_key = paste(miR, paste0("hr=", hr), sep = "_")
          ) %>%
          arrange(miR, suppressWarnings(as.numeric(hr)), panel_key)
        sub$panel_key <- factor(sub$panel_key, levels = unique(sub$panel_key))

        ymax <- max(sub$mean_fc + ifelse(is.na(sub$sd_fc), 0, sub$sd_fc), na.rm = TRUE)
        if (!is.finite(ymax) || ymax <= 0) ymax <- 1

        fallback_strain_plots[[st]] <-
          ggplot(sub, aes(x = panel_key, y = mean_fc)) +
          geom_col(width = 0.6, fill = "#d57500") +
          geom_errorbar(
            aes(ymin = pmax(0, mean_fc - ifelse(is.na(sd_fc), 0, sd_fc)),
                ymax = mean_fc + ifelse(is.na(sd_fc), 0, sd_fc)),
            width = 0.15,
            linewidth = 0.8
          ) +
          scale_y_continuous(limits = c(0, ymax * 1.2), breaks = scales::pretty_breaks(n = 5)) +
          labs(
            x = "miR and Time",
            y = "2^-deltaCt",
            title = paste("Fallback (Strain)", toupper(st))
          ) +
          theme_bw() + gg_theme +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }

      for (nm in names(fallback_strain_plots)) {
        nm_safe <- safe_file_component(nm)
        save_plot_if_nonempty(
          fallback_strain_plots[[nm]],
          file.path(figure_dir, paste0("fallback_strain_", nm_safe, ".png")),
          width = 10, height = 7
        )
      }

      ordered_strains <- unique(c(sample_order, names(fallback_strain_plots)))
      ordered_strains <- ordered_strains[ordered_strains %in% names(fallback_strain_plots)]
      merge_plots <- fallback_strain_plots[ordered_strains]
      merge_plots <- merge_plots[!vapply(merge_plots, is.null, logical(1))]
      if (length(merge_plots) > 0) {
        merged_plot <- make_square_plot_grid(merge_plots)
        save_plot_if_nonempty(
          merged_plot,
          file.path(figure_dir, "fallback_strain_merged_cowplot.png"),
          width = 16, height = 9
        )
      }

      message("[ok] Wrote fallback strain figures: ", length(fallback_strain_plots))
    } else {
      message("[warn] Fallback plotting skipped: no valid rows after Ct filters.")
    }

    message("[done] Pipeline completed.")
    message("[done] Excel outputs in: ", out_dir)
    message("[done] Figure outputs in: ", figure_dir)
    return(invisible(NULL))
  }

  plot_df_list <- list()
  for (i in names(result_time_list)) {
    plot_df <- result_time_list[[i]] %>%
      group_by(p) %>%
      summarise(
        mean_fc = mean(fold_change, na.rm = TRUE),
        sd_fc = sd(fold_change, na.rm = TRUE),
        .groups = "drop"
      )

    if (nrow(plot_df) < 2) next
    plot_df$sd_fc[is.na(plot_df$sd_fc)] <- 0

    stars <- NA_character_
    if (nrow(plot_df) == 2) {
      pval <- safe_t_test_p(result_time_list[[i]], group_col = "p")
      stars <- p_to_stars(pval)
    }

    plot_df$y_max <- max(plot_df$mean_fc + plot_df$sd_fc, na.rm = TRUE)
    plot_df$bracket_y <- plot_df$y_max * 1.05
    plot_df$stars_y <- plot_df$y_max * 1.1
    plot_df$stars <- stars

    plot_df_list[[i]] <- plot_df
  }

  plot_list <- list()
  y_limits_by_mir <- list()
  for (nm in names(plot_df_list)) {
    meta_nm <- header_cleaning(nm, "_")
    if (nrow(meta_nm) == 0 || !("V2" %in% colnames(meta_nm))) next
    mir_nm <- as.character(meta_nm$V2[1])
    y_top <- max(plot_df_list[[nm]]$y_max, na.rm = TRUE)
    y_bottom <- min(plot_df_list[[nm]]$mean_fc - plot_df_list[[nm]]$sd_fc, na.rm = TRUE)
    if (!is.finite(y_top)) y_top <- 1
    if (!is.finite(y_bottom)) y_bottom <- -0.5
    if (is.null(y_limits_by_mir[[mir_nm]])) {
      y_limits_by_mir[[mir_nm]] <- c(y_bottom, y_top)
    } else {
      y_limits_by_mir[[mir_nm]] <- c(
        min(y_limits_by_mir[[mir_nm]][1], y_bottom, na.rm = TRUE),
        max(y_limits_by_mir[[mir_nm]][2], y_top, na.rm = TRUE)
      )
    }
  }

  for (i in names(plot_df_list)) {
    plot_df <- plot_df_list[[i]]
    if (is.null(plot_df) || nrow(plot_df) < 2) next

    mir <- i %>% header_cleaning("_") %>% pull(V2)
    mir <- as.character(mir[1])
    limits <- y_limits_by_mir[[mir]]
    if (is.null(limits)) {
      ymax <- max(plot_df$mean_fc + plot_df$sd_fc, na.rm = TRUE)
      min_y <- min(plot_df$mean_fc - plot_df$sd_fc, na.rm = TRUE)
    } else {
      min_y <- limits[1]
      ymax <- limits[2]
    }
    if (is.na(min_y) || min_y > -0.5) min_y <- -0.5
    if (!is.finite(ymax) || ymax <= 0) ymax <- max(plot_df$mean_fc + plot_df$sd_fc, na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1

    plot_df <- plot_df %>%
      mutate(
        time_num = suppressWarnings(as.numeric(str_extract(p, "[0-9]+")))
      ) %>%
      arrange(ifelse(str_detect(p, paste0("_", control_time, "$")), -Inf, time_num), p)
    plot_df$p <- factor(plot_df$p, levels = plot_df$p)

    p <-
      ggplot(plot_df, aes(x = p, y = mean_fc)) +
      geom_col(width = 0.6, fill = "#d57500") +
      geom_errorbar(
        aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc),
        width = 0.15,
        linewidth = 0.8
      ) +
      scale_y_continuous(limits = c(min_y, ymax * 1.2), breaks = scales::pretty_breaks(n = 5)) +
      labs(x = "Sample", y = "Fold Change", title = paste(i, "Fold Change", sep = ", ")) +
      theme_bw() + gg_theme +
      theme(axis.title.x = element_blank())

    if (nrow(plot_df) == 2 && !is.na(plot_df$stars[1])) {
      p <- p +
        annotate("segment", x = 1, xend = 2, y = plot_df$bracket_y[1], yend = plot_df$bracket_y[1], linewidth = 0.7) +
        annotate("text", x = 1.5, y = plot_df$stars_y[1], label = plot_df$stars[1], size = 12)
    } else if (nrow(plot_df) > 2) {
      control_label <- plot_df$p[which(str_detect(as.character(plot_df$p), paste0("_", control_time, "$")))[1]]
      if (!is.na(control_label)) {
        stars_map <- c()
        for (p_lab in as.character(plot_df$p)) {
          if (identical(p_lab, as.character(control_label))) next
          sub_df <- result_time_list[[i]] %>% filter(p %in% c(as.character(control_label), p_lab))
          pval <- safe_t_test_p(sub_df, group_col = "p")
          stars_map[p_lab] <- p_to_stars(pval)
        }
        if (length(stars_map) > 0) {
          stars_df <- plot_df %>%
            filter(as.character(p) %in% names(stars_map)) %>%
            mutate(stars = unname(stars_map[as.character(p)]))
          p <- p + annotate(
            "text",
            x = stars_df$p,
            y = stars_df$stars_y + 0.05,
            label = stars_df$stars,
            size = 12
          )
        }
      }
    }

    plot_list[[i]] <- p
  }

  if (length(plot_list) == 0) {
    message("[warn] No valid figures were generated from all_time_results (insufficient replicates after filtering).")
    message("[done] Pipeline completed.")
    message("[done] Excel outputs in: ", out_dir)
    message("[done] Figure outputs in: ", figure_dir)
    return(invisible(NULL))
  }

  for (nm in names(plot_list)) {
    nm_safe <- safe_file_component(nm)
    save_plot_if_nonempty(plot_list[[nm]], file.path(figure_dir, paste0(nm_safe, "_result.png")), width = 16, height = 9)
  }

  miR_for_combine <- names(plot_list) %>% header_cleaning("_") %>% pull(V2) %>% unique()
  combine_plot_by_mir <- list()
  combine_plot_by_mir_n <- list()
  for (m in miR_for_combine) {
    plots <- plot_list[grep(m, names(plot_list))]
    names(plots) <- panel_name_from_result_key(names(plots), mode = analysis_mode, perspective = "miR")
    if (length(sample_order) > 0) {
      order_final <- unique(c(sample_order, names(plots)))
      order_final <- order_final[order_final %in% names(plots)]
      plots <- plots[order_final]
    }
    combine_plot_by_mir_n[[m]] <- length(plots)
    combine_plot_by_mir[[m]] <- make_vertical_plot_grid(plots)
    if (is.null(combine_plot_by_mir[[m]])) next
  }

  for (nm in names(combine_plot_by_mir)) {
    nm_safe <- safe_file_component(nm)
    h <- max(4, 4 * combine_plot_by_mir_n[[nm]])
    save_plot_if_nonempty(combine_plot_by_mir[[nm]], file.path(figure_dir, paste0("Gene_", nm_safe, "_combined.png")), width = 16, height = h)
  }

  strain_for_combine <- names(plot_list) %>% header_cleaning("_") %>% pull(V1) %>% unique()
  combine_plot_by_strain <- list()
  combine_plot_by_strain_n <- list()
  for (s in strain_for_combine) {
    plots <- plot_list[grep(paste0("^", s, "_"), names(plot_list))]
    names(plots) <- panel_name_from_result_key(names(plots), mode = analysis_mode, perspective = "strain")
    combine_plot_by_strain_n[[s]] <- length(plots)
    combine_plot_by_strain[[s]] <- make_vertical_plot_grid(plots)
    if (is.null(combine_plot_by_strain[[s]])) next
  }

  for (nm in names(combine_plot_by_strain)) {
    nm_safe <- safe_file_component(nm)
    h <- max(4, 4 * combine_plot_by_strain_n[[nm]])
    save_plot_if_nonempty(combine_plot_by_strain[[nm]], file.path(figure_dir, paste0("strain_", nm_safe, "_combined.png")), width = 16, height = h)
  }

  message("[done] Pipeline completed.")
  message("[done] Excel outputs in: ", out_dir)
  message("[done] Figure outputs in: ", figure_dir)
}

main()
