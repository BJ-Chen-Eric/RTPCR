ensure_package <- function(package) {
  # Default: do not auto-install packages in pipeline runs.
  # Set RT_AUTO_INSTALL_PKGS=true to enable install.packages fallback.
  auto_install <- tolower(Sys.getenv("RT_AUTO_INSTALL_PKGS", "false")) == "true"
  if (!requireNamespace(package, quietly = TRUE)) {
    if (auto_install) {
      install.packages(package, dependencies = TRUE)
    }
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(
        sprintf(
          "Missing required package '%s'. Install it manually or run with RT_AUTO_INSTALL_PKGS=true.",
          package
        )
      )
    }
  }
}

ensure_packages <- function(packages) {
  invisible(lapply(packages, ensure_package))
}


header_cleaning <- function(header_vector, pattern='/') {
  if (length(header_vector) == 0) {
    return(data.frame())
  }

  split_vec <- str_split(header_vector, pattern = pattern)
  max_cols <- max(lengths(split_vec))
  if (max_cols == 0) {
    return(data.frame())
  }

  padded <- lapply(split_vec, function(x) {
    length(x) <- max_cols
    x[is.na(x)] <- ""
    x
  })
  test <- do.call(what = rbind, padded) %>% as.data.frame(stringsAsFactors = FALSE)
  A <- data.frame(
    n_parts = lengths(split_vec),
    index = seq_len(nrow(test))
  )

  for (i in unique(A$n_parts)) {
    if (identical(i, ncol(test))) {next}
    test[A[A$n_parts %in% i, 'index'], (i + 1):ncol(test)] <- ''
  }
  return(test)
}



read_as_list <- function(path, sep='auto', prefix='', header='auto', cols=NA, file_type='txt') {
  out_list <- list()
  file_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
  # file_type <- sub(str_extract(file_path[1], pattern = '\\..*'), pattern = '.', replacement = '')
  file_path <- as_fs_path(file_path)
  if(file_type == 'tree')  {
    for(i in file_path) {
      out_list[[i]] <- ape::read.nexus(i)
    }
  }
  if(file_type == 'fasta')  {
    out_list <- file_path %>% map(.f = function(path){read.fasta(path, as.string = T, whole.header = T)})
  }
  if(file_type == 'RData')  {
    out_list <- file_path %>% map(.f = function(path){readRDS(path)})
  }
  if(file_type == 'csv')  {
    out_list <- file_path %>% map(.f = function(path){read.csv(path, header = T, sep = ',')})
  }
  if(file_type == 'txt')  {
    for(i in file_path) {
      out_list[[i]] <- fread(file = i, header = header, stringsAsFactors = F, sep = sep,  fill=TRUE) %>% as.data.frame()
    }
  }
  if(file_type == 'byline')  {
    for(i in file_path) {
      out_list[[i]] <- read.table(file = i, header = header, sep = "\n",) %>% as.data.frame()
    }
  }
  names(out_list) <- file_path
  pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
  for(i in 1:2)  {
    names(out_list) <- sub(names(out_list), pattern = pattern[i], replacement = '')
  }
  return(out_list)
}


table_DF <- function(x, prop=F) {
  A <- table(x) %>% as.data.frame(stringsAsFactors=F) %>% arrange(desc(Freq))
  if(isTRUE(prop))  {
    A <- table(x) %>% prop.table() %>% as.data.frame(stringsAsFactors=F) %>% arrange(desc(Freq))
  }
  return(A)
}



list2df <- function(list_input) {
  df <- do.call(list_input, what=rbind) %>% as.data.frame()
  return(df)
}



write_fasta <- 
  function(sequence, header, file_out) {
    out <- matrix(0,0,0) %>% as.data.frame()
    if(identical(stringr::str_extract(header[1], pattern = '>'), '>')) {
      next
    }else{header <- paste('>', header, sep = '')}
    
    out[seq(1, length(sequence)*2, 2), 1] <- header
    out[seq(2, length(sequence)*2, 2), 1] <- sequence
    
    write.table(out, file = file_out, quote = F, col.names = F, row.names = F)
  }


gg_theme <- ggplot2::theme(
  title = ggplot2::element_text(size = 20),
  axis.title.x = ggplot2::element_text(size = 28),
  axis.text.x = ggplot2::element_text(size = 22), # -6, angle=90,hjust=0.95,vjust=0.2
  axis.title.y = ggplot2::element_text(size = 28),
  axis.text.y = ggplot2::element_text(size = 22),
  plot.subtitle = ggplot2::element_text(size = 16),
  plot.caption = ggplot2::element_text(size = 30),
  legend.text = ggplot2::element_text(size = 20),
  legend.key.size = grid::unit(4, "lines"),
  legend.key.height = grid::unit(1, "cm"),
  strip.text = ggplot2::element_text(size = 20),
  strip.background = ggplot2::element_blank()
)



# helper: remove one outlier per group (farthest from the group median)
remove_outlier_by_median <- function(dt, group_col, value_col, min_group_size = 3) {
  stopifnot(is.data.table(dt))
  tmp <- copy(dt)

  tmp[, row_id := .I]
  tmp[, grp_size := .N, by = group_col]
  tmp[, dist_from_med := abs(get(value_col) - median(get(value_col), na.rm = TRUE)), by = group_col]
  tmp[, idx_to_drop := row_id[which.max(dist_from_med)], by = group_col]

  keep <- tmp[!(grp_size >= min_group_size & row_id == idx_to_drop)]
  keep[, c("grp_size", "dist_from_med", "idx_to_drop", "row_id") := NULL]
  keep[]
}
# geometric mean that skips non-positive values
geom_mean_pos <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  exp(mean(log(x)))
}



make_wide_reps <- function(df, control_gene, gene, include_data_resoure = TRUE, data_resoure_value = NULL) {
  # df=to_calc; include_data_resoure = TRUE; gene='MIR319'
  wide_out_gene <- data.frame()
  wide_out_control <- data.frame()
  id_cols <- c("p", "rep")
  if (include_data_resoure) {
    id_cols <- c(id_cols, "data_resoure")
  }
  if ((df$Gene %>% unique() %>% length()) == 1) {return(data.frame())}

  if (include_data_resoure && !("data_resoure" %in% colnames(df))) {
    return(data.frame())
  }

  wide_input <- df
  if (include_data_resoure) {
    wide_input <- wide_input %>% mutate(p = paste0(p, '|', rep, data_resoure))
  } else {
    wide_input <- wide_input %>% mutate(p = paste0(p, '|', rep))
  }

  wide <- wide_input %>%
    # filter(Ct < 40) %>%
    tidyr::pivot_wider(
      id_cols = tidyselect::all_of(id_cols),
      names_from = Gene,
      values_from = Ct
    ) %>%
    as.data.frame() %>%
    mutate(
      p  = str_remove(p, pattern = '\\|.*'),
      hr = str_extract(p, pattern = '_[0-9]+') %>% str_remove(pattern = '_')
    ) %>%
    arrange(p)

  if (!("data_resoure" %in% colnames(wide))) {
    wide$data_resoure <- if (!is.null(data_resoure_value)) data_resoure_value else "merge"
  }

  required_cols <- c("p", "rep", "data_resoure", control_gene, gene, "hr")
  if (!all(required_cols %in% colnames(wide))) {
    return(data.frame())
  }

  wide <- wide[, required_cols, drop = FALSE]
  # colnames(wide) <- c("p", "rep", "data_resoure", control_gene, gene, "hr")
  colnames(wide)[5] <- 'gene'
  colnames(wide)[5] <- gene
  colN <- colnames(wide)

  complete_flag <- FALSE
  if(all(!is.na(wide[, 5])) && all(!is.na(wide[, 4]))) {
    wide <- wide[complete.cases(wide), ]
    wide <- wide %>% as.data.frame() %>%
              mutate(dup_resource='complete data', note='complete data')
    colnames(wide) <- c("p", "rep", "data_resoure", "con", "gene", "hr", "dup_resource", "note")
    complete_flag <- TRUE

    # if(wide$data_resoure %>% unique() %>% length() != 1) {
    #   cross_wide <- list()
    #   wide_list <- wide %>% group_split(data_resoure) %>% as.list() %>% lapply(as.data.frame)
    #   for(w1 in seq_along(wide_list)) {
    #     for(w2 in seq_along(wide_list)) {
    #         if(w1 == w2) next
    #         cross_wide[[length(cross_wide)+1]] <- 
    #           cbind(wide_list[[w1]][, 1:4], wide_list[[w2]][, 5:8]) %>%
    #           mutate(data_resoure=wide_list[[w1]]$data_resoure[1], 
    #                  dup_resource=wide_list[[w1]]$data_resoure[1])
    #         cross_wide[[length(cross_wide)+1]] <- 
    #           cbind(wide_list[[w1]][, 1:4], wide_list[[w2]][, 5:8]) %>%
    #           mutate(data_resoure=wide_list[[w1]]$data_resoure[1], 
    #                  dup_resource=wide_list[[w1]]$data_resoure[1])
    #       }
    #     }
    #     wide <- c(wide_list, cross_wide) %>% list2df() 
    #   }
      
    }
  
  
  if(any(is.na(wide[, 5]))) {
    g <- wide$p %>% unique()
    sub <- wide 
    if(nrow(sub) <=1 ) return(data.frame())
    if(length(sub[which(is.na(sub[, 5])), 5]) < length(sub[which(!is.na(sub[, 5])), 5])) {
      withna <- sub[which(is.na(sub[, 5])), ]
      nonna <- sub[which(!is.na(sub[, 5])), ]
      sub <- rbind(nonna, withna[rep(seq_len(nrow(withna)), length.out = nrow(nonna)), ])
    }
    withna <- sub[which(is.na(sub[, 5])), ] 
    nonna <- sub[which(!is.na(sub[, 5])), ]
    sub <- merge(withna[, -5], nonna[, c('rep', 'data_resoure', gene)] %>% dplyr::rename(gene_resource=data_resoure), 
                by='rep', all.x=TRUE)[, c(colN, 'gene_resource')] %>% arrange(data_resoure)
    # sub[which(is.na(sub[, 5])), 5] <- sub[which(!is.na(sub[, 5])), 5]
    # final <- c(final, list(sub))
    sub <- rbind(sub, nonna %>% mutate(gene_resource=data_resoure))
    # wide_out_gene <- list2df(final)
    wide_out_gene <- sub[complete.cases(sub), ]
    # wide_out_gene <- wide_out_gene %>% as.data.frame() %>% mutate(note='duplicate gene')
    wide_out_gene <- wide_out_gene %>% as.data.frame() %>%
                          dplyr::rename(dup_resource=gene_resource) %>%
                          mutate(dup_resource=paste(dup_resource, gene, sep='_'), note='duplicate gene', pg=paste0(data_resoure, dup_resource, rep)) %>% 
                          distinct(pg, .keep_all = TRUE) %>% select(-pg)
  }else if (!any(is.na(wide[, 5]))){
    wide_out_gene <- data.frame()
  }
  

   
  if(any(is.na(wide[, 4]))) {
    final <- list()
    g <- wide$p %>% unique()
    sub <- wide 
    if(nrow(sub) <=1 ) return(data.frame())
    if(length(sub[which(is.na(sub[, 4])), 4]) < length(sub[which(!is.na(sub[, 4])), 4])) {
      withna <- sub[which(is.na(sub[, 4])), ]
      nonna <- sub[which(!is.na(sub[, 4])), ]
      sub <- rbind(nonna,   withna[rep(seq_len(nrow(withna)), length.out = nrow(nonna)), ])
    }

    withna <- sub[which(is.na(sub[, 4])), ] 
    nonna <- sub[which(!is.na(sub[, 4])), ]
    sub <- merge(withna[, -4], nonna[, c(2, 3, 4)] %>% dplyr::rename(control_resource=data_resoure), 
                by='rep', all.x=TRUE)[, c(colN, 'control_resource')] %>% arrange(data_resoure)
    # sub[which(is.na(sub[, 4])), 4] <- sub[which(!is.na(sub[, 4])), 4]
    # final <- c(final, list(sub))
    # wide_out_control <- list2df(final)
    sub <- rbind(sub, nonna %>% mutate(control_resource=data_resoure))
    wide_out_control <- sub[complete.cases(sub), ]
    # wide_out_control <- wide_out_control %>%  %>% mutate(note='duplicate gene')
    wide_out_control <- wide_out_control %>% as.data.frame() %>%
                          dplyr::rename(dup_resource=control_resource) %>%
                          mutate(dup_resource=paste(dup_resource, 'control', sep='_'), note='duplicate control', pg=paste0(data_resoure, dup_resource, rep)) %>% 
                          distinct(pg, .keep_all = TRUE) %>% select(-pg)
  }else if (!any(is.na(wide[, 4]))){
    wide_out_control <- data.frame()
  }
  

  if (!include_data_resoure && !is.null(data_resoure_value)) {
    wide$data_resoure <- data_resoure_value
  }
  if(isFALSE(complete_flag)) {
    wide <- rbind(wide_out_gene, wide_out_control) %>% as.data.frame()
  }
  return(wide)
      # p rep               data_resoure  5.8S MIR319 hr
}


write_workbook_from_list <- function(data_list, file_path) {
  wb <- openxlsx::createWorkbook()
  for (nm in names(data_list)) {
    openxlsx::addWorksheet(wb, sheetName = nm)
    openxlsx::writeData(wb, sheet = nm, x = data_list[[nm]])
  }

  openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE)
}


drop_outliers_by_group <- function(dt, max_sample=3) {
  # stopifnot(is.data.table(dt))
  # dt <- wide_reps_dt_c

  for(i in 1:20) {
    dt <- 
      dt %>%
        mutate(
          center  = median(delta),
          mad     = mad(delta, constant = 1),
          z_score = abs((delta - center) / mad),
          n       = sum(!is.na(delta))
        ) %>%
        slice(if (first(n) <= max_sample) row_number() else -which.max(z_score)) %>%
        as.data.frame()

  }
  
  dt
}

make_plot_list <- function(result_list, skip_small = FALSE) {
  plot_list <- list()
  for (i in names(result_list)) {
    plot_df <- result_list[[i]] %>%
      group_by(p) %>%
      summarise(
        mean_fc = mean(fold_change, na.rm = TRUE),
        sd_fc   = sd(fold_change, na.rm = TRUE),
        .groups = "drop"
      )

    if (skip_small && plot_df %>% nrow() < 2) next

    tt <- t.test(fold_change ~ p, data = result_list[[i]])
    pval <- tt$p.value

    stars <- case_when(
      pval <= 0.001 ~ "***",
      pval <= 0.01  ~ "**",
      pval <= 0.05  ~ "*",
      TRUE          ~ "ns"
    )

    y_max     <- max(plot_df$mean_fc + plot_df$sd_fc, na.rm = TRUE)
    bracket_y <- y_max * 1.05
    stars_y   <- y_max * 1.1

    plot_list[[i]] <-
      ggplot(plot_df, aes(x = p, y = mean_fc)) +
      geom_col(width = 0.6, fill = "#5fcd57ff") +
      geom_errorbar(
        aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc),
        width = 0.15,
        linewidth = 0.8
      ) +
      annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.7) +
      annotate("text", x = 1.5, y = stars_y, label = stars, size = 12) +
      labs(x = "Sample", y = "Fold Change", title = paste(i, "Fold Change", sep = ", ")) +
      theme_bw() + gg_theme +
      theme(axis.title.x = element_blank())
  }

  plot_list
}

save_plot_list <- function(plot_list, figure_dir) {
  for (i in names(plot_list)) {
    png(filename = file.path(figure_dir, paste0(i, "_result.png")), width = 2000, height = 1800, res = 300)
    print(plot_list[[i]])
    dev.off()
  }
}

make_combine_plot <- function(plot_list, sample_order) {
  miR <- names(plot_list) %>% header_cleaning('_') %>% pull(V2) %>% unique()
  time <- names(plot_list) %>% header_cleaning('_') %>% pull(V3) %>% unique()

  combine_plot <- list()
  for (t in time) {
    for (i in miR) {
      plots <- plot_list[grep(i, names(plot_list))]
      plots <- plots[grep(t, names(plots))]
      names(plots) <- names(plots) %>% header_cleaning('_') %>% pull(V1) %>% unique()
      plots <- plots[sample_order]
      plots <- plots[!vapply(plots, is.null, logical(1))]

      combine_plot[[paste(i, t, sep = "_")]] <- cowplot::plot_grid(plotlist = plots, nrow = 2)
    }
  }

  combine_plot
}

save_combine_plot <- function(combine_plot, figure_dir) {
  for (i in names(combine_plot)) {
    png(file.path(figure_dir, paste0('a_', i, "_combined.png")), width = 9000, height = 5000, res = 300)
    print(combine_plot[[i]])
    dev.off()
  }
}
