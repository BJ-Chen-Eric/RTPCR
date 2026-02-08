#!/usr/bin/env Rscript
source('/Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/functions_rt.R')


parent <- '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/output'
raw_dir    <- "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw"
# files <- list.files(raw_dir, pattern = '20251219localmir319|20251222mir159local|20251223mir159lical2')
# files <- list.files(raw_dir, pattern = '20251216ws319|20251217|20251217mir159')
# files <- list.files(raw_dir, pattern = '20251218mir319target|20251219opr2_lox10_chd_pr1')
# files <- list.files(raw_dir, pattern = '20260108tcp43|20251222b73target_sys|20251218mir319target|20251219opr2_lox10_chd_pr1')
files <- list.files(raw_dir, pattern = '20260206mir319b73.csv|20260206mir3192.csv')
# files <- list.files(raw_dir, pattern = '20260108tcp43')
out_dir <- file.path(parent, "/output_20260206mir319b73_merge/")
control_gene <- toupper("5.8s")


figure_dir <- file.path(out_dir, "figure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)



filter_csv <- list()
for (i in files) {
  A <- read.table(file = file.path(raw_dir, i), header = TRUE, sep = "\n") %>% as.data.frame()
  colnames(A) <- 'line'
  A <- A %>%
    mutate(
      line = str_replace(line, pattern = 'OPR6,7', replacement = 'OPR67'),
      line = str_replace(line, pattern = 'OPR6/7', replacement = 'OPR67'),
      line = str_replace(line, pattern = 'LOX4,5', replacement = 'LOX45'),
      line = str_replace(line, pattern = 'LOX4/5', replacement = 'LOX45'),
      line = str_replace(line, pattern = 'hr$',    replacement = 'h'),
      line = str_replace(line, pattern = 'H,',    replacement = 'h,')
    )

  A <- A$line %>% header_cleaning(',')
  colnames(A) <- c('well', 'V2', 'Sample_name', 'Sample_type', 'sg', 'Gene', 'Ct', 'mean')
  A <- A %>%
    filter(Sample_name != '', !Sample_type %in% c('Empty', 'NTC')) %>%
    mutate(
      p          = paste(Sample_name, Sample_type, sep = '') %>% str_remove(pattern = 'Unknown'),
      Time       = str_extract(p, pattern = '[0-9]h$|[0-9]h\\-2$') %>% str_remove(pattern = '\\-2$'),
      Sample_name = str_remove(p, pattern = '_[0-9]h$|_[0-9]h\\-2$|_-2$'),
      Sample_name = toupper(Sample_name),
      # Sample_name = ifelse(Sample_name %in% c(''), c(''), tolower(Sample_name)),
      Sample_name = ifelse(Sample_name %in% c('B73'), c('B73'), tolower(Sample_name)),
      Sample_name = ifelse(Sample_name == 'opr67', 'opr78', Sample_name),
      Sample_name = str_remove(Sample_name, pattern = ' '),
      Gene       = toupper(Gene),
      Ct         = ifelse(Ct == 'No Ct', 100, Ct) %>% as.numeric(),
      Time       = ifelse(is.na(Time), '0h', Time),
      p          = paste(Sample_name, Time, sep = '_'),
      p_filter   = paste(Sample_name, Time, Gene, sep = '_')
    ) %>% 
    select(Sample_name, Time, Ct, Gene, p, p_filter)

    A <- A %>%
          group_by(p_filter) %>%
          mutate(rep = row_number()) %>%
          ungroup()

  # A_filtered <- remove_outlier_by_median(as.data.table(A), group_col = "p_filter", value_col = "Ct")
  filter_csv[[i]] <- A %>% as.data.frame()
}
names(filter_csv) <- str_remove(names(filter_csv), pattern = '.*\\/') %>% str_remove(pattern = '\\.[a-z]+$')


calli <- lapply(filter_csv, FUN=function(x) {x %>% filter(Sample_name == '5standard', Time=='0h', Gene=='MIR319', Ct!=100) %>% pull(Ct) %>% mean()})  # check standard
# calli <- lapply(filter_csv, FUN=function(x) {x %>% filter(Sample_name == 'lox4', Gene=='5.8S', Ct!=100) %>% pull(Ct) %>% mean()})  # check standard
calli <- calli[[2]]-calli[[1]]
filter_csv[[2]]$Ct <- filter_csv[[2]]$Ct - calli  # adjust standard difference


cname <- Reduce(intersect, filter_csv %>% lapply(colnames))
for (i in names(filter_csv)) {
  filter_csv[[i]] <- filter_csv[[i]][, cname]
}


all <- filter_csv %>% list2df() %>% tibble::rownames_to_column('data_resoure') %>% 
      mutate(data_resoure = str_remove(data_resoure, pattern = '\\..*') %>% str_remove(pattern = '\\.[a-z]+$'), 
             p_filter = paste(p_filter, data_resoure, rep, sep = '_'))

table(all$Sample_name, all$Time, all$Gene)
sample <- all$Sample_name %>% unique()

# con    <- all %>% filter(Gene == control_gene)
# con <- rbind(con, con %>% mutate(data_resoure = '20251219opr2_lox10_chd_pr1'))  # add merged control

miR_df <- all %>% filter(Gene != control_gene)
check <- table(miR_df$Sample_name, miR_df$Gene) %>% as.data.frame.matrix() %>% as.data.frame() 
table(miR_df$Sample_name, miR_df$Gene)
table(miR_df$Sample_name, miR_df$Time, miR_df$Gene)
miR <- miR_df$Gene %>% unique()
miR <- miR[!miR %in% c(toupper(control_gene))]  # remove mRNA from miR list

sample <- sample[!sample %in% c('1standard', '1', '11standard', '', 'b73standard', '5standard')]  # remove control from sample list
miR <- miR[!miR %in% c('5.8S', 'U6', '')]  # remove control from miR list
time <- all$Time %>% unique()



g='MIR159'
s1='opr2'
t='1h'
# opr2_mir319_0h
result_raw <- result_list <- list()
for(t in time) {
  # t <- '0h'
  for (s1 in sample) {
    # s1 <- sample[1]
    for (g in miR) {
      # g <- miR[3]
      to_calc <- all %>% filter(p %in% c(paste0(s1, '_', t)), Gene %in% c(g, control_gene)) %>% 
                 select(p, Ct, Gene, rep, data_resoure)
      # to_calc <- rbind(
      #     con    %>% filter(p %in% c(paste0('B73', '_', t), paste0(s1, '_', t))) %>% select(p, Ct, Gene, rep, data_resoure),
      #     miR_df %>% filter(p %in% c(paste0('B73', '_', t), paste0(s1, '_', t)), Gene %in% c(g)) %>% select(p, Ct, Gene, rep, data_resoure)
          # con    %>% filter(p %in% c(paste0('', '_', t), paste0(s1, '_', t))) %>% select(p, Ct, Gene, rep, data_resoure),
          # miR_df %>% filter(p %in% c(paste0('', '_', t), paste0(s1, '_', t)), Gene %in% c(g)) %>% select(p, Ct, Gene, rep, data_resoure)
        # )
      # if(nrow(to_calc) != 12) next

      wide_reps <- make_wide_reps(to_calc, control_gene, g, include_data_resoure = TRUE)
      if(wide_reps %>% nrow() == 0) next
      if(nrow(wide_reps) < (nrow(to_calc)-3)/2) stop("Error in wide reps")
      # if(unique(wide_reps$p) %>% length() < 2) {
      #   if(nrow(wide_reps) != 0) {
      #     fine <- wide_reps$p[1]
      #     wide_reps_re <- to_calc %>% filter(p != fine) %>%
      #       select(p, rep, Gene, Ct) %>%
      #       make_wide_reps(include_data_resoure = FALSE, data_resoure_value = "merge")
      #     wide_reps_re <- wide_reps_re[, colnames(wide_reps)]
      #     wide_reps <- rbind(wide_reps, wide_reps_re)
      #   }else {
      #      wide_reps_re <- to_calc %>%
      #       select(p, rep, Gene, Ct) %>%
      #       make_wide_reps(include_data_resoure = FALSE, data_resoure_value = "merge")
          
      #     wide_reps <- wide_reps_re
      #   }
      # }
      wide_reps$pg <- paste0(wide_reps$rep, wide_reps[, 4], wide_reps[, 5])
      wide_reps <- wide_reps %>% distinct(pg, .keep_all = TRUE) %>% select(-pg)
      wide_reps$hr <- ifelse(is.na(wide_reps$hr), 0, wide_reps$hr)

      # needed_cols <- c("p", "rep", control_gene, g, "hr", 'data_resoure', 'note')
      # if (!all(needed_cols %in% colnames(wide_reps))) next

      # wide_reps <- wide_reps[, needed_cols]
      colnames(wide_reps)[c(4, 5)] <- c("con", "gene")
      wide_reps <- wide_reps %>% mutate(hr = as.numeric(hr))

      wide_reps <- wide_reps %>%
        mutate(
          species    = str_extract(p, pattern = '.*_'),
          delta       = gene - con
        ) %>% select(-species)

      wide_reps_dt <- as.data.frame(wide_reps)

      result_raw[[paste(s1, g, t, sep = "_")]]  <- wide_reps_dt %>% as.data.frame()

      # wide_reps_dt <- wide_reps_dt %>% filter(con < 40 & gene < 40)  # remove Ct > 35
      # if(length(unique(wide_reps_dt$p)) < 2) next
      
      # wide_reps_dt_c <- as.data.table(wide_reps_dt[grep(wide_reps_dt$p, pattern = '^_'), ])
      # wide_reps_dt_c$dis <- wide_reps_dt_c$delta - wide_reps_dt_c$delta %>% median()
      # wide_reps_dt_c <- wide_reps_dt_c %>% arrange((abs(dis))) %>% distinct(rep, .keep_all = TRUE) %>% select(-dis)

      # wide_reps_dt_s <- as.data.table(wide_reps_dt[-grep(wide_reps_dt$p, pattern = '^_'), ])
      # wide_reps_dt_s$dis <- wide_reps_dt_s$delta - wide_reps_dt_s$delta %>% median()
      # wide_reps_dt_s <- wide_reps_dt_s %>% arrange((abs(dis))) %>% distinct(rep, .keep_all = TRUE) %>% select(-dis)
      
      # red=4
      # p <- paste(s1, t, sep='_')
      # if(g %in% c("xxx")) {
      #   # wide_reps_dt_c <- wide_reps_dt_c[c(3,4,5), ]
      #   # red=1
      #   }else{
      #     wide_reps_dt_c <- drop_outliers_by_group(wide_reps_dt_c, group_col = "p", value_col = "delta", iterations = red)
      #   }
      
    
      
      # red=3
      # p <- paste(s1, t, sep='_')
      # # if(p %in% c("lox10_0h")) {red=4}
      # wide_reps_dt_s <- drop_outliers_by_group(wide_reps_dt_s, group_col = "p", value_col = "delta", iterations = red)
      
      
      # wide_reps_dt[, delta_raw := NULL]
      # wide_reps <- as.data.frame(rbind(wide_reps_dt_c, wide_reps_dt_s)) 
  

      # wide_reps <- wide_reps %>%
      #   mutate(
      #     species     = str_extract(p, pattern = '.*_'),
      #     delta       = gene - con,
      #     # control     = geom_mean_pos(delta[species == 'c_']),
      #     control     = mean(delta[species == '_'], na.rm = TRUE),
      #     ddelta      = delta - control,
      #     fold_change = 2^(-ddelta)
      #   ) %>% select(-species)
        
      # result <- wide_reps %>%
      #   group_by(p) %>%
      #   mutate(
      #     mean_fc = if_else(row_number() == 1, mean(fold_change, na.rm = TRUE), NA_real_)
      #   ) %>%
      #   ungroup() %>%
      #   as.data.frame()

      # result_list[[paste(s1, g, t, sep = "_")]] <- 
      #       result[, c("p", 'data_resoure', "rep", "con", "gene", "hr", "delta", "control", "ddelta", "fold_change", "mean_fc")]
    }
  }
}



result_raw %>% lapply(nrow) %>% unlist() %>% table()

all_time <- result_raw %>%  list2df()


meta <- rownames(all_time) %>% header_cleaning('_') %>% mutate(
  miR    = V2,
) %>% select(-c(V1, V2, V3))
all_time <- cbind(meta, all_time)
all_time$strain <- all_time$p %>% str_extract(pattern = '.*_') %>% str_remove(pattern = '_')
# all_time$strain <- all_time$p %>% str_extract(pattern = '.*_') %>% str_remove(pattern = '_')
head(all_time)
species <- all_time$strain %>% unique()
miR <- all_time$miR %>% unique()
time <- all_time$hr %>% unique()
# miR <- "TCP43" 

control <- all_time %>% filter(hr == 0)

mean_control <- TRUE
g='MIR159'
s1='B73'
t='1h'
result_time_list <- result_raw_time <- list()
for (s1 in species) {
  # s1 <- species[1]
  for (g in miR) {
    # g <- miR[1]

    for(t in time) {
      if(t == 0) next
      wide_reps_dt <- all_time %>% filter(strain == s1, miR %in% g, hr %in% c(0, t)) %>%
        mutate(p_filter = paste(p, data_resoure, rep, sep = '_')) %>%
        # distinct(p_filter, .keep_all = TRUE) %>%
        select(p, data_resoure, dup_resource, note, rep, con, gene, hr, )
    
      wide_reps_dt <- wide_reps_dt %>%
          mutate(
            species    = str_extract(p, pattern = '.*_'),
            delta       = gene - con
          ) %>% select(-species)
      
      result_raw_time[[paste(s1, g, t, 'hr', sep = "_")]]  <- wide_reps_dt %>% as.data.frame()

      wide_reps_dt <- wide_reps_dt %>% filter(con < 40 & gene < 40)  # remove Ct > 35
      if(length(unique(wide_reps_dt$p)) < 2) next

      wide_reps_dt_c <- as.data.table(wide_reps_dt[grep(wide_reps_dt$p, pattern = '0h'), ])
      # if(s1 == 's' & g == 'MIR319' & t == 1) {
      #   wide_reps_dt_c <- wide_reps_dt_c %>% top_n(5, delta)
      # }
      wide_reps_dt_c$dis <- wide_reps_dt_c$delta - wide_reps_dt_c$delta %>% median()
      wide_reps_dt_c <- wide_reps_dt_c %>% arrange((abs(dis))) #%>% distinct(rep, .keep_all = TRUE) %>% select(-dis)
      wide_reps_dt_c <- wide_reps_dt_c %>%
        mutate(
          center  = median(delta),
          mad     = mad(delta, constant = 1),
          z_score = ((delta - center) / mad),
          n       = sum(!is.na(delta))
        ) %>% arrange(abs(z_score)) %>% distinct(rep, .keep_all = TRUE) %>% select(-c(center, mad, z_score, n, dis))
      

      wide_reps_dt_s <- as.data.table(wide_reps_dt[-grep(wide_reps_dt$p, pattern = '0h'), ])
      # if(s1 == 's' & g == 'MIR319' & t == 1) {
      #   wide_reps_dt_s <- wide_reps_dt_s %>% top_n(-5, delta)
      # }
      wide_reps_dt_s$dis <- wide_reps_dt_s$delta - wide_reps_dt_s$delta %>% median()
      wide_reps_dt_s <- wide_reps_dt_s %>% arrange((abs(dis))) #%>% distinct(rep, .keep_all = TRUE) %>% select(-dis)
      wide_reps_dt_s <- wide_reps_dt_s %>%
        mutate(
          center  = median(delta),
          mad     = mad(delta, constant = 1),
          z_score = abs((delta - center) / mad),
          n       = sum(!is.na(delta))
        ) %>% arrange(abs(z_score)) %>% distinct(rep, .keep_all = TRUE) %>% select(-c(center, mad, z_score, n, dis))
      

      p <- paste(s1, g, sep='_')
      if(p %in% c('xxx')) {
        wide_reps_dt_c <- wide_reps_dt_c %>% filter(rep %in% c(1:3))
        }else{
          wide_reps_dt_c <- drop_outliers_by_group(wide_reps_dt_c, max_sample=3)
        }
      
    
      
      p <- paste(s1, g, sep='_')
      if(p %in% c('xxx')) {
        wide_reps_dt_s <- wide_reps_dt_s %>% filter(rep %in% c(1:3))
        }else{
          wide_reps_dt_s <- drop_outliers_by_group(wide_reps_dt_s, max_sample=3)
        }
      
      
      
      # wide_reps_dt[, delta_raw := NULL]
      wide_reps <- as.data.frame(rbind(wide_reps_dt_c, wide_reps_dt_s)) 


      wide_reps <- wide_reps %>%
        mutate(
          species     = str_extract(p, pattern = '.*_'),
          delta       = gene - con,
          # control     = geom_mean_pos(delta[species == 'c_']),
          pre_fold_change = 2^(-delta), 
          control     = mean(pre_fold_change[hr == '0'], na.rm = TRUE),
          fold_change = 2^(-delta) / control
        ) %>% select(-species)
      # wide_reps <- wide_reps %>%
      #   mutate(
      #     species     = str_extract(p, pattern = '.*_'),
      #     delta       = gene - con,
      #     # control     = geom_mean_pos(delta[species == 'c_']),
      #     control     = mean(delta[hr == 0], na.rm = TRUE),
      #     ddelta      = delta - control,
      #     fold_change = 2^(-ddelta)
      #   ) %>% select(-species)
      
      result <- wide_reps %>%
        group_by(p) %>%
        mutate(
          mean_fc = if_else(row_number() == 1, mean(fold_change, na.rm = TRUE), NA_real_)
        ) %>%
        ungroup() %>%
        as.data.frame()

      result <- replace(result, is.na(result), '')
      result_time_list[[paste(s1, g, t, 'hr', sep = "_")]] <- result
    }
    
  }
}

lapply(result_time_list, nrow)



library(openxlsx)
# wb <- createWorkbook()
# for (i in names(result_list)) {
#   addWorksheet(wb, sheetName = i)
#   writeData(wb, sheet = i, x = result_list[[i]])
# }

# saveWorkbook(wb, file = file.path(out_dir, "all_results.xlsx"), overwrite = TRUE)

write_workbook_from_list(result_raw,      file.path(out_dir, "all_raw_results.xlsx"))
write_workbook_from_list(result_raw_time, file.path(out_dir, "all_raw_time_results.xlsx"))
write_workbook_from_list(result_time_list, file.path(out_dir, "all_time_results.xlsx"))


sample_order <- c("B73", "lox2", "lox4", "lox5", "lox45", "lox6", "lox810", "lox10", "opr2", "opr78")
# plot_list <- make_plot_list(result_list, skip_small = TRUE)
# save_plot_list(plot_list, figure_dir)
# combine_plot <- make_combine_plot(plot_list, sample_order)
# save_combine_plot(combine_plot, figure_dir)

# end of file

species <- all_time$strain %>% unique()

for(spe in species) {
  
    sub_list <- result_time_list[grep(names(result_time_list), pattern= paste0('^', spe))]

    sub <- rbind((sub_list %>% list2df)[1:6, ], sub_list[-1] %>% lapply(function(x) {x %>% filter(hr != 0)}) %>% list2df)
    sub <- sub_list %>% list2df
    sub <- sub %>% mutate(gene_name=rownames(sub) %>% header_cleaning('_') %>% pull(V2) %>% str_remove('\\.[0-9]$'), 
                          pg=paste(p, gene_name, sep = '_'))

    plot_list <- list()

    for(g in unique(sub$gene_name)) {
      # g <- unique(sub$gene_name)[1]
      sub_plot_df <- sub %>% filter(gene_name == g) 
      plot_df <- sub_plot_df %>%
        group_by(pg) %>%
        summarise(
          mean_fc = mean(fold_change, na.rm = TRUE),
          sd_fc   = sd(fold_change, na.rm = TRUE),
          .groups = "drop"
        )

      plot_df$y_max     <- plot_df$mean_fc + plot_df$sd_fc
      plot_df$bracket_y <- plot_df$y_max * 1.05
      plot_df$stars_y   <- plot_df$y_max * 1.1


      ymax <- plot_df$y_max %>% max

      pv_vec <- c()
      for(pgs in c(plot_df$pg[-1])) {
        tt <- t.test(fold_change ~ pg, data = sub_plot_df %>% filter(pg %in% c(plot_df$pg[1], pgs)), var.equal = FALSE, alternative = "two.sided")
        pval <- tt$p.value
        pv_vec <- c(pv_vec, pval) 
      }

      names(pv_vec) <- plot_df$pg[-1]


      stars <- case_when(
        pv_vec <= 0.001 ~ "***",
        pv_vec <= 0.01  ~ "**",
        pv_vec <= 0.05  ~ "*",
        TRUE          ~ "ns"
      )
      names(stars) <- plot_df$pg[-1]

      plot_df$pg <- factor(plot_df$pg, levels = c(plot_df$pg))
      min_y <- min(plot_df$mean_fc - plot_df$sd_fc, na.rm = TRUE)
      min_y <- -0.1
      ymax = max(plot_df$mean_fc + plot_df$sd_fc)
      plot_list[[g]] <- 
        ggplot(plot_df, aes(x = pg, y = mean_fc)) +
              geom_col(width = 0.6, fill = "#d57500") +
              geom_errorbar(
                aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc),
                width = 0.15,
                linewidth = 0.8
              ) +
              # annotate("segment", x = 1, xend = 2, y = plot_df$bracket_y, yend = plot_df$bracket_y, linewidth = 0.7) +
              annotate("text", x = plot_df$pg[-1], y = plot_df$stars_y[-1]+0.05, label = stars, size = 12) +
              scale_y_continuous(limits = c(min_y, ymax * 1.2), breaks = scales::pretty_breaks(n = 5)) +
              geom_text(
                aes(label = round(mean_fc, 2), y=0),
                # vjust = -2.5,
                size = 6
              ) +
              labs(
                x = "Genotype",
                y = 'Fold Change',
                title = paste( toupper(spe), ', ', g, ", mean ± sd", sep = "")
              ) +
              theme_bw() + gg_theme +
              theme(axis.title.x = element_blank(), axis.text.x = element_text(size=16))
              # theme(axis.title.y = element_text(size=), axis.text.y = element_text(size=14))

      }


      ggsave(
          filename = file.path(figure_dir, paste0(spe, "_bar.png")),
          plot     = plot_list %>% cowplot::plot_grid(plotlist = ., ncol = 1)+ggtitle("L"),
          width    = 20,
          height   = 12,
          dpi      = 300
        )


}








# plot_list <- plot_df_list <- list()
# for (i in names(result_time_list)) {
#   # i <- names(result_time_list)[6]
#   # mir <- i %>% header_cleaning('_') %>% pull(V2)
#   # result_time_list[grep(mir, names(result_time_list))]
#   plot_df <- result_time_list[[i]] %>%
#     group_by(p) %>%
#     summarise(
#       mean_fc = mean(fold_change, na.rm = TRUE),
#       sd_fc   = sd(fold_change, na.rm = TRUE),
#       .groups = "drop"
#     )
#   if(any(is.na(plot_df$sd_fc))) next
#   tt <- t.test(fold_change ~ p, data = result_time_list[[i]], var.equal = FALSE, alternative = "two.sided")
#   pval <- tt$p.value

#   stars <- case_when(
#     pval <= 0.001 ~ "***",
#     pval <= 0.01  ~ "**",
#     pval <= 0.05  ~ "*",
#     TRUE          ~ "ns"
#   )

#   plot_df$y_max     <- max(plot_df$mean_fc + plot_df$sd_fc, na.rm = TRUE)
#   plot_df$bracket_y <- plot_df$y_max * 1.05
#   plot_df$stars_y   <- plot_df$y_max * 1.1

#   plot_df_list[[i]] <- plot_df
# }


# for(i in names(result_time_list)) {
#   # i <- names(result_time_list)[2]
#   mir <- i %>% header_cleaning('_') %>% pull(V2)
#   ymax <- plot_df_list[grep(mir, names(plot_df_list))] %>% lapply(FUN=function(x) {x[, 4] %>% max()}) %>% unlist() %>% max()
#   plot_df <- plot_df_list[[i]]
#   if(is.null(plot_df)) next
#     tt <- t.test(fold_change ~ p, data = result_time_list[[i]], var.equal = FALSE, alternative = "two.sided")
#   pval <- tt$p.value

#   stars <- case_when(
#     pval <= 0.001 ~ "***",
#     pval <= 0.01  ~ "**",
#     pval <= 0.05  ~ "*",
#     TRUE          ~ "ns"
#   )
  
#   min_y <- min(plot_df$mean_fc - plot_df$sd_fc, na.rm = TRUE)
#   if(min_y < -0.5) {min_y <- min_y} else {min_y <- -0.5}
#   ymax = max(plot_df$mean_fc + plot_df$sd_fc)
#   plot_list[[i]] <-
#       ggplot(plot_df, aes(x = p, y = mean_fc)) +
#       geom_col(width = 0.6, fill = "#d57500") +
#       geom_errorbar(
#         aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc),
#         width = 0.15,
#         linewidth = 0.8
#       ) +
#       annotate("segment", x = 1, xend = 2, y = plot_df$bracket_y[1], yend = plot_df$bracket_y[1], linewidth = 0.7) +
#       annotate("text", x = 1.5, y = plot_df$stars_y[1], label = stars, size = 12) +
#       scale_y_continuous(limits = c(min_y, ymax * 1.2), breaks = scales::pretty_breaks(n = 5)) +
#       labs(x = "Sample", y = "Fold Change", title = paste(i, "Fold Change", sep = ", ")) +
#       theme_bw() + gg_theme +
#       theme(axis.title.x = element_blank())
# }
  

# # for (i in names(plot_list)) {
# #   png(filename = file.path(figure_dir, paste0(i, "_result.png")), width = 1600, height = 1900, res = 300)
# #   print(plot_list[[i]])
# #   dev.off()
# # }

# miR <- names(plot_list) %>% header_cleaning('_') %>% pull(V2) %>% unique()
# combine_plot <- list()
# for (i in miR) {
#   plots <- plot_list[grep(i, names(plot_list))]
#   names(plots) <- names(plots) %>% header_cleaning('_') %>% pull(V1) %>% unique()
#   plots <- plots[c("B73", "lox2", "lox4", "lox5", "lox45", "lox6", "lox810", "lox10", "opr2", "opr78")]
#   plots <- plots[!vapply(plots, is.null, logical(1))]

#   combine_plot[[i]] <- cowplot::plot_grid(plotlist = plots, nrow = 2)
# }


# for (i in names(combine_plot)) {
#   png(file.path(figure_dir, paste0(i, "_combined.png")), width = 9000, height = 5000, res = 300)
#   print(combine_plot[[i]])
#   dev.off()
# }



# miR <- names(plot_list) %>% header_cleaning('_') %>% pull(V1) %>% unique()
# combine_plot <- list()
# for (i in miR) {
#   plots <- plot_list[grep(i, names(plot_list))]
#   names(plots) <- names(plots) %>% header_cleaning('_') %>% pull(V2) %>% unique()
#   # plots <- plots[c("B73", "lox2", "lox4", "lox5", "lox45", "lox6", "lox810", "lox10", "opr2", "opr78")]
#   # plots <- plots[!vapply(plots, is.null, logical(1))]

#   combine_plot[[i]] <- cowplot::plot_grid(plotlist = plots, nrow = 2)
# }


# for (i in names(combine_plot)) {
#   png(file.path(figure_dir, paste0(i, "_combined.png")), width = 9000, height = 5000, res = 300)
#   print(combine_plot[[i]])
#   dev.off()
# }



# # strain <- names(plot_list) %>% header_cleaning('_') %>% pull(V1) %>% unique()
# # combine_plot <- list()
# # for (i in strain) {
# #   plots <- plot_list[grep(i, names(plot_list))]
# #   names(plots) <- names(plots) %>% header_cleaning('_') %>% pull(V1) %>% unique()
# #   # plots <- plots[c("", "lox2", "lox4", "lox5", "lox45", "lox6", "lox810", "lox10", "opr2", "opr78")]
# #   # plots <- plots[!vapply(plots, is.null, logical(1))]

# #   combine_plot[[i]] <- cowplot::plot_grid(plotlist = plots, nrow = 2)
# # }

# # for (i in names(combine_plot)) {
# #   png(file.path(figure_dir, paste0(i, "_combined1227.png")), width = 9000, height = 5000, res = 300)
# #   print(combine_plot[[i]])
# #   dev.off()
# # }


