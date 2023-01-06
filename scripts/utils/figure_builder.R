###########################################################
# SCRIPT CONTAINING THE CODE TO GENERATE ALL MAIN FIGURES #
###########################################################

if (!require(tidyverse)) install.packages("tidyverse")
if (!require(magrittr)) install.packages("magrittr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggpubr)) install.packages("ggpubr")
if (!require(ggrepel)) install.packages("ggrepel")
if (!require(BiocManager)) install.packages("BiocManager")
if (!require(ComplexHeatmap)) BiocManager::install("ComplexHeatmap")
if (!require(EnhancedVolcano)) BiocManager::install("EnhancedVolcano")
if (!require(SummarizedExperiment)) BiocManager::install("SummarizedExperiment")
if (!require(ggtext)) install.packages("ggtext")
if (!require(ggrastr)) install.packages("ggrastr")
if (!require(gridExtra)) install.packages("gridExtra")
if (!require(kSamples)) install.packages("kSamples")
if (!require(reshape2)) install.packages("reshape2")
if (!require(ggbreak)) install.packages("ggbreak")
if (!require(rstatix)) install.packages("rstatix")
if (!require(ggvenn)) install.packages("ggvenn")
if (!require(scales)) install.packages("scales")
if (!require(ggpmisc)) install.packages("ggpmisc")
if (!require(cowplot)) install.packages("cowplot")
if (!require(DESeq2)) install.packages("DESeq2")
if (!require(gridtext)) remotes::install_github("wilkelab/gridtext")



# load plotting config
source("scripts/utils/config.R")

#
# This script contains only the code to create the figures and subfigures.
# It requires the main analysis files (.Rmd) to be run prior to its execution.
#


##################
##################
## MAIN FIGURES ##
##################
##################

############
# FIGURE 1 #
############
############

#####
# A #
#####


build_fig_1a_1 <- function () {
  
  comparisons <- list(
    c('GN', 'B'),
    c('GN', 'DC'),
    c('GN', 'Mono'),
    c('GN', 'NK'),
    c('GN', 'T')
  )
  
  human_n <- readRDS("data/processed/data_fig_2a_1.rds")
  
  human_n <- human_n %>%
    mutate(lineage_short = ifelse(
      lineage_clean == "B cells",
      "B",
      ifelse(lineage_clean == "Dendritic cells",
             "DC",
             ifelse(lineage_clean == "Monocytes",
                    "Mono",
                    ifelse (lineage_clean == "Neutrophils",
                            "GN",
                            ifelse (lineage_clean == "NK cells",
                                    "NK",
                                    "T")
                            )
                    )
             )
      )
    )

  
  
  plt <- ggplot(human_n, mapping = aes(x=lineage_short, y=frac_top, color = lineage_short)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.5) +
    stat_compare_means(method = "wilcox.test",
                       comparisons = comparisons,
                       exact = F,
                       label = 'p.signif',
                       size = 1.5,
                       bracket.size = 0.25,
                       vjust = 0.2
    ) +
    theme_rgb() +
    theme(axis.title = element_text()) +
    # theme_classic() +
    # theme(legend.position = "none",
    #       panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #       axis.line = element_blank()#,
    #       #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    # ) +
    scale_color_manual(values = lineage_pal_short) +
    ylim(0.45, 0.93) +
    xlab("Lineage") +
    ylab("Fraction of reads mappable")
  
  return (plt)
  
}

build_fig_1a_2 <- function () {
  
  comparisons <- list(
    c('GN', 'B'),
    c('GN', 'DC'),
    c('GN', 'Mono'),
    c('GN', 'NK'),
    c('GN', 'T')
  )
  
  mouse_n <- readRDS("data/processed/data_fig_2a_2.rds")
  
  mouse_n <- mouse_n %>%
    mutate(lineage_short = ifelse(
      lineage_clean == "B cells",
      "B",
      ifelse(lineage_clean == "Dendritic cells",
             "DC",
             ifelse(lineage_clean == "Monocytes",
                    "Mono",
                    ifelse (lineage_clean == "Neutrophils",
                            "GN",
                            ifelse (lineage_clean == "NK cells",
                                    "NK",
                                    "T")
                            )
                    )
             )
      )
    )

  plt <- ggplot(mouse_n, mapping = aes(x=lineage_short, y=frac_top, color = lineage_short)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.5) +
    stat_compare_means(method = "wilcox.test",
                       comparisons = comparisons,
                       exact = F,
                       label = 'p.signif',
                       size = 1.5,
                       bracket.size = 0.25,
                       vjust = 0.2
    ) +
    theme_rgb() +
    theme(axis.title = element_text()) +
    # theme_classic() +
    # theme(legend.position = "none",
    #       panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #       axis.line = element_blank()#,
    #       # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    # ) +
    scale_color_manual(values = lineage_pal_short) +
    ylim(0.45, 0.93) +
    xlab("Lineage") +
    ylab("")
  
  return (plt)
  
}

build_fig_1a_3 <- function(return_color_lgd = FALSE) {
  load("data/processed/data_fig_1b.rda")
  if(return_color_lgd == FALSE){
    ggplot (plot_data, mapping = aes(x = PC1, y = PC2, color = lineage, shape = species)) +
      geom_point (size = 1.5,
                  alpha = 0.8) +
      theme_rgb () +
      labs (x = paste0("PC1, ", pc1_var_expl, "% of variance explained"),
            y = paste0("PC2, ", pc2_var_expl, "% of variance explained"),
            color = "Lineage",
            shape = "Species") +
      theme(plot.subtitle = element_text(hjust = 0.5),
            legend.position = "bottom",
            axis.title.x = element_text(),
            axis.title.y = element_text(),
            aspect.ratio = 1) +
      scale_color_manual(values = lineage_pal,
                         breaks = names(lineage_pal), guide = "none") +
      scale_shape_manual(values = c(19, 17),
                         breaks = c("Homo sapiens", "Mus musculus"),
                         labels = c("Human", "Mouse"),
                         guide = "none") +
      ggtitle("")
  } else{
    ggplot (plot_data, mapping = aes(x = PC1, y = PC2, color = lineage, shape = species)) +
      geom_point (size = 2.5,
                  alpha = 0.8) +
      theme_rgb () +
      labs (x = paste0("PC1, ", pc1_var_expl, "% of variance explained"),
            y = paste0("PC2, ", pc2_var_expl, "% of variance explained"),
            color = "Lineage",
            shape = "Species") +
      theme(plot.subtitle = element_text(hjust = 0.5),
            legend.position = "bottom",
            legend.box = "vertical",
            legend.margin = margin(),
            legend.title = element_text(face = "bold", size = 8),
            legend.text = element_text(size = 6),
            legend.title.align = 0.5,
            axis.title.x = element_text(),
            axis.title.y = element_text()) +
      scale_color_manual(values = lineage_pal,
                         breaks = names(lineage_pal)) +
      scale_shape_manual(values = c(19, 17),
                         breaks = c("Homo sapiens", "Mus musculus"),
                         labels = c("Human", "Mouse")) +
      guides(
        color = guide_legend(ncol=2,
                             byrow = TRUE,
                             title.position = "top"),
        shape = guide_legend(override.aes = list(color = species_pal_long),
                             ncol=2,byrow = TRUE,
                             title.position = "top"))
   }

}


#####
# B #
#####

build_fig_1b_1 <- function () {
  
  ind_frac_top_genes_human_seq_df <- readRDS("data/processed/data_fig_2b_1.rds")
  
  plt <- ggplot (ind_frac_top_genes_human_seq_df,
                 mapping = aes(x = n_genes,
                               y = lineage_n_genes_mean,
                               ymax = lineage_n_genes_max,
                               ymin = lineage_n_genes_min,
                               fill = lineage_clean)) +
    geom_line(aes(color = lineage_clean)) +
    geom_ribbon(alpha = 0.1) +
    theme_rgb() +
    scale_color_manual(values = lineage_pal, guide = "none") +
    scale_fill_manual(values = lineage_pal) +
    ylim(0, 0.8) +
    #ggtitle("Human") +
    xlab("N genes") +
    ylab("Fraction of total reads represented") +
    labs(color = "Lineage", fill = "Lineage") +
    theme_rgb() +
    theme(axis.title = element_text(), legend.title = element_text(face = "bold", size = 8), legend.margin = margin(c(0,0,0,0))) +
    guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2)))
    # theme_classic() +
    # theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #       axis.line = element_blank(),
    #       legend.position = "none")
  
  return (plt)
}

build_fig_1b_2 <- function () {
  
  ind_frac_top_genes_mouse_seq_df <- readRDS("data/processed/data_fig_2b_2.rds")
  
  plt <- ggplot (ind_frac_top_genes_mouse_seq_df,
                 mapping = aes(x = n_genes,
                               y = lineage_n_genes_mean,
                               ymax = lineage_n_genes_max,
                               ymin = lineage_n_genes_min,
                               fill = lineage_clean)) +
    geom_line(aes(color = lineage_clean)) +
    geom_ribbon(alpha = 0.1) +
    theme_rgb() +
    scale_color_manual(values = lineage_pal) +
    scale_fill_manual(values = lineage_pal) +
    ylim(0, 0.8) +
    #ggtitle("Human") +
    xlab("N genes") +
    ylab("") +
    labs(color = "Lineage") +
    theme_rgb() +
    theme(axis.title = element_text())
    # theme_classic() +
    # theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #       axis.line = element_blank(),
    #       legend.position = "none")
  
  return (plt)
}

build_fig_1b_3 <- function () {
  
  ind_frac_top_genes_human_seq_df_rpkm <- readRDS("data/processed/data_fig_2b_3.rds")
  
  plt <- ggplot (ind_frac_top_genes_human_seq_df_rpkm,
                 mapping = aes(x = n_genes,
                               y = lineage_n_genes_mean,
                               ymax = lineage_n_genes_max,
                               ymin = lineage_n_genes_min,
                               fill = lineage_clean)) +
    geom_line(aes(color = lineage_clean)) +
    geom_ribbon(alpha = 0.1) +
    theme_rgb() +
    scale_color_manual(values = lineage_pal) +
    scale_fill_manual(values = lineage_pal) +
    #ggtitle("Human") +
    xlab("") +
    ylab("") +
    labs(color = "Lineage") +
    theme_rgb()
  
  return (plt)
}

build_fig_1b_4 <- function () {
  
  ind_frac_top_genes_mouse_seq_df_rpkm <- readRDS("data/processed/data_fig_2b_4.rds")
  
  plt <- ggplot (ind_frac_top_genes_mouse_seq_df_rpkm,
                 mapping = aes(x = n_genes,
                               y = lineage_n_genes_mean,
                               ymax = lineage_n_genes_max,
                               ymin = lineage_n_genes_min,
                               fill = lineage_clean)) +
    geom_line(aes(color = lineage_clean)) +
    geom_ribbon(alpha = 0.1) +
    theme_rgb() +
    scale_color_manual(values = lineage_pal) +
    scale_fill_manual(values = lineage_pal) +
    #ggtitle("Human") +
    xlab("") +
    ylab("") +
    labs(color = "Lineage") +
    theme_rgb()
  
  return (plt)
}

build_fig_1b_5 <- function () {
  
  ind_frac_top_genes_human_seq_df <- readRDS("data/processed/data_fig_2b_5.rds")
  
  plt <- ggplot (ind_frac_top_genes_human_seq_df,
                 mapping = aes(x = n_genes,
                               y = lineage_n_genes_mean,
                               ymax = lineage_n_genes_max,
                               ymin = lineage_n_genes_min,
                               fill = lineage_clean)) +
    geom_line(aes(color = lineage_clean)) +
    geom_ribbon(alpha = 0.1) +
    theme_rgb() +
    scale_color_manual(values = lineage_pal, guide = "none") +
    scale_fill_manual(values = lineage_pal) +
    ylim(0, 0.8) +
    #ggtitle("Human") +
    xlab("N genes") +
    ylab("Fraction of total reads represented") +
    labs(color = "Lineage", fill = "Lineage") +
    theme_rgb() +
    theme(axis.title = element_text(), legend.title = element_text(face = "bold", size = 8), legend.margin = margin(c(0,0,0,0))) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))
  # theme_classic() +
  # theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  #       axis.line = element_blank(),
  #       legend.position = "none")
  
  return (plt)
}

build_fig_1b_6 <- function () {
  
  ind_frac_top_genes_human_seq_df <- readRDS("data/processed/data_fig_2b_6.rds")
  
  plt <- ggplot (ind_frac_top_genes_human_seq_df,
                 mapping = aes(x = n_genes,
                               y = lineage_n_genes_mean,
                               ymax = lineage_n_genes_max,
                               ymin = lineage_n_genes_min,
                               fill = lineage_clean)) +
    geom_line(aes(color = lineage_clean)) +
    geom_ribbon(alpha = 0.1) +
    theme_rgb() +
    scale_color_manual(values = lineage_pal, guide = "none") +
    scale_fill_manual(values = lineage_pal) +
    ylim(0, 0.8) +
    #ggtitle("Human") +
    xlab("N genes") +
    ylab("Fraction of total reads represented") +
    labs(color = "Lineage", fill = "Lineage") +
    theme_rgb() +
    theme(axis.title = element_text(), legend.title = element_text(face = "bold", size = 8), legend.margin = margin(c(0,0,0,0))) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))
  # theme_classic() +
  # theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  #       axis.line = element_blank(),
  #       legend.position = "none")
  
  return (plt)
}

#####
# C #
#####

# Helper for 2C
get_ad_results <- function(fractions_df){
  lineage <- unique(fractions_df$lineage)
  ad.results <- data.frame(cbind(lineage), P=NA)  
  for(i in 1:length(lineage)){
    current.lineage <- lineage[i]
    for(j in 1:(length(lineage)-1)){
      current.lineage2 <- lineage[lineage != current.lineage][j]
      current.ad <- ad.test(fractions_df[fractions_df$lineage == current.lineage,]$frac_mappable,
                            fractions_df[fractions_df$lineage == current.lineage2,]$frac_mappable)$ad[[5]]
      if(is.na(ad.results[lineage == current.lineage, "P"])){
        ad.results[lineage == current.lineage, "P"] <- current.ad
      } else if(ad.results[lineage == current.lineage, "P"] < current.ad){
        ad.results[lineage == current.lineage, "P"] <- current.ad
      }
    }
  }
  ad.results$Ptext <- ifelse(ad.results$P > 0.05, "ns", 
                             ifelse(ad.results$P > 0.01 & ad.results$P <= 0.05, "*",
                                    ifelse(ad.results$P > 0.001 & ad.results$P <= 0.01, "**",
                                           ifelse(ad.results$P > 0.0001 & ad.results$P <= 0.001, "***",
                                                  "****"))))
  return(ad.results)
}

build_fig_1c_1 <- function () {
  
  avg_mappable_fracs_human_seq_df <- readRDS("data/processed/fig_2c_1.rds")

  ad.results <- get_ad_results(avg_mappable_fracs_human_seq_df)
  
  avg_mappable_fracs_human_seq_df <-  merge(avg_mappable_fracs_human_seq_df, ad.results[,c("lineage", "Ptext")], all.x=T)
  avg_mappable_fracs_human_seq_df[duplicated(avg_mappable_fracs_human_seq_df$lineage) == TRUE, "Ptext"] <- NA
  
  plt <- ggplot (avg_mappable_fracs_human_seq_df,
                 mapping = aes(x = n,
                               y = frac_mappable,
                               group = sample,
                               color = lineage)) +
    xlab("N genes compared") +
    ylab("Fraction of genes mappable") +
    geom_line() +
    facet_wrap(~ lineage) +
    geom_text(aes(label = Ptext), x = 160, y = 0.85, color = "black", size = 3) +
    theme_rgb() +
    theme(axis.title = element_text()) +
    scale_color_manual(values = lineage_pal)
    
  
  return (plt)
}

build_fig_1c_2 <- function () {
  
  avg_mappable_fracs_mouse_seq_df <- readRDS("data/processed/fig_2c_2.rds")
  
  plt <- ggplot (avg_mappable_fracs_mouse_seq_df,
                 mapping = aes(x = n,
                               y = frac_mappable,
                               group = sample,
                               color = lineage)) +
    xlab("N genes compared") +
    ylab("Fraction of genes mappable") +
    geom_line() +
    facet_wrap(~ lineage) +
    # geom_text(aes(label = Ptext), x = 150, y = 0.9, color = "black", size = 3) +
    #geom_text(aes(label = Ptext), x = 150, y = 0.9, color = "black", size = 3) +
    theme_rgb() +
    theme(axis.title = element_text()) +
    scale_color_manual(values = lineage_pal)
  
  return (plt)
}


build_fig_1c_3 <- function () {
  
  avg_mappable_fracs_mouse_seq_df <- readRDS("data/processed/fig_2c_3.rds")
  
  plt <- ggplot (avg_mappable_fracs_mouse_seq_df,
                 mapping = aes(x = n,
                               y = frac_mappable,
                               group = sample,
                               color = lineage)) +
    xlab("N genes compared") +
    ylab("Fraction of genes mappable") +
    geom_line() +
    facet_wrap(~ lineage) +
    # geom_text(aes(label = Ptext), x = 150, y = 0.9, color = "black", size = 3) +
    #geom_text(aes(label = Ptext), x = 150, y = 0.9, color = "black", size = 3) +
    theme_rgb() +
    theme(axis.title = element_text()) +
    scale_color_manual(values = lineage_pal)
  
  return (plt)
}

build_fig_1c_4 <- function () {
  
  avg_mappable_fracs_mouse_seq_df <- readRDS("data/processed/fig_2c_4.rds")
  
  plt <- ggplot (avg_mappable_fracs_mouse_seq_df,
                 mapping = aes(x = n,
                               y = frac_mappable,
                               group = sample,
                               color = lineage)) +
    xlab("N genes compared") +
    ylab("Fraction of genes mappable") +
    geom_line() +
    facet_wrap(~ lineage) +
    # geom_text(aes(label = Ptext), x = 150, y = 0.9, color = "black", size = 3) +
    #geom_text(aes(label = Ptext), x = 150, y = 0.9, color = "black", size = 3) +
    theme_rgb() +
    theme(axis.title = element_text()) +
    scale_color_manual(values = lineage_pal)
  
  return (plt)
}


#####
# D #
#####

build_fig_1d_1 <- function (return_lgd = FALSE) {
  
  plt_ranking <- readRDS("data/processed/data_1_fig_1c.rds")
  plt_wide <- readRDS("data/processed/data_2_fig_1c.rds")
  

  
  set_sizes <- plt_ranking %>%
    group_by(lineage_clean) %>%
    summarise(size = length(lineage_clean)) %>%
    pull(size) %>%
    as.character()
  
  col_anno <- HeatmapAnnotation(test = anno_block(gp = gpar(fill = lineage_pal,
                                                            col = "transparent"),
                                                  labels = c("B", "DC", "Mono", "GN", "NK", "T"),
                                                  labels_gp = gpar(col = "white",
                                                                   fontsize = 6),
                                                  height = unit(2.5, "mm")),
                                # Lineage = str_remove(colnames(plt_wide), "_.*$"),
                                col = list(
                                  Species = species_pal_long#,
                                  # Lineage = lineage_pal
                                ),
                                Species = str_remove(colnames(plt_wide), "^.*_"),
                                # Species = anno_block(gp = gpar(fill = rep(species_pal_long, 6),
                                #                                col = "transparent"),
                                #                      # labels = str_remove(colnames(plt_wide), "^.*_"),
                                #                      labels_gp = gpar(col = "white",
                                #                                       fontsize = 6),
                                #                      ),
                                border = F,
                                show_legend = c(FALSE, FALSE),
                                annotation_name_gp= gpar(fontsize = 6),
                                annotation_legend_param = list(),
                                simple_anno_size = unit(2.5, "mm"))
  
  # Define which gene to label (same genes as in scatter plot)
  plt_data <- readRDS("data/processed/data_fig_1e.rds")
  gene_annotation_data <- plt_data %>%
    mutate(
      expr_sum = `Homo sapiens` + `Mus musculus`,
      annotation_color = case_when(
             lineage_clean %in% c("B cells", "Monocytes", "NK cells") ~ "black",
             TRUE ~ "darkgrey"
        )) %>%
    group_by(lineage_clean) %>%
    arrange(desc(expr_sum)) %>%
    mutate(
      rank = 1:length(lineage_clean),
      label = ifelse(rank <= 8, symbol, "")
      ) %>%
    ungroup() %>%
    filter(label !="") %>%
    mutate(position = match(label, rownames(plt_wide)))
  
  stopifnot(all(!duplicated(gene_annotation_data$label)))
  
  
  row_anno_genes <- rowAnnotation(
    gene = anno_mark(at = gene_annotation_data$position,
                     labels = gene_annotation_data$label,
                     labels_gp = gpar(col = gene_annotation_data$annotation_color,
                                      fontsize = 5,
                                      fontface = "italic"), 
                     link_gp = gpar(col = gene_annotation_data$annotation_color,
                                    lwd = 0.2),
                     padding = 0.2,
                     side ="right")
  )
  
  row_anno <- rowAnnotation(
    `Marker gene set` = anno_block(gp = gpar(fill = lineage_pal,
                                             col = "transparent"),
                                   labels = set_sizes,
                                   labels_gp = gpar(col = "white",
                                                    fontsize = 6),
                                   width = unit(2.5, "mm")),
    # col = list(
    #   `Marker gene set` = lineage_pal
    # ),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(2.5, "mm"))
  
  
  # RdBu Brewer (imported from config)
  colfun <- generate_rd_bu_colfun(-5, 5)
  
  plt <- Heatmap(plt_wide,
                 top_annotation = col_anno,
                 left_annotation = row_anno,
                 right_annotation = row_anno_genes,
                 cluster_rows = F,
                 cluster_columns = F,
                 show_row_names = F,
                 show_column_names = F,
                 col = colfun,
                 column_split = str_remove(colnames(plt_wide), "_.*$"),
                 row_gap = unit(0, "mm"),
                 row_split = plt_ranking$lineage_clean,
                 # name = "Centered\nlog2 cpm",
                 heatmap_legend_param = list(title = expression(bold("Centered log"["2"](CPM)))),
                 column_title = "Average expression across lineage samples",
                 row_title = "Lineage specific gene sets",
                 column_title_gp = gpar(fontsize = 6),
                 row_title_gp = gpar(fontsize = 6),
                 show_heatmap_legend = return_lgd,
                 border = F)
  
  if (!return_lgd) {
    return (plt)
  } else {
    lgd <- Legend(col_fun = colfun,
                  title = expression(bold("Centered log"["2"](CPM))),
                  title_gp = gpar(fontsize = 8, fontface = "bold"),
                  title_position = "topcenter",
                  labels_gp = gpar(fontsize = 8),
                  grid_width = unit(3, "mm"),
                  legend_height = unit(20, "mm"),
                  border = T)
    return (lgd)
  }
  
}

build_fig_1d_2 <- function (return_lgd = FALSE) {
  
  plt_ranking <- readRDS("data/processed/data_1_fig_2d_2.rds")
  plt_wide <- readRDS("data/processed/data_2_fig_2d_2.rds")
  
  plt_ranking_ordered <- plt_ranking %>%
    mutate(lineage_clean = ifelse(
      lineage_clean == "Neutrophils",
      "A_Neutrophils",
      ifelse(
        lineage_clean == "T cells",
        "C_T cells",
        ifelse(
          lineage_clean == "Monocytes",
          "B_Monocytes",
          ifelse(
            lineage_clean == "B cells",
            "D_B cells",
            ifelse(
              lineage_clean == "NK cells",
              "E_NK cells",
              "F_Dendritic cells"
            )
          )
        )
      )
    ))
  
  plt_wide_ordered <- plt_wide
  colnames(plt_wide_ordered) <- paste0(c("A_",
                                         "A_",
                                         "B_",
                                         "B_",
                                         "C_",
                                         "C_",
                                         "D_",
                                         "D_",
                                         "E_",
                                         "E_",
                                         "F_",
                                         "F_"),
                                       colnames(plt_wide))
  lineage_pal_ordered <- lineage_pal
  names(lineage_pal_ordered) <- paste0(c("D_",
                                         "F_",
                                         "B_",
                                         "A_",
                                         "E_",
                                         "C_"),
                                       names(lineage_pal))
  
  set_sizes <- plt_ranking_ordered %>%
    group_by(lineage_clean) %>%
    summarise(size = length(lineage_clean)) %>%
    pull(size) %>%
    as.character()
  
  col_anno <- HeatmapAnnotation(Lineage = anno_block(gp = gpar(fill = lineage_pal_ordered[c(4,3,6,1,5,2)],
                                                               col = "transparent"),
                                                     labels = c("GN", "Mono", "T", "B", "NK", "DC"),
                                                     labels_gp = gpar(col = "white",
                                                                      fontsize = 6),
                                                     height = unit(2.5, "mm")),
                                # Lineage = str_remove(colnames(plt_wide), "_.*$"),
                                Species = str_remove(colnames(plt_wide_ordered), "^.*_"),
                                col = list(
                                  Species = species_pal_long#,
                                  # Lineage = lineage_pal
                                ),
                                show_legend = c(FALSE, FALSE),
                                annotation_name_gp= gpar(fontsize = 6),
                                simple_anno_size = unit(2.5, "mm"))
  
  row_anno <- rowAnnotation(`Marker gene set` = anno_block(gp = gpar(fill = lineage_pal_ordered[c(4,3,6,1,5,2)],
                                                                     col = "transparent"),
                                                           labels = set_sizes,
                                                           labels_gp = gpar(col = "white",
                                                                            fontsize = 6),
                                                           width = unit(2.5, "mm")),
                            # col = list(
                            #   `Marker gene set` = lineage_pal_ordered
                            # ),
                            show_annotation_name = FALSE,
                            show_legend = FALSE,
                            simple_anno_size = unit(2.5, "mm"))
  
  # RdBu Brewer (imported from config)
  colfun <- generate_rd_bu_colfun(-5, 5)
  
  plt <- Heatmap(plt_wide_ordered,
                 top_annotation = col_anno,
                 left_annotation = row_anno,
                 cluster_rows = F,
                 cluster_columns = F,
                 show_row_names = F,
                 show_column_names = F,
                 col = colfun,
                 column_split = str_remove(colnames(plt_wide_ordered), "_.*$"),
                 row_gap = unit(0, "mm"),
                 row_split = plt_ranking_ordered$lineage_clean,
                 name = "Centered\nlog2 cpm",
                 column_title = "Average expression in samples",
                 row_title = "Lineage specific gene sets",
                 column_title_gp = gpar(fontsize = 6),
                 row_title_gp = gpar(fontsize = 6),
                 heatmap_legend_param = list(),
                 show_heatmap_legend = return_lgd)
  
  
  
  if (!return_lgd) {
    return (plt)
  } else {
    lgd <- Legend(col_fun = colfun,
                  title = "Centered\nlog2 cpm",
                  title_gp = gpar(fontsize = 8, fontface = "bold"),
                  labels_gp = gpar(fontsize = 8),
                  grid_width = unit(3, "mm"),
                  legend_height = unit(20, "mm"))
    return (lgd)
  }
  
}


#####
# E #
#####

build_fig_1e <- function (return_lgd = FALSE) {
  
  cors_cross <- readRDS("data/processed/data_1_fig_1d.rds")
  annotation <- readRDS("data/processed/data_2_fig_1d.rds")
  
  anno_col <- HeatmapAnnotation(Lineage = annotation$lineage_clean,
                                col = list(
                                  Species = species_pal_long,
                                  Lineage = lineage_pal
                                ),
                                Species = annotation$species,
                                # Species = anno_simple(col = species_pal_long, annotation$species, pch = ifelse(annotation$species == "Homo sapiens", 21, 2),
                                #                       pt_size = unit(0.5, "mm"),
                                #                       pt_gp = gpar(col = "white", fill = "white"),
                                #                       simple_anno_size = unit(2.5, "mm")),
                                show_legend = c(FALSE, FALSE),
                                annotation_name_gp = gpar(fontsize = 6),
                                annotation_name_side = "left",
                                simple_anno_size = unit(2.5, "mm"),
                                border = T)
  
  if (any(annotation$sample_id != rownames(cors_cross))) {
    warning ("Mismatch!")
  }
  
  anno_row <- rowAnnotation(Lineage = annotation$lineage_clean,
                            col = list(
                              Species = species_pal_long,
                              Lineage = lineage_pal
                            ),
                            Species = annotation$species,
                            show_legend = c(FALSE, FALSE),
                            show_annotation_name = c(FALSE, FALSE),
                            annotation_name_gp = gpar(fontsize = 6),
                            simple_anno_size = unit(2.5, "mm"),
                            border = T)
  
  # blues into oranges Brewer
  colfun <- generate_blues_w_oranges_colfun(-1, 1)
  
  
  plt <- Heatmap(cors_cross,
                 top_annotation = anno_col,
                 left_annotation = anno_row,
                 col = colfun,
                 row_dend_width = unit(2, "mm"),
                 column_dend_height = unit(2, "mm"),
                 show_row_names = F,
                 show_column_names = F,
                 name = "Pearson",
                 row_title = "Samples",
                 column_title = "Samples",
                 row_title_gp = gpar(fontsize = 6),
                 column_title_gp = gpar(fontsize = 6),
                 show_heatmap_legend = return_lgd,
                 column_dend_gp = gpar(lwd = 0.5),
                 row_dend_gp = gpar(lwd = 0.5),
                 border = T
                 )
  
  # sorry for that
  if (!return_lgd) {
    return (plt)
  } else {
    lgd <- Legend(col_fun = colfun, 
                  title = expression(bold(paste("Pearson","'", "s"~bolditalic(r), sep=""))),
                  title_gp = gpar(fontsize = 8, fontface = "bold"),
                  title_position = "topcenter",
                  labels_gp = gpar(fontsize = 8),
                  grid_width = unit(3, "mm"),
                  legend_height = unit(20, "mm"),
                  border = T)
    return (lgd)
  }
  
}


#####
# F #
#####
build_fig_1f <- function () {
  
  plt_data <- readRDS("data/processed/data_fig_1e.rds")
  
  # define which labels to show (below defines top 10 most expressed genes
  
  plt_data %<>%
    mutate(expr_sum = `Homo sapiens` + `Mus musculus`) %>%
    group_by(lineage_clean) %>%
    arrange(desc(expr_sum)) %>%
    mutate(rank = 1:length(lineage_clean),
           label = ifelse(rank <= 10, symbol, ""))
  set.seed(42)
  plt <- ggplot (plt_data, aes (x = `Homo sapiens`, y = `Mus musculus`, label = label)) +
    geom_point(aes(color = lineage_clean),
               size = 0.5) + #color = "Darkred"

    theme_rgb()+
    facet_wrap(~lineage_clean) +
    #scale_color_brewer(palette = "Set1")
    geom_text_repel(max.overlaps = Inf,
                    force = 42,
                    seed = 2,
                    size = 2,
                    # box.padding = 0.35,
                    min.segment.length = 0,
                    segment.size = 0.1,
                    aes (color = NULL), 
                    nudge_y = -1,
                    max.iter = 1e5,
                    fontface = "italic") +
    stat_cor(
      cor.coef.name = "r",
      method = "pearson",
             label.x.npc = "left",
             label.y.npc = "top",
             aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~"),color = NULL),
             size = 2
             ) +
    xlab("Mean-centered log<sub>2</sub>(CPM), homo sapiens") +
    ylab("Mean-centered log<sub>2</sub>(CPM), mus musculus") +
    xlim(c(min(plt_data$`Homo sapiens`), 15)) +
    ylim(c(min(plt_data$`Mus musculus`), 15)) +
    coord_fixed() +
    scale_color_manual(values = lineage_pal) +
    theme(legend.position = "none",
          axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
    # theme_classic() +
    # theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #       axis.line = element_blank(),
    #       legend.position = "none")
  
  return (plt)
  
}



############
# FIGURE 2 #
############
############

#####
# A #
#####

build_fig_2a_1 <- function () {
  
  plot_data <- readRDS("data/processed/data_5_fig_3a.rds")
  
  plot_data %<>%
    ungroup() %>%
    mutate(p_masked = ifelse(wilcox_p <= 0.05, -log10(wilcox_p), NA),
           label = ifelse(
             hs_avg + mm_avg > (sort(hs_avg + mm_avg, decreasing = T)[10]) |
               hs_avg > (0.88 * mm_avg + 8) |
               hs_avg < (0.88 * mm_avg - 8),
             symbol,
             ""),
           )
  
  plt <- ggplot(data = plot_data, mapping = aes(x = hs_avg,
                                         y = mm_avg,
                                         label = label)) +
    rasterize(geom_point(shape = 16,
                         size = 0.05), dpi = 300) +
    geom_density_2d(alpha = 0.8) +
    coord_fixed() +
    theme_rgb() +
    # geom_label_repel(max.overlaps = Inf) +
    stat_cor(size = 2,
             method = "pearson",
             cor.coef.name = "r",
             aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~"),color = NULL)) +
    xlab("Average log<sub>2</sub>(TPM+1), Human") +
    ylab("Average log<sub>2</sub>(TPM+1), Mouse") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
    
  
  # OLD CODE
  
  # plot_data_highlight_conserved <- plot_data %>%
  #   mutate(sum = hs_avg + mm_avg) %>%
  #   slice_max(order_by = sum, n = 40)
  # plot_data_highlight_divergent <- plot_data %>%
  #   #mutate(diff = abs(abs(hs_avg) - abs(mm_avg))) %>%
  #   #slice_max(order_by = diff, n = 80)
  #   # y > x + 1
  #   # y < x - 1
  #   filter(hs_avg > (0.88 * mm_avg + 6) |
  #            hs_avg < (0.88 * mm_avg - 6))
  # 
  # plot_data_no_highlight <- plot_data %>%
  #   filter(!gene_symbol %in% plot_data_highlight_conserved$gene_symbol,
  #          !gene_symbol %in% plot_data_highlight_divergent$gene_symbol)
  
  # plt <- ggplot(data = plot_data_no_highlight, mapping = aes(x = hs_avg, y = mm_avg, label = gene_symbol)) +
  #   geom_point(size = 0.4) +
  #   geom_point(data = plot_data_highlight_conserved, color = "darkgreen") +
  #   geom_point(data = plot_data_highlight_divergent, color = "darkred") +
  #   coord_fixed() +
  #   theme_rgb() +
  #   geom_label_repel(data = plot_data_highlight_conserved,
  #                    max.overlaps = 35) +
  #   geom_label_repel(data = plot_data_highlight_divergent,
  #                    max.overlaps = 20) +
  #   stat_cor() +
  #   xlab("Average LogRPKM, Homo sapiens") +
  #   ylab("Average LogRPKM, Mus musculus")
  
  
  
  return (plt)
  
}

build_fig_2a_2 <- function () {
  
  plot_data <- readRDS("data/processed/data_5_fig_3a.rds")
  plt_ranking <- readRDS("data/processed/data_1_fig_1c.rds")
  chea_res_up_df_all <- readRDS("data/processed/data_fig_5a_1.rds")
  
  tf_label <- c("JUNB", 
                "JUND", 
                "SPI1", 
                "KLF2", 
                "ATF3", 
                "CEBPA", 
                "CEBPB", 
                "CEBPE")
    
  # plot the data
  
  plot_data_highlight <- plot_data %>%
    filter(symbol %in% unique(chea_res_up_df_all$TF)) %>%
    mutate(sum_avg = hs_avg+mm_avg) %>%
    arrange(desc(sum_avg)) %>%
    dplyr::select(-sum_avg)
  top5_tf <- plot_data_highlight$symbol[1:5]
  
  
  
  plot_data_label <- plot_data %>%
    filter(symbol %in% c(unique(tf_label), top5_tf))
  
  unique(plot_data_label[plot_data_label$symbol %in% tf_label & !(plot_data_label$symbol %in% top5_tf),"symbol"])
  
  plot_data_highlight %<>%
    filter(!symbol %in% plot_data_label$symbol)
  
  plot_data_no_highlight <- plot_data %>%
    filter(!symbol %in% plot_data_highlight$symbol)
  
  plt <- ggplot(data = plot_data_no_highlight, mapping = aes(x = hs_avg, y = mm_avg, label = symbol)) +
    rasterize(geom_point(shape = 16,
                         size = 0.05,
                         color = "grey70"), dpi = 300) +
    # geom_density_2d() +
    geom_point(data = plot_data_highlight,
               color = "darkgreen",
               size = 0.1) +
    geom_point(data = plot_data_label,
               color = "darkred",
               size = 0.1) +
    coord_fixed() +
    theme_rgb() +
    geom_label_repel(data = plot_data_label,
                     size = 1.5,
                     # force = 5,
                     nudge_y = -0.2,
                     alpha = 0.7,
                     label.padding = 0.075,
                     box.padding = 0.1,
                     label.size = 0.1,
                     segment.color = "darkred",
                     segment.linetype = "dashed",
                     segment.size = 0.35,
                     min.segment.length = 0,
                     max.overlaps = Inf,
                     seed = 42,
                     fontface = "italic") +
    stat_cor(data = plot_data_highlight, size = 2,
             method = "pearson",
             cor.coef.name = "r",
             aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~"),color = NULL)) +
    xlab("Average log<sub>2</sub>(TPM+1), Human") +
    ylab("Average log<sub>2</sub>(TPM+1), Mouse") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
  
  return(plt)
  
}


build_fig_2a_3 <- function () {
  
  plot_data <- readRDS("data/processed/data_5_fig_3a.rds")
  plt_ranking <- readRDS("data/processed/data_1_fig_1c.rds")
  
  neutrophil_genes <- plt_ranking %>%
    filter(lineage_clean == "Neutrophils") %>%
    pull(symbol) 
    # `[`(1:50)
  
  # plot the data
  plot_data_highlight <- plot_data %>%
    filter(symbol %in% neutrophil_genes[11:length(neutrophil_genes)])
  
  plot_data_no_highlight <- plot_data %>%
    filter(!symbol %in% plot_data_highlight$symbol)
  
  plot_data_label <- plot_data %>%
    filter(symbol %in% neutrophil_genes[1:10])
  
  plt <- ggplot(data = plot_data_no_highlight, mapping = aes(x = hs_avg, y = mm_avg, label = symbol)) +
    rasterize(geom_point(shape = 16,
                         size = 0.05,
                         color = "grey70"), dpi = 300) +
    geom_point(data = plot_data_highlight,
               color = "darkgreen",
               size = 0.1) +
    geom_point(data = plot_data_label,
               color = "darkred",
               size = 0.1) +
    coord_fixed() +
    theme_rgb() +
    geom_label_repel(data = plot_data_label,
                     size = 1.5,
                     # force = 5,
                     nudge_y = -0.2,
                     alpha = 0.7,
                     label.padding = 0.075,
                     box.padding = 0.1,
                     label.size = 0.1,
                     segment.color = "darkred",
                     segment.linetype = "dashed",
                     segment.size = 0.35,
                     min.segment.length = 0,
                     max.overlaps = Inf,
                     seed = 42,
                     fontface = "italic") +
    stat_cor(data = plot_data_highlight, size = 2,
             method = "pearson",
             cor.coef.name = "r",
             aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~"),color = NULL)) +
    xlab("Average log<sub>2</sub>(TPM+1), Human") +
    ylab("Average log<sub>2</sub>(TPM+1), Mouse") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
  
  return(plt)
  
}


#####
# B #
#####

build_fig_2b <- function () {
  
  plot_data <- readRDS("data/processed/data_5_fig_3a.rds")
  
  # TODO: Ccl6 - no high confidence orthology
  neutrophil_genes <- c("S100A8", 
                        # "S100A9", 
                        "ELANE", 
                        "MPO", 
                        # "CEACAM8", 
                        # "AZU1", 
                        # "LYZ", 
                        "IL1B", 
                        # "CD177", 
                        "MME", 
                        "PADI4", 
                        "CAMP")
  
  # plot the data
  plot_data_highlight <- plot_data %>%
    filter(symbol %in% neutrophil_genes)
  
  plot_data_no_highlight <- plot_data %>%
    filter(!symbol %in% plot_data_highlight$symbol)
  
  plt <- ggplot(data = plot_data_no_highlight, mapping = aes(x = hs_avg, y = mm_avg, label = symbol)) +
    rasterize(geom_point(shape = 16,
                         size = 0.05,
                         color = "grey70"), dpi = 300) +
    geom_point(data = plot_data_highlight,
               color = "darkred",
               size = 0.1) +
    coord_fixed() +
    theme_rgb() +
    geom_label_repel(data = plot_data_highlight,
                     size = 1.5,
                     # force = 5,
                     # nudge_x = -1,
                     # nudge_y = 2,
                     alpha = 0.7,
                     label.padding = 0.075,
                     label.size = 0.1,
                     segment.color = "darkred",
                     segment.linetype = "dashed",
                     segment.size = 0.35,
                     min.segment.length = 0,
                     max.overlaps = Inf,
                     seed = 42,
                     fontface = "italic") +
    # stat_cor(data = plot_data_highlight) +
    xlab("Average log<sub>2</sub>(TPM+1), Human") +
    ylab("Average log<sub>2</sub>(TPM+1), Mouse") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
  
  return(plt)
  
}


####################################
# old C, now integrated in heatmap #
####################################
# 
# build_fig_2c_1 <- function() {
#   
#   data_cons <- readRDS("data/processed/data_fig_3d.rds")
#   
#   stat.test <- data_cons %>%
#     group_by(symbol) %>%
#     wilcox_test(count ~ species) %>%
#     adjust_pvalue(method = "holm") %>%
#     add_significance() %>% add_xy_position(x = "species")
#   
#   bxplt <- ggplot(data_cons, aes(x = species, y = count)) +
#     geom_violin() +
#     geom_jitter(size = .3) +
#     ylim(NA, max(data_cons$count)+2) +
#     facet_wrap(~symbol,
#                strip.position = "bottom") +
#     # stat_compare_means(label = "p.format",
#     #                    label.x.npc = "left",
#     #                    label.y.npc = 0.95, size = 3) +
#     stat_pvalue_manual(stat.test) +
#     theme_rgb() +
#     theme(strip.placement = "outside",
#           panel.spacing = unit(0, 'mm'),
#           axis.title.x = element_markdown(),
#           axis.title.y = element_markdown()) +
#     labs(title = "High concordance") +
#     xlab("Species") +
#     ylab("Log<sub>2</sub>(TPM+1)")
#   
#   return(bxplt)
#   
#   
# }
# 
# build_fig_2c_2 <- function() {
#   
#   data_hi_hs <- readRDS("data/processed/data_fig_3e.rds")
#   
#   stat.test <- data_hi_hs %>%
#     group_by(symbol) %>%
#     wilcox_test(count ~ species) %>%
#     adjust_pvalue(method = "holm") %>%
#     add_significance() %>% add_xy_position(x = "species")
#   
#   bxplt <- ggplot(data_hi_hs, aes(x = species, y = count)) +
#     geom_violin(
#       # outliers are already shown in jitter plot, prevent overplotting
#       # outlier.alpha = 0
#     ) +
#     geom_jitter(size = .3) +
#     ylim(NA, max(data_hi_hs$count)+2) +
#     facet_wrap(~symbol,
#                strip.position = "bottom") +
#     # stat_compare_means(label = "p.format",
#     #                    label.x.npc = "left",
#     #                    label.y.npc = 0.95, size = 3) +
#     stat_pvalue_manual(stat.test) +
#     theme_rgb() +
#     theme(strip.placement = "outside",
#           panel.spacing = unit(0, 'mm'),
#           axis.title.x = element_markdown(),
#           axis.title.y = element_markdown()) +
#     labs(title = "Low concordance: High in human") +
#     xlab("Species") +
#     ylab("Log<sub>2</sub>(TPM+1)")
#   
#   return(bxplt)
#   
# }
# 
# build_fig_2c_3 <- function() {
#   
#   data_hi_mm <- readRDS("data/processed/data_fig_3f.rds")
#   
#   stat.test <- data_hi_mm %>%
#     group_by(symbol) %>%
#     wilcox_test(count ~ species) %>%
#     adjust_pvalue(method = "holm") %>%
#     add_significance() %>% add_xy_position(x = "species")
#   
#   bxplt <- ggplot(data_hi_mm, aes(x = species, y = count)) +
#     geom_violin(
#       # outliers are already shown in jitter plot, prevent overplotting
#       # outlier.alpha = 0
#     ) +
#     geom_jitter(size = .3) +
#     ylim(NA, max(data_hi_mm$count)+2) +
#     facet_wrap(~symbol,
#                strip.position = "bottom") +
#     # stat_compare_means(label = "p.format",
#     #                    label.x.npc = "left",
#     #                    label.y.npc = 0.95, size = 3) +
#     stat_pvalue_manual(stat.test) +
#     theme_rgb() +
#     theme(strip.placement = "outside",
#           panel.spacing = unit(0, 'mm'),
#           axis.title.x = element_markdown(),
#           axis.title.y = element_markdown()) +
#     labs(title = "Low concordance: High in mouse") +
#     xlab("Species") +
#     ylab("Log<sub>2</sub>(TPM+1)")
#   
#   return(bxplt)
#   
# }

#####
# C #
#####

build_fig_2c <- function(plt_genes_hi = c("CSF3R", "ACTB"),
                         plt_genes_hs = c("MME", "NLRP6"),
                         plt_genes_mm = c("CTSE", "NEURL3")) {
  load("data/processed/data_fig_3c.rda")
  cross_se <- readRDS("data/processed/cross_se.rds")
  enrichr_go_bp <- readRDS("data/processed/addition_fig3c.rds")
  lin_spec_genes <- readRDS("data/processed/data_1_fig_1c.rds") %>%
    filter(lineage_clean == "Neutrophils") %>%
    pull(symbol)
  
  
  plt <- expr_data_all %>%
    pivot_wider(id_cols = c("symbol"),
                values_from = "count",
                names_from = "sample_id") %>%
    column_to_rownames("symbol")
  plt <- plt[c(all_sets, set_nm_hs, set_nm_mm), ]
  
  colfun <- generate_rd_bu_colfun(min = quantile(plt, 0.01, na.rm = T), max = quantile(plt, 0.99, na.rm = T))
  
  # HEATMAP ANNOTATIONS
  # manual annotations -- old
  # annotate_genenames <- c(
  #   # both
  #   "CSF3R", "MCL1", "ACTB", "B2M", # no JUNB
  #   # human
  #   "LYZ", "ALPL", "DYSF", "S100A4", "SLC19A1",
  #   # mouse
  #   "NEURL3", "PI16", "TMEM72", "TMSB4Y", # no CD101
  #   # homan nm
  #   "CXCL8", "FCGR3A",
  #   # mouse nm
  #   "Ccl6")
  # manual annotations, discussed in manuscript, not in top 6
  annotate_genenames <- c(
    # both
    "MCL1", # no JUNB, B2M; ACTB, CSF3R in top genes anyways
    # human
    "NLRP6", "GSDME", # discussed in manuscript
    # mouse
     # no CD101
    # homan nm
    "CXCL8",
    # mouse nm
    "Ccl6")
  
  annotate_genes <- match(annotate_genenames, rownames(plt))
  
  # automatic annotations
  annotate_auto_positions <- c(1:6)
  # for each geneset
  annotate_auto <- c(
    annotate_auto_positions,
    # human
    annotate_auto_positions + length(set_hi),
    # mouse
    annotate_auto_positions + length(c(set_hi, set_hs)),
    # human nm
    annotate_auto_positions + length(c(set_hi, set_hs, set_mm)),
    # mouse nm
    annotate_auto_positions + length(c(set_hi, set_hs, set_mm, set_nm_hs))
  )
  
  # merge
  annotate_all <- c(annotate_genes, annotate_auto) %>%
    unique() %>%
    sort()
  
annotate_auto_gennames <- rownames(plt)[annotate_auto]

manual_not_in_auto <- annotate_genenames[!annotate_genenames %in% annotate_auto_gennames]
print(paste0("Annotated manually and not in auto annotation: ", manual_not_in_auto))

  # MANAGE HEATMAP SPLITTING
  # rows
  row_split <- c(rep("a", length(set_hi)),
                 rep("b", length(set_hs)),
                 rep("c", length(set_mm)),
                 rep("d", length(set_nm_hs)),
                 rep("e", length(set_nm_mm)))
  
  # columns
  col_split <- cross_se %>%
    colData() %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    filter(sample_id %in% colnames(plt)) %>%
    dplyr::select(sample_id, species)
  samples_hs <- col_split %>%
    filter(species == "Hs") %>%
    pull(sample_id)
  samples_mm <- col_split %>%
    filter(species == "Mm") %>%
    pull(sample_id)
  
  
  # BUILD LEFT ANNOTATIONS
  desc_mapping <- c(
    `High in both` = "a",
    `High in human` = "b",
    `High in mouse` = "c",
    `Non mappable human` = "d",
    `Non mappable mouse` = "e"
  )
  
  anno_go_data <- Reduce(function(df1, df2) rbind(df1, df2), enrichr_go_bp) %>%
    filter(adj_p_value <= 0.05, rank <= 5) %>%
    transmute(name = str_remove(term_name, " \\(.*$"),
              rank,
              adj_p_value,
              overlap = as.numeric(lapply(overlapping_genes, length)),
              # overlapping_genes,
              # query_description,
              group = desc_mapping[query_description])
  
  # define a function to plot the GO results
  go_bp_fun <- function(index, nm){
    select_group <- unique(nm)
    gt <- (
      ggplot(anno_go_data %>%
             filter(group == select_group),
           aes(x = 1,
               y = reorder(name, desc(rank)),
               size = overlap,
               color = adj_p_value)) +
      geom_point() +
      scale_size(range = c(0, 2.5),
                 limits = c(min(anno_go_data$overlap), max(anno_go_data$overlap))) +
      scale_color_continuous(limits = c(min(anno_go_data$adj_p_value), max(anno_go_data$adj_p_value))) +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 6),
            axis.title.y = element_blank(),
            panel.grid = element_blank(),
            #panel.border = element_rect(),
            aspect.ratio = 5,
            legend.position = "none")
      ) %>%
      ggplot_build() %>%
      ggplot_gtable()
    pushViewport(viewport())
    #grid.rect()
    vp_width <- sum(gt$widths) + unit(0.05, "npc")
    pushViewport(viewport(x = 1, y = 0, 
             width = vp_width, height = 1,
             just = c("right", "bottom")))
    #grid.draw(roundrectGrob())
    grid.draw(gt)
    popViewport(2)
  }
  
  # define the annotation itself
  # gene_group_go_annotation <- rowAnnotation(
  #   go = anno_link(
  #     align_to = row_split,
  #     which = "row",
  #     side = "left",
  #     panel_fun = go_bp_fun,
  #     size = unit(2, "cm"),
  #     gap = unit(0, "cm"),
  #     width = unit(10, "cm"),
  #     link_width = unit(0, "cm"),
  #     link_gp = gpar(lwd = unit(0, "cm"))
  #   ),
  #   group = anno_block(labels = c(paste0("High in\nboth (", length(set_hi), ")"),
  #                                 paste0("High in\nhuman (", length(set_hs), ")"),
  #                                 paste0("High in\nmouse (", length(set_mm), ")"),
  #                                 paste0("High without\nmouse ortho (", length(set_nm_hs), ")"),
  #                                 paste0("High without\nhuman ortho (", length(set_nm_mm), ")")),
  #                      labels_gp = gpar(fontsize = 8)),
  #   gene = anno_mark(at = annotate_all,
  #                    labels = rownames(plt)[annotate_all],
  #                    labels_gp = gpar(fontsize = 6), side ="left")
  # )
  
  gene_group_go_annotation <- rowAnnotation(
    gene = anno_mark(at = annotate_all,
                     labels = rownames(plt)[annotate_all],
                     labels_gp = gpar(fontsize = 6), side ="left"),
    group = anno_block(labels = c(paste0("High in\nboth (", length(set_hi), ")"),
                                  paste0("High in\nhuman (", length(set_hs), ")"),
                                  paste0("High in\nmouse (", length(set_mm), ")"),
                                  paste0("High without\nmouse ortho (", length(set_nm_hs), ")"),
                                  paste0("High without\nhuman ortho (", length(set_nm_mm), ")")),
                       labels_gp = gpar(fontsize = 8)),
    inlingenes = anno_simple(
      as.character(rownames(plt) %in% lin_spec_genes),
      col = c(`TRUE` = "black", `FALSE` = "white"),
      border = TRUE
    ),
    show_annotation_name = c(gene = FALSE, group = FALSE ,inlingenes = FALSE)
  )

  
  
  # BUILD RIGHT ANNOTATIONS
  
  # Boxplot integration
  # get the indices to annotate
  annotate_high_both <- match(plt_genes_hi, rownames(plt))
  annotate_high_hs <- match(plt_genes_hs, rownames(plt))
  annotate_high_mm <- match(plt_genes_mm, rownames(plt))
  
  # define a list containing groups of genes to annotate
  annotate_zoom <- list(
    both = annotate_high_both,
    hs = annotate_high_hs,
    mm = annotate_high_mm
  )
  
  # define a function drawing the plot annotation
  panel_fun <- function(index, nm) {
    stopifnot(nm %in% c("both", "hs", "mm"))
    
    genes <- rownames(plt)[index]
    data_plt <- expr_data_all %>%
      filter(symbol %in% genes)
    
    # if (nm == "both") {
    #   # gt <- ggplot_gtable(ggplot_build(fig_2c_1))
    #   gt <- grid.grabExpr(print(build_fig_2c_1()))
    # } else if (nm == "hs") {
    #   # gt <- ggplot_gtable(ggplot_build(fig_2c_2))
    #   gt <- grid.grabExpr(print(build_fig_2c_2()))
    # } else if (nm == "mm") {
    #   # gt <- ggplot_gtable(ggplot_build(fig_2c_3))
    #   gt <- grid.grabExpr(print(build_fig_2c_3()))
    # }
    fig_title <- ifelse(
      all(genes %in% plt_genes_hi),
      "High concordance",
      ifelse(
        all(genes %in% plt_genes_hs),
        "Low concordance: High in human",
        ifelse(all(genes %in% plt_genes_mm),
               "Low concordance: High in mouse",
               "!!NO GROUP!!"
        )
      )
    )
    
    gt <- (ggplot(data_plt, aes(x = species, y = count)) +
             geom_violin() +
             geom_jitter(size = .3) +
             ylim(NA, max(data_plt$count)+2) +
             facet_wrap(~symbol,
                        strip.position = "bottom") +
             stat_compare_means(label = "p.format",
                                label.x.npc = "left",
                                label.y.npc = 0.95, size = 3) +
             theme_rgb() +
             theme(strip.placement = "outside",
                   panel.spacing = unit(0, 'mm'),
                   axis.title.x = element_markdown(),
                   axis.title.y = element_markdown()) +
             labs(title = fig_title) +
             xlab("Species") +
             ylab("Log<sub>2</sub>(TPM+1)")) %>%
      ggplot_build() %>%
      ggplot_gtable()
    
    pushViewport(viewport())
    grid.rect()
    grid.draw(gt)
    #grid.text(paste(genes, collapse = ", "))
    popViewport()
  }
  
  # define the annotation
  anno_all <- rowAnnotation(
    genes = anno_zoom(align_to = annotate_zoom,
                      which = "row",
                      panel_fun = panel_fun,
                      size = unit(4, "cm"),
                      gap = unit(0.5, "cm"),
                      width = unit(6, "cm")))
  
  
  # BUILD HEATMAPS
  
  # set the correct heatmap widths
  # 5 inch == 127mm (however this does not seem to be correct)
  total_width <- unit(5.5, "inch")
  
  hs_anno_width <- width.AnnotationFunction(gene_group_go_annotation)
  mm_anno_width <- width.AnnotationFunction(anno_all)
  
  hm_width <- (total_width - hs_anno_width - mm_anno_width)/2
  
  set.seed(42)
  hm_hs <- Heatmap(plt[, samples_hs],
                   show_row_names = F,
                   column_order = sample(1:ncol(plt[,samples_hs])),
                   show_column_names = F,
                   cluster_rows = F,
                   cluster_columns = F,
                   row_split = row_split,
                   row_title = NULL,
                   column_title = paste("Human, n =", length(samples_hs)),
                   column_title_gp = gpar(fontsize = 8),
                   col = colfun,
                   #right_annotation = gene_annotation,
                   left_annotation = gene_group_go_annotation,
                   # heatmap_legend_param = list(title = expression(bold("Log"["2"](TPM+1)))),
                   heatmap_legend_param = list(title = expression(bold("Log"["2"]("TPM+1")))),
  # )
                   # heatmap_legend_param = list(title = gt_render("Log<sub>2</sub>(RPKM)")),
                   width = hm_width)
  
  set.seed(42)
  hm_mm <- Heatmap(plt[, samples_mm],
                   column_order = sample(1:ncol(plt[,samples_mm])),
                   show_row_names = F,
                   show_column_names = F,
                   cluster_rows = F,
                   cluster_columns = F,
                   row_split = row_split,
                   row_title = NULL,
                   column_title = paste("Mouse, n =", length(samples_mm)),
                   column_title_gp = gpar(fontsize = 8),
                   col = colfun,
                   right_annotation = anno_all,
                   # heatmap_legend_param = list(title = expression(bold("Log"["2"](RPKM)))),
                   show_heatmap_legend = F,
                   width = hm_width)
  
  return(hm_hs+hm_mm)
}



#############
# C Updated #
#############

build_fig_2c_updated <- function(
  #plt_genes_hi = c("S100A9 : S100a9", "S100A8 : S100a8"),
  plt_genes_hi = c("CSF3R : Csf3r", "CXCR2 : Cxcr2"),
  plt_genes_hs = c("C5AR1 : C5ar1", "CXCR1 : Cxcr1"),
  # plt_genes_mm = c("UBB : Ubb", "CHI3L1 : Chil1")
  plt_genes_mm = c("IL1B : Il1b", "RETNLB : Retnlg")
  ) {
  
  # GET THE DATA
  
  joined_plt_all <- readRDS("data/processed/data_fig_2c_updated.rds") %>%
    mutate(match_uniq_rn = str_replace(match_uniq_rn, "__", " : "))
  
  # matrix containing all numerical cols for both species (currently only used for coloring computations)
  joined_hm_mat <- joined_plt_all %>%
    dplyr::select(c(
      match_uniq_rn,
      starts_with("Hs_SRX"),
      starts_with("Mm_SRX")
    )) %>%
    column_to_rownames("match_uniq_rn") %>%
    as.matrix()
  
  joined_hm_mat_hs <- joined_plt_all %>%
    dplyr::select(c(
      match_uniq_rn,
      starts_with("Hs_SRX")
    )) %>%
    column_to_rownames("match_uniq_rn") %>%
    as.matrix()
  
  joined_hm_mat_mm <- joined_plt_all %>%
    dplyr::select(c(
      match_uniq_rn,
      starts_with("Mm_SRX")
    )) %>%
    column_to_rownames("match_uniq_rn") %>%
    as.matrix()
  
  # be extra sure, the order matches
  stopifnot(all(rownames(joined_hm_mat) == joined_plt_all$match_uniq_rn))
  stopifnot(all(rownames(joined_hm_mat_hs) == joined_plt_all$match_uniq_rn))
  stopifnot(all(rownames(joined_hm_mat_mm) == joined_plt_all$match_uniq_rn))
  

  # BUILD ANNOTAIONS
  
  # BUILD LEFT ANNOTATION
  # - gene name labels
  # - block anno for title, set size
  # - simple color anno for ortho type
  # - simple color anno for origin species
  
  
  # manual annotations, discussed in manuscript, not in top 6
  # CSF3R, CXCR2, NCF4, MCL1, SPI1, JUNB
  # FCGR3A, FCGR3B, C5AR1, MCXCR1
  # Mmp9, Camp, Il1b, Retnlg
  annotate_genenames_manual <- c(
    # both
    "ACTB : Actb",
    "CSF3R : Csf3r",
    "CXCR2 : Cxcr2",
    "NCF4 : Ncf4",
    "MCL1 : Mcl1",
    "SPI1 : Spi1",
    "JUNB : Junb",
    "SELL : Sell",
    "S100A11 : S100a11",
    "TYROBP : Tyrobp",
    # human
    "FCGR3A : Fcgr4",
    "FCGR3B : Fcgr4",
    "C5AR1 : C5ar1",
    "CXCR1 : Cxcr1",
    "FOS : Fos",
    "LITAF : Litaf",
    "TNFAIP2 : Tnfaip2",
    "CD63 : Cd63",
    "MMP25 : Mmp25",
    "LYN : Lyn",
    # mouse
    "MMP9 : Mmp9",
    "CAMP : Camp",
    "IL1B : Il1b",
    "RETNLB : Retnlg",
    "CLEC4D : Clec4d",
    "SAMHD1 : Samhd1",
    "MMP8 : Mmp8",
    # no CD101
    # homan nm
    "HLA-A : NA",
    "HLA-B : NA",
    "HLA-C : NA",
    "HLA-E : NA",
    "ICAM3 : NA",
    "PRR13 : NA",
    "S100P : NA",
    "FCGR2A : NA",
    "CXCL8 : NA",
    # mouse nm
    "NA : Ccl6",
    "NA : Cd177",
    "NA : Cd52",
    "NA : Cd33",
    "NA : Fcgr3"
    )
  
  # get indices
  annotate_genenames_manual_indices <- match(annotate_genenames_manual,
                                             joined_plt_all$match_uniq_rn)
  
  # TODO: define
  annotate_genenames_auto <- c()
  annotate_genenames_auto_indices <- c()
  
  # merge
  annotate_genenames_all_indices <- c(annotate_genenames_manual_indices,
                                      annotate_genenames_auto_indices) %>%
    unique() %>%
    sort()
  
  orig_pal <- c(species_pal, "#FFFFFF")
  names(orig_pal) <- c("human", "mouse", "both")
  orthotype_pal <- c(
    `ortholog_NA` = "#FBF7F4",
    `ortholog_one2one` = "#6C9A8B",
    `ortholog_one2many` = "#EED2CC",
    `ortholog_many2many` = "#E8998D"
    
  )
  
  print(paste0("% of one2many in cat 1-3: ", sum(joined_plt_all$orthotype == "ortholog_one2many")/ sum(joined_plt_all$orthotype %in% c("ortholog_one2one", "ortholog_one2many", "ortholog_many2many")) * 100))
  
  left_annotation <- rowAnnotation(
    gene = anno_mark(at = annotate_genenames_all_indices,
                     labels = joined_plt_all$match_uniq_rn[annotate_genenames_all_indices],
                     labels_gp = gpar(fontsize = 6, fontface = "italic"),
                     side ="left"),
    group = anno_block(labels = c(paste0("High in\nboth (", sum(joined_plt_all$group_split == "a_hi_both"), ")"),
                                  paste0("High in\nhuman (", sum(joined_plt_all$group_split == "b_hi_hs"), ")"),
                                  paste0("High in\nmouse (", sum(joined_plt_all$group_split == "c_hi_mm"), ")"),
                                  paste0("Only in\nhuman (", sum(joined_plt_all$group_split == "d_not_mappable_human"), ")"),
                                  paste0("Only in\nmouse (", sum(joined_plt_all$group_split == "d_not_mappable_mouse"), ")")),
                       labels_gp = gpar(fontsize = 8)),
    orthotype = joined_plt_all$orthotype,
    orig_group = joined_plt_all$orig_group,
    col = list(
      orthotype = orthotype_pal,
      orig_group = orig_pal
    ),
    border = T,
    show_annotation_name = c(gene = FALSE,
                             group = FALSE,
                             orthotype = FALSE,
                             orig_group = FALSE),
    annotation_legend_param = list(
      orthotype = list(
        title = "Type",
        border = "black",
        at = c("ortholog_one2one",
               "ortholog_one2many",
               "ortholog_many2many",
               "ortholog_NA"),
        labels = c("1 to 1",
                   "1 to many",
                   "many to many",
                   "NA"),
        title_gp = gpar(fontsize = 8,
                        fontface = "bold"),
        labels_gp = gpar(fontsize = 6)
        ),
      orig_group = list(
        title = "Species",
        border = "black",
        at = c("human",
               "mouse",
               "both"),
        labels = c("Human",
                   "Mouse",
                   "Both"),
        title_gp = gpar(fontsize = 8,
                        fontface = "bold"),
        labels_gp = gpar(fontsize = 6)
        )
      )
    )
  
  
  
  # BUILD RIGHT ANNOTATIONS
  # violinplots for selected genes
  
  # violinplots integration
  # get the indices to annotate
  plot_indices_high_both <- match(plt_genes_hi,
                                  joined_plt_all$match_uniq_rn)
  plot_indices_high_hs <- match(plt_genes_hs,
                                joined_plt_all$match_uniq_rn)
  plot_indices_high_mm <- match(plt_genes_mm,
                                joined_plt_all$match_uniq_rn)
  
  # define a list containing groups of genes to annotate
  plot_indices_zoom <- list(
    both = plot_indices_high_both,
    hs = plot_indices_high_hs,
    mm = plot_indices_high_mm
  )
  
  
  # define a function drawing the plot annotation
  panel_fun <- function(index, nm) {
    stopifnot(nm %in% c("both", "hs", "mm"))
    
    genes <- joined_plt_all$match_uniq_rn[index]
    
    data <- joined_plt_all %>%
      filter(match_uniq_rn %in% genes) %>%
      dplyr::select(c(match_uniq_rn,
                      model_padj,
                      starts_with("Hs_SRX"),
                      starts_with("Mm_SRX"))) %>%
      pivot_longer(cols = -c(match_uniq_rn, model_padj),
                   names_to = "sample",
                   values_to = "count") %>%
      mutate(
        species = case_when(
          str_detect(sample, "^Hs_")~"Human",
          str_detect(sample, "^Mm_")~"Mouse",
          TRUE ~ "PROBLEM"
          )
        )
    
    # if (nm == "both") {
    #   # gt <- ggplot_gtable(ggplot_build(fig_2c_1))
    #   gt <- grid.grabExpr(print(build_fig_2c_1()))
    # } else if (nm == "hs") {
    #   # gt <- ggplot_gtable(ggplot_build(fig_2c_2))
    #   gt <- grid.grabExpr(print(build_fig_2c_2()))
    # } else if (nm == "mm") {
    #   # gt <- ggplot_gtable(ggplot_build(fig_2c_3))
    #   gt <- grid.grabExpr(print(build_fig_2c_3()))
    # }
    fig_title <- ifelse(
      all(genes %in% plt_genes_hi),
      "High concordance",
      ifelse(
        all(genes %in% plt_genes_hs),
        "Low concordance: High in human",
        ifelse(all(genes %in% plt_genes_mm),
               "Low concordance: High in mouse",
               "!!NO GROUP!!"
        )
      )
    )
    
    data %<>%
      # get. the. label. to. print. nicely.
      mutate(label = sapply(format(model_padj,
                                   digits = 3,
                                   nsmall = 2),
                            function(p) bquote(italic(P)[adjusted]==.(p))),
             label = replace(
               label,
               which(duplicated(label)),
               "")#,
             # match_uniq_rn = str_replace(match_uniq_rn, "__", "  ")
             )
    
    gt <- (ggplot(data, aes(x = species, y = count)) +
             geom_violin(scale = "width",
                         trim = T,
                         fill = "lightblue",
                         draw_quantiles = c(0.5)) +
             geom_jitter(size = .5,
                         shape = 21,
                         fill = "white",
                         stroke = 0.25) +
             ylim(0, max(data$count)+2) +
             facet_wrap(~match_uniq_rn,
                        strip.position = "bottom") +
             # stat_compare_means(label = "p.format",
             #                    label.x.npc = "left",
             #                    method = "wilcox",
             #                    label.y.npc = 0.95, size = 3) +
             geom_text(aes(x = 1.5,
                           y = max(data$count)+2,
                           label = label),
                       hjust = 0.5,
                       vjust = 1,
                       size = 2,
                       parse = T#,
                       # lwd = 0.1
                       ) +
             theme_rgb() +
             theme(strip.placement = "outside",
                   panel.spacing = unit(0, 'mm'),
                   strip.text = element_text(face = "italic"),
                   axis.title.x = element_markdown(),
                   axis.title.y = element_markdown()) +
             labs(title = fig_title) +
             xlab("Species") +
             ylab("Log<sub>2</sub>(TPM+1)")) %>%
      ggplot_build() %>%
      ggplot_gtable()
    
    pushViewport(viewport())
    grid.rect()
    grid.draw(gt)
    #grid.text(paste(genes, collapse = ", "))
    popViewport()
  }
  
  
  # define the annotation based on the function definition and gene info
  right_annotation <- rowAnnotation(
    genes = anno_zoom(align_to = plot_indices_zoom,
                      which = "row",
                      panel_fun = panel_fun,
                      size = unit(4, "cm"),
                      gap = unit(0.5, "cm"),
                      width = unit(6, "cm")
                      )
    )
  
  
  # BUILD HEATMAPS
  # define colors, remove outliers from coloring scale
  colfun <- generate_rd_bu_colfun(min = quantile(joined_hm_mat, 0.01, na.rm = T),
                                  max = quantile(joined_hm_mat, 0.99, na.rm = T))
  
  # set the correct heatmap widths
  # 5 inch == 127mm (however this does not seem to be correct)
  total_width <- unit(5.5, "inch")
  
  anno_width <- width.AnnotationFunction(left_annotation) + width.AnnotationFunction(right_annotation)
  
  hm_width <- convertUnit((total_width - anno_width)/2, "inch")
  
  set.seed(42)
  hm_hs <- Heatmap(joined_hm_mat_hs,
                   show_row_names = F,
                   column_order = sample(1:ncol(joined_hm_mat_hs)),
                   show_column_names = F,
                   cluster_rows = F,
                   cluster_columns = F,
                   row_split = joined_plt_all$group_split,
                   row_title = NULL,
                   column_title = paste("Human, n =", ncol(joined_hm_mat_hs)),
                   column_title_gp = gpar(fontsize = 8),
                   col = colfun,
                   left_annotation = left_annotation,
                   # heatmap_legend_param = list(title = expression(bold("Log"["2"](TPM+1)))),
                   heatmap_legend_param = list(
                     title = expression(bold("Log"["2"]("TPM+1"))),
                     border = T,
                     title_gp = gpar(fontsize = 8,
                                     fontface = "bold"),
                     labels_gp = gpar(fontsize = 6)
                     ),
                   width = hm_width,
                   border = T)
  
  set.seed(42)
  hm_mm <- Heatmap(joined_hm_mat_mm,
                   column_order = sample(1:ncol(joined_hm_mat_mm)),
                   show_row_names = F,
                   # show_row_names = T,
                   # row_names_side = "right",
                   # row_names_gp = gpar(fontsize = 2),
                   show_column_names = F,
                   cluster_rows = F,
                   cluster_columns = F,
                   row_split = joined_plt_all$group_split,
                   row_title = NULL,
                   column_title = paste("Mouse, n =", ncol(joined_hm_mat_mm)),
                   column_title_gp = gpar(fontsize = 8),
                   col = colfun,
                   right_annotation = right_annotation,
                   # heatmap_legend_param = list(title = expression(bold("Log"["2"](RPKM)))),
                   show_heatmap_legend = F,
                   width = hm_width,
                   border = T)
  
  return(hm_hs+hm_mm)
}



############
# FIGURE 3 #
############
############

#####
# A #
#####

build_fig_3a_1 <- function() {
  
  study_stats <- readRDS("data/processed/data_fig_4a.rds")
  
  return(tableGrob(study_stats))
  
  
}

build_fig_3a_2 <- function () {
  
  study_stats <- readRDS("data/processed/data_fig_4a.rds")
  
  plt_stats <- study_stats %>%
    arrange(species, tissues) %>%
    column_to_rownames("study_id") %>%
    dplyr::transmute(species = c("Hs" = "Human", "Mm" = "Mouse")[species],
                     model = model,
                     tissue = tissues,
                     stimuli = stimuli %>%
                       str_remove_all("_in_vitro|_in_vivo|_MTX_treatment|, JIA_untreated") %>%
                       str_replace_all("Y.pestis", "Y. pestis") %>%
                       str_replace_all("_", "-"),
                     resting = n_samp_resting,
                     stimulated = n_samp_stim
    )
  # Rename F.alocis stimulus to simplify plot
  plt_stats[rownames(plt_stats) == "Miralda I 2020","stimuli"] <- "F. alocis"
  colnames(plt_stats) <- str_to_title(colnames(plt_stats))
  
  plt_hm <- matrix(
    expand_grid(rows = 1:nrow(plt_stats), cols = 1:ncol(plt_stats)) %>%
      transmute(var = paste(rows, cols, sep = "_")) %>%
      pull(var),
    byrow = T,
    ncol = ncol(plt_stats)
  )
  rownames(plt_hm) <- rownames(plt_stats)
  colnames(plt_hm) <- colnames(plt_stats)
  
  
  
  samp_max <- max(c(plt_stats$Resting, plt_stats$Stimulated))
  samp_min <- min(c(plt_stats$Resting, plt_stats$Stimulated))
  
  rescale_hm <- function(x, x_min, x_max, a=0.5, b=1) {
    
    n = a + ((x-x_min))*(b-a)/(x_max-x_min)
    
    return (n)
  }
  
  colfun <- generate_blues_colfun(0, samp_max)
  
  plt <- Heatmap(
    plt_hm,
    # remove generated cells
    rect_gp = gpar(type = "none"),
    cluster_rows = F,
    cluster_columns = F,
    show_heatmap_legend = F,
    column_names_rot = 0,
    column_names_side = "top",
    column_names_centered = T,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x = x,
                y = y,
                width = width,
                height = height,
                gp = gpar(col = NA, fill = NA))
      if (colnames(plt_hm)[j] == "Species") {
        if(plt_stats[i, j] == "Human") fill = species_pal[["Hs"]] else fill = species_pal[["Mm"]]
        grid.rect(x = x,
                  y = y,
                  height = 0.9*height,
                  width = 0.9*height,
                  gp = gpar(fill = fill, col = NA, alpha = 0.5))
      }
      if (colnames(plt_hm)[j] == "Tissue") {
        if(plt_stats[i, j] == "PB") fill = "darkred"
        else if(plt_stats[i, j] == "BM") fill = "dodgerblue"
        else fill = "darkgreen"
        grid.rect(x = x,
                  y = y,
                  height = 0.9*height,
                  width = 1.7*height,
                  gp = gpar(fill = fill, col = NA, alpha = 0.5))
      }
      if (colnames(plt_hm)[j] %in% c("Resting", "Stimulated")) {
        # rescale values to n on interval [0.5, 1]
        n = rescale_hm(as.numeric(plt_stats[i, j]), samp_min, samp_max)
        grid.circle(x = x,
                    y = y,
                    r = n*0.5*min(unit.c(width, height)),
                    gp = gpar(fill = colfun(plt_stats[i, j]), col = NA, alpha = 0.5))
      }
      
      grid.text(plt_stats[i, j], x = x, y = y, gp = gpar(fontsize = 7))
    }
    
  )
  return(plt)
}

#####
# B #
#####
build_fig_3b <- function(return_lgd = FALSE) {
  
  volcano_stats <- readRDS("data/processed/data_fig_4b.rds")
  
  validated_genes <- c("CD40", "CD274", "CD14", "CD69", "IL4R")
  
  volcano_stats %<>%
    arrange(desc(lfc)) %>%
    mutate(rank = 1:length(symbol),
           label = ifelse(symbol %in% validated_genes | (rank <= 10 | rank >= max(rank)-9) & pval <= 0.05 & abs(lfc) >= 1 , symbol, NA))
  
  
  lab_italics <- ifelse(is.na(volcano_stats$label), NA, paste0( "italic('", volcano_stats$label, "')"))
  if(return_lgd == FALSE){
    set.seed(42)
    plt <- enhancedVolcanoFAR(
      volcano_stats,
      x = "lfc",
      y = "pval",
      lab = lab_italics,
      parseLabels = TRUE,
      drawConnectors = TRUE,
      max.overlaps = Inf,
      maxoverlapsConnectors = Inf,
      pointSize = 0.3,
      labSize = 2,
      pCutoff = 1.4e-39,
      FCcutoff = 0.5,
      FCcutoff_label = 0,
      title = "",
      titleLabSize = 0,
      legendLabSize = 7,
      axisLabSize = 5,
      subtitle = "",
      subtitleLabSize = 0,
      caption = "",
      captionLabSize = 0,
      legendPosition = "none"
    ) + theme_classic() + theme(axis.line = element_line(colour = 'black', size = 0.3),
              axis.ticks = element_line(colour = "black", size = 0.3),
              legend.position = "none",
              axis.title = element_text(size = 8)
              )
    return(plt)
  } else{
    plt <-EnhancedVolcano(
      volcano_stats,
      x = "lfc",
      y = "pval",
      lab = volcano_stats$symbol,
      pointSize = 0.7,
      labSize = 2,
      title = "",
      titleLabSize = 0,
      legendLabSize = 6,
      axisLabSize = 5,
      subtitle = "",
      subtitleLabSize = 0,
      caption = "",
      captionLabSize = 0,
      legendPosition = "bottom", 
    ) + 
      theme(legend.spacing.x = unit(0, "cm"), legend.background=element_rect(fill = alpha("white", 0)),
              legend.key=element_rect(fill = alpha("white", 0), color = alpha("white", 0))) +
      guides(color = guide_legend(override.aes = list(size = 2, fill = NA)))
    
    
    
    
    return(get_legend(plt))
    
  }
  
}

#####
# C #
#####

build_fig_3c_1 <- function(return_labels = FALSE) {
  
  heatmap_matrix <- readRDS("data/processed/data_fig_4c_1.rds")
  fisher_up_genes <- readRDS("data/processed/fisher_up_genes.rds")
  fisher_dn_genes <- readRDS("data/processed/fisher_dn_genes.rds")
  study_stats <- readRDS("data/processed/data_fig_4a.rds")
  
  validated_genes <- c(
    "CD14",
    "CD40",
    "CD69",
    "CD274",
    "IL4R"
    )
  label_genes_manual <- c(
    "TLR5",
    "MAP3K15"
  )
  
  merge_anno <- rownames(heatmap_matrix) %>%
    str_split(., pattern = "_") %>%
    sapply(., "[[",1) %>%
    data.frame() %>%
    dplyr::rename("study_id" = ".") %>%
    left_join(study_stats[,c("study_id", "model", "species")]) %>%
    mutate(orig = rownames(heatmap_matrix))
  
  #Confirm correct order
  stopifnot(identical(merge_anno$study_id, sapply(str_split(rownames(heatmap_matrix), pattern = "_"), "[[",1)))
  
  annotation_row <- rowAnnotation(Species = str_remove(rownames(heatmap_matrix), "^.*_"),
                                  Model = merge_anno$model,
                                  col = list (Species = species_pal,
                                              Model = model_pal),
                                  annotation_name_gp = gpar(fontsize = 6),
                                  annotation_legend_param = list(
                                    Species = list(
                                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                                      grid_width = unit(0.3, "cm"),
                                      labels = c("Human", "Mouse"),
                                      at = c("Hs", "Mm"),
                                      labels_gp = gpar(fontsize = 6),
                                      border = T
                                    ),
                                    Model = list(
                                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                                      grid_width = unit(0.3, "cm"),
                                      labels_gp = gpar(fontsize = 6),
                                      border = T
                                    )
                                  ),
                                  border = T)
  rownames(heatmap_matrix) <- str_remove(rownames(heatmap_matrix), "_(Hs|Mm)$")
  rownames(heatmap_matrix) <- str_remove(rownames(heatmap_matrix), "condition_")
  rownames(heatmap_matrix) <- str_replace_all(rownames(heatmap_matrix), "_", " ")
  
  col_split_df <- tibble(
    gene = colnames(heatmap_matrix),
    group = case_when(
      gene %in% fisher_up_genes ~ "up",
      gene %in% fisher_dn_genes ~ "dn"
    )
  )
  
  l <- length(colnames(heatmap_matrix))
  label_genes <- unique(
    c(
      colnames(heatmap_matrix)[c(1:15, l-(0:14))],
      validated_genes,
      label_genes_manual
      )
    )
  if(return_labels == TRUE){
    return(label_genes)
  }else{
  
  gene_annotation <- columnAnnotation(genes = anno_mark(at = which(colnames(heatmap_matrix) %in% label_genes),
                                                        labels = colnames(heatmap_matrix)[colnames(heatmap_matrix) %in% label_genes],
                                                        side = "bottom",
                                                        labels_gp = gpar(fontsize = 6, fontface = "italic"),
                                                        link_gp = gpar(lwd = 0.5)))
  
  split_annotation <- columnAnnotation(
    group = anno_block(
      labels = c(paste0("Downregulated \n (N = ", sum(col_split_df == "dn"), ")"), 
                 paste0("Upregulated \n (N = ", sum(col_split_df == "up"), ")")),
      labels_gp = gpar(fontsize = 6),
      height = unit(6, "mm")
    )
  )
  
  
  # RdBu Brewer
  colfun <- generate_rd_bu_colfun(-5, 5)
  
  plt <- Heatmap(heatmap_matrix,
                 heatmap_legend_param = list(title = expression(bold("Log"["2"]("FC"))),
                                             legend_height = unit(1.19, "cm"),
                                             grid_width = unit(0.25, "cm"),
                                             title_gp = gpar(fontsize = 8),
                                             labels_gp = gpar(fontsize = 6),
                                             border = T
                                             ), 
                 left_annotation = annotation_row,
                 bottom_annotation = gene_annotation,
                 show_column_names = FALSE,
                 cluster_rows = TRUE,
                 cluster_columns = FALSE,
                 use_raster = FALSE,
                 col = colfun,
                 column_split = col_split_df$group,
                 top_annotation = split_annotation,
                 column_title = c(),
                 row_names_max_width = unit(100, "cm"),
                 row_names_gp = gpar(fontsize = 6),
                 row_dend_width = unit(4, "mm"),
                 border = T)
  return (plt)
  }
  
}

build_fig_3c_2 <- function() {
  
  heatmap_matrix_scaled <- readRDS("data/processed/data_fig_4c_2.rds")
  
  
  annotation_row <- rowAnnotation(Species = str_remove(rownames(heatmap_matrix_scaled), "^.*_"),
                                  col = list (Species = species_pal))
  rownames(heatmap_matrix_scaled) <- str_remove(rownames(heatmap_matrix_scaled), "_(Hs|Mm)$")
  rownames(heatmap_matrix_scaled) <- str_remove(rownames(heatmap_matrix_scaled), "condition_")
  rownames(heatmap_matrix_scaled) <- str_replace_all(rownames(heatmap_matrix_scaled), "_", " ")
  
  #label_genes <- c("IL1A", "IL1B", "NFKBIZ", "JUNB", "JAK3")
  l <- length(colnames(heatmap_matrix_scaled))
  label_genes <- colnames(heatmap_matrix_scaled)[c(1:5, 10, 15, 20, 25, 30,
                                                   l-(0:4), l-9, l-14, l-19, l-24, l-29)]
  
  gene_annotation <- columnAnnotation(genes = anno_mark(at = which(colnames(heatmap_matrix_scaled) %in% label_genes),
                                                        labels = label_genes,
                                                        side = "bottom"))
  
  
  
  # RdBu Brewer
  colfun <- generate_rd_bu_colfun(-2, 2)
  
  
  
  plt <- Heatmap(heatmap_matrix_scaled,
                 name = "Log2FC",
                 left_annotation = annotation_row,
                 bottom_annotation = gene_annotation,
                 show_column_names = FALSE,
                 cluster_rows = TRUE,
                 cluster_columns = FALSE,
                 use_raster = FALSE,
                 col = colfun,
                 row_names_max_width = unit(100, "cm"))
  
  return (plt)
  
}

#####
# D #
#####
build_fig_3d_1 <- function() {
  
  all_vsd_df <- readRDS("data/processed/data_1_fig_4d_1.rds")
  annotation <- readRDS("data/processed/data_2_fig_4d_1.rds")
  
  anno_col <- HeatmapAnnotation(Species = annotation$species,
                                Condition = annotation$condition,
                                col = list (
                                  Species = species_pal,
                                  Condition = hc_infl_pal
                                ))
  
  # RdBu Brewer
  colfun <- generate_rd_bu_colfun(-3, 3)
  
  plt <- Heatmap(all_vsd_df,
                 top_annotation = anno_col,
                 show_row_names = F,
                 clustering_method_columns = "ward.D2",
                 clustering_method_rows = "ward.D2",
                 col = colfun)
  
  return (plt)
  
}

# this function draws fig 4d on the current graphics device
draw_fig_3d_2 <- function() {
  
  all_vsd_df_sorted_rank <- readRDS("data/processed/data_1_fig_4d_2.rds")
  all_vsd_df_sorted_wide <- readRDS("data/processed/data_2_fig_4d_2.rds")
  annotation_sorted <- readRDS("data/processed/data_3_fig_4d_2.rds") %>%
    mutate(group = paste(condition, species, sep = "_"))
  
  anno_col <- HeatmapAnnotation(super_block =
                                  anno_empty(
                                    border = FALSE
                                  ),
                                block =
                                  anno_block(
                                    labels = c("Human",
                                               "Mouse",
                                               "Human", 
                                               "Mouse")
                                    )#,
                                #Species = annotation_sorted$species,
                                #Condition = annotation_sorted$condition,
                                #col = list (
                                #  Species = species_pal,
                                #  Condition = hc_infl_pal
                                #),
                                #show_legend = c(FALSE, FALSE)
                                )
  
  # RdBu Brewer
  colfun <- generate_rd_bu_colfun(-3, 3)
  
  all_vsd_df_sorted_wide <- all_vsd_df_sorted_wide[
    match(all_vsd_df_sorted_rank$symbol, rownames(all_vsd_df_sorted_wide)),
    match(annotation_sorted %>%
            pull(sample_id), colnames(all_vsd_df_sorted_wide))]

  
  
  Heatmap(all_vsd_df_sorted_wide,
               name = "Relative\nexpression",
               top_annotation = anno_col,
               show_row_names = F,
               show_column_names = F,
               cluster_columns = F,
               cluster_rows = F,
               col = colfun,
               column_split = annotation_sorted$group,
          column_title = c()) %>% draw()
  
  seekViewport("annotation_super_block_1")
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport("annotation_super_block_2")
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("annotation_super_block_3")
  loc3 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport("annotation_super_block_4")
  loc4 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gpar(fill = hc_infl_pal[["HC"]]))
  grid.text("Healthy", x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5)
  
  grid.rect(loc3$x, loc3$y, width = loc4$x - loc3$x, height = loc4$y - loc3$y, 
            just = c("left", "bottom"), gp = gpar(fill = hc_infl_pal[["INFL"]]))
  grid.text("Inflamed", x = (loc3$x + loc4$x)*0.5, y = (loc3$y + loc4$y)*0.5)
  
}

# this function draws fig 4d on a virtual pdf device (thats what grid.grabExpr does)
# and returns the resulting gTree object which can be drawn with grid.draw()
# the w and h args define the size of the virtual pdf device and are important
# for correct positioning of the top healthy/inflamed annotaion as well as
# row level gene annotations
build_fig_3d_2 <- function(w = 7, h = 7) {
  
  annotate_manual_list <- c("TRAF1",
                            "JUNB",
                            "NFKBIZ",
                            "IL1A",
                            "IL1B",
                            "IL1RN",
                            "CDK5R1",
                            "CXCR4")
  
  labels_from_3c <- build_fig_3c_1(return_labels = TRUE)
  annotate_manual_list[!annotate_manual_list %in% labels_from_3c]
  
  annotate_manual <- unique(c(annotate_manual_list, labels_from_3c))
  
  p_lfc_df <- readRDS("data/processed/fisher_p_lfc_df.rds")
  lcpm_sorted_rank <- readRDS("data/processed/data_1_fig_4d_2.rds")
  lcpm_sorted_wide <- readRDS("data/processed/data_2_fig_4d_2.rds")
  annotation_sorted <- readRDS("data/processed/data_3_fig_4d_2.rds") %>%
    mutate(group = paste(condition, species, sep = "_"))
  
  anno_col <- HeatmapAnnotation(super_block =
                                  anno_empty(
                                    border = FALSE
                                  ),
                                block =
                                  anno_block(
                                    labels = c("Hs",
                                               "Mm",
                                               "Hs", 
                                               "Mm"),
                                    labels_gp = gpar(fontsize = 6)
                                  ), height = unit(1, "cm")#,
                                #Species = annotation_sorted$species,
                                #Condition = annotation_sorted$condition,
                                #col = list (
                                #  Species = species_pal,
                                #  Condition = hc_infl_pal
                                #),
                                #show_legend = c(FALSE, FALSE)
  )
  

  
  # RdBu Brewer
  colfun <- generate_rd_bu_colfun(
    quantile(lcpm_sorted_wide %>% as.matrix(), 0.01),
    quantile(lcpm_sorted_wide %>% as.matrix(), 0.99)
    )
  
  set.seed(42)
  permuted_cols <- 
    c(
      annotation_sorted %>%
        filter(group == "HC_Hs") %>%
        pull(sample_id) %>%
        sample,
      annotation_sorted %>%
        filter(group == "HC_Mm") %>%
        pull(sample_id) %>%
        sample,
      annotation_sorted %>%
        filter(group == "INFL_Hs") %>%
        pull(sample_id) %>%
        sample,
      annotation_sorted %>%
        filter(group == "INFL_Mm") %>%
        pull(sample_id) %>%
        sample
    )
  
  stopifnot(all(permuted_cols %in% annotation_sorted$sample_id))
  stopifnot(all(annotation_sorted$sample_id %in% permuted_cols))
  
  lcpm_sorted_wide <- lcpm_sorted_wide[
    match(lcpm_sorted_rank$symbol, rownames(lcpm_sorted_wide)),
    match(permuted_cols, colnames(lcpm_sorted_wide))]
  
  
  

  # old approach
  # p_values <- lcpm_sorted_wide %>% 
  #   data.frame %>% 
  #   rownames_to_column("symbol") %>% 
  #   melt %>%
  #   merge(annotation_sorted[,c("sample_id", "condition")], by.x = "variable", by.y = "sample_id", all.x=T) %>%
  #   group_by(symbol) %>%
  #   wilcox_test(value ~ condition) %>%
  #   adjust_pvalue(method = "fdr") %>%
  #   add_significance %>%
  #   arrange(match(symbol, rownames(lcpm_sorted_wide))) %>%
  #   dplyr::select(p.adj) %>%
  #   unlist
  
  # new, more consistent fisher annotation
  p_values <- p_lfc_df %>%
    filter(symbol %in% rownames(lcpm_sorted_wide)) %>%
    arrange(match(symbol, rownames(lcpm_sorted_wide))) %>%
    pull (fisher_adjusted)
  
  annotate_genes <- p_lfc_df %>%
    filter(symbol %in% rownames(lcpm_sorted_wide)) %>%
    slice_min(order_by = fisher_adjusted, n = 20, with_ties = F) %>%
    pull(symbol)
  
  annotate_all <- unique(c(annotate_genes, annotate_manual))
  
  anno_right <- rowAnnotation(
    Gene = anno_mark(
      at = match(annotate_all, rownames(lcpm_sorted_wide)),
      labels = annotate_all,
      labels_gp = gpar(fontsize = 6, fontface = "italic"),
      link_gp = gpar(lwd = 0.5)
      )
    )
  
  
  # pvalue_col_fun = colorRamp2(c(1.30103, 3), c("white", "green"))
  pvalue_col_fun <- generate_or_rd_colfun(min(-log10(p_values)), max(-log10(p_values[p_values != 0])))
  anno_p <- rowAnnotation(
    Sig. = anno_simple(ifelse(p_values <= 0.05, -log10(p_values), NA),
                       col = pvalue_col_fun, na_col = "grey70"), 
    annotation_name_side = "top", show_annotation_name = T, width = unit(0.5, "mm"),
    annotation_name_gp = gpar(fontsize = 8, fontface = "bold"))
  
  # for some reason, gt_render is needed for title_gp to work
  lgd_pvalue = Legend(title = expression(bold(paste("Fisher","'", "s"~bolditalic(P), sep=""))),
                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                      col_fun = pvalue_col_fun,
                      at = c(50,
                             100,
                             150,
                             200,
                             250),
                      labels = c("1e-50",
                                 "1e-100",
                                 "1e-150",
                                 "1e-200",
                                 "1e-250"),
                      labels_gp = gpar(fontsize = 8),
                      border = T
                      )
  # lgd_notsig <- Legend(labels = c("Not significant"), legend_gp = gpar(fill = "grey70"), title = "")
  
  # inflamm_legend <- packLegend(lgd_pvalue, lgd_notsig)
  inflamm_legend <- lgd_pvalue
  
  lcpm_sorted_wide
  
  plt <- grid.grabExpr({
    
    Heatmap(lcpm_sorted_wide,
            # name = "Z-score",
            top_annotation = anno_col,
            right_annotation = c(anno_p, anno_right),
            show_row_names = F,
            show_column_names = F,
            cluster_columns = F,
            cluster_rows = F,
            col = colfun,
            column_split = annotation_sorted$group,
            heatmap_legend_param = list(
              title = expression(bold(paste(bolditalic(z), "-score", sep=""))),
              legend_direction = "vertical", legend_height = unit(1.9, "cm"),
                                        grid_width = unit(0.25, "cm"), #title_position = "lefttop-rot",
                                        title_gp = gpar(fontsize = 8, fontface = "bold"),
                                        labels_gp = gpar(fontsize = 8),
                                        border = T),
            column_title = c()) %>% draw(annotation_legend_list = list(inflamm_legend), heatmap_legend_side = "right", merge_legend = T)
    
    seekViewport("annotation_super_block_1")
    loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
    seekViewport("annotation_super_block_2")
    loc2 = deviceLoc(x = unit(1, "npc"), y = unit(0.7, "npc"))
    
    seekViewport("annotation_super_block_3")
    loc3 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
    seekViewport("annotation_super_block_4")
    loc4 = deviceLoc(x = unit(1, "npc"), y = unit(0.7, "npc"))
    
    seekViewport("global")
    grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
              just = c("left", "bottom"), gp = gpar(fill = hc_infl_pal[["HC"]]))
    grid.text("Healthy", x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = gpar(fontsize = 7))
    
    grid.rect(loc3$x, loc3$y, width = loc4$x - loc3$x, height = loc4$y - loc3$y, 
              just = c("left", "bottom"), gp = gpar(fill = hc_infl_pal[["INFL"]]))
    grid.text("Inflamed", x = (loc3$x + loc4$x)*0.5, y = (loc3$y + loc4$y)*0.5, gp = gpar(fontsize = 7))
    
  }, width = w, height = h)
  return(plt)
} 


# landscape version

build_fig_3d_3 <- function() {
  
  annotate_manual_list <- c("TRAF1",
                            "JUNB",
                            "NFKBIZ",
                            "IL1A",
                            "IL1B",
                            "IL1RN",
                            "CDK5R1",
                            "CXCR4")
  
  labels_from_3c <- build_fig_3c_1(return_labels = TRUE)
  annotate_manual_list[!annotate_manual_list %in% labels_from_3c]
  
  annotate_manual <- unique(c(annotate_manual_list, labels_from_3c))
  
  p_lfc_df <- readRDS("data/processed/fisher_p_lfc_df.rds")
  lcpm_sorted_rank <- readRDS("data/processed/data_1_fig_4d_2.rds")
  lcpm_sorted_wide <- readRDS("data/processed/data_2_fig_4d_2.rds") %>% t()
  annotation_sorted <- readRDS("data/processed/data_3_fig_4d_2.rds") %>%
    mutate(group = paste(condition, species, sep = "_"))
  
  # RdBu Brewer
  colfun <- generate_rd_bu_colfun(
    quantile(lcpm_sorted_wide %>% as.matrix(), 0.01),
    quantile(lcpm_sorted_wide %>% as.matrix(), 0.99)
  )
  
  set.seed(42)
  permuted_rows <- 
    c(
      annotation_sorted %>%
        filter(group == "HC_Hs") %>%
        pull(sample_id) %>%
        sample,
      annotation_sorted %>%
        filter(group == "HC_Mm") %>%
        pull(sample_id) %>%
        sample,
      annotation_sorted %>%
        filter(group == "INFL_Hs") %>%
        pull(sample_id) %>%
        sample,
      annotation_sorted %>%
        filter(group == "INFL_Mm") %>%
        pull(sample_id) %>%
        sample
    )
  
  stopifnot(all(permuted_rows %in% annotation_sorted$sample_id))
  stopifnot(all(annotation_sorted$sample_id %in% permuted_rows))
  
  stopifnot(all(rownames(lcpm_sorted_wide) %in% permuted_rows))
  stopifnot(all(colnames(lcpm_sorted_wide) %in% lcpm_sorted_rank$symbol))
  
  lcpm_sorted_wide <- lcpm_sorted_wide[permuted_rows, lcpm_sorted_rank$symbol]
  
  
  p_values <- p_lfc_df %>%
    filter(symbol %in% colnames(lcpm_sorted_wide)) %>%
    arrange(match(symbol, colnames(lcpm_sorted_wide))) %>%
    pull (fisher_adjusted)
  
  annotate_genes <- p_lfc_df %>%
    filter(symbol %in% colnames(lcpm_sorted_wide)) %>%
    slice_min(order_by = fisher_adjusted, n = 20, with_ties = F) %>%
    pull(symbol)
  
  annotate_all <- unique(c(annotate_genes, annotate_manual))
  
  unique(annotate_manual[!annotate_manual %in% annotate_genes])
  
  anno_bottom <- columnAnnotation(
    Gene = anno_mark(
      at = match(annotate_all, colnames(lcpm_sorted_wide)),
      labels = annotate_all,
      labels_gp = gpar(fontsize = 6,
                       fontface = "italic"),
      padding = unit(0.2, "mm"),
      link_gp = gpar(lwd = 0.5),
      side = "bottom"
    )
  )
  
  
  # pvalue_col_fun = colorRamp2(c(1.30103, 3), c("white", "green"))
  pvalue_col_fun <- generate_or_rd_colfun(min(-log10(p_values)), max(-log10(p_values[p_values != 0])))
  anno_p <- columnAnnotation(
    fisher = anno_simple(
      ifelse(
        p_values <= 0.05,
        -log10(p_values),
        NA
        ),
      col = pvalue_col_fun,
      na_col = "grey70",
      border = T
      ), 
    # annotation_name_side = "top",
    show_annotation_name = F,
    height = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 8,
                              fontface = "bold")
    )
  
  # for some reason, gt_render is needed for title_gp to work
  lgd_pvalue = Legend(#title = expression(bold(paste("Fisher","'", "s"~bolditalic(P), sep=""))),
                      title = expression(bold("-log"[10]*"(Fisher's"~bolditalic(P)*")")),
                      title_gp = gpar(fontsize = 8,
                                      fontface = "bold"),
                      title_position = "lefttop",
                      col_fun = pvalue_col_fun,
                      at = c(50,
                             100,
                             150,
                             200,
                             250),
                      # labels = c("1e-50",
                      #            "1e-100",
                      #            "1e-150",
                      #            "1e-200",
                      #            "1e-250"),
                      labels = c("50",
                                 "100",
                                 "150",
                                 "200",
                                 "250"),
                      labels_gp = gpar(fontsize = 6),
                      labels_rot = 90,
                      legend_width = unit(2, "cm"),
                      grid_height = unit(0.25, "cm"),
                      border = T,
                      direction = "horizontal"
  )
  # lgd_notsig <- Legend(labels = c("Not significant"), legend_gp = gpar(fill = "grey70"), title = "")
  
  # inflamm_legend <- packLegend(lgd_pvalue, lgd_notsig)
  inflamm_legend <- lgd_pvalue
    
  plt <- Heatmap(lcpm_sorted_wide,
          # name = "Z-score",
          top_annotation = anno_p,
          bottom_annotation = anno_bottom,
          show_row_names = F,
          show_column_names = F,
          cluster_columns = F,
          cluster_rows = F,
          col = colfun,
          row_split = annotation_sorted$group,
          heatmap_legend_param = list(
            title = expression(bold(paste(bolditalic(z),
                                          "-score", sep=""))),
            legend_direction = "horizontal",
            legend_width = unit(2, "cm"),
            grid_height = unit(0.25, "cm"),
            # grid_width = unit(0.25, "cm"),
            title_position = "lefttop",
            title_gp = gpar(fontsize = 8,
                            fontface = "bold"),
            labels_gp = gpar(fontsize = 6),
            labels_rot = 90,
            border = T),
          row_title = c(
            "HC_Hs" = "Human\nHC",
            "HC_Mm" = "Mouse\nHC",
            "INFL_Hs" = "Human\nStim",
            "INFL_Mm" = "Mouse\nStim"
            )[unique(annotation_sorted$group)],
          row_title_gp = gpar(fontsize = 6,
                              fontface = "bold"),
          row_title_rot = 0,
          row_dend_side = "right",
          border = T) %>%
    draw(annotation_legend_list = list(inflamm_legend),
         heatmap_legend_side = "bottom",
         merge_legend = T)
  
  return(plt)
} 



#####
# E #
#####

build_fig_3e_1 <- function () {
  
  plt_nes_hm <- readRDS("data/processed/data_fig_4e_1.rds")
  
  col_annotation <- data.frame(names = colnames(plt_nes_hm)) %>%
    mutate(Species = str_remove(colnames(plt_nes_hm), "^.* ")) %>%
    column_to_rownames("names")
  
  colnames(plt_nes_hm) <- colnames(plt_nes_hm) %>%
    str_remove(" (Hs|Mm)$")
  
  rownames(col_annotation) <- rownames(col_annotation) %>%
    str_remove(" (Hs|Mm)$")
  
  plt_nes_hm <- plt_nes_hm[, match(colnames(plt_nes_hm), rownames(col_annotation))]
  
  plt <- ComplexHeatmap::pheatmap(plt_nes_hm,
                                  cellwidth = 30,
                                  cellheight = 30,
                                  border_color = NA,
                                  name = "NES",
                                  show_colnames = T,
                                  show_rownames = T,
                                  legend = T,
                                  cluster_rows = F,
                                  cluster_cols = F,
                                  annotation_col = col_annotation,
                                  annotation_colors = list(Species = species_pal),
                                  annotation_legend = F,
                                  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(25))
  
  return (plt)
  
}


build_fig_3e_2 <- function () {
  
  plt_sig_nes_hm <- readRDS("data/processed/data_fig_4e_2.rds") %>%
    t()
  
  row_annotation <- data.frame(names = rownames(plt_sig_nes_hm)) %>%
    mutate(Species = str_remove(rownames(plt_sig_nes_hm), "^.* ")) %>%
    column_to_rownames("names")
  
  rownames(plt_sig_nes_hm) <- rownames(plt_sig_nes_hm) %>%
    str_remove(" (Hs|Mm)$")
  
  rownames(row_annotation) <- rownames(row_annotation) %>%
    str_remove(" (Hs|Mm)$")
  
  plt_sig_nes_hm <- plt_sig_nes_hm[match(rownames(plt_sig_nes_hm), rownames(row_annotation)), ]

  
  # plt <- ComplexHeatmap::pheatmap(
  #   plt_sig_nes_hm,
  #   cellwidth = 6,
  #   cellheight = 6,
  #   border_color = NA,
  #   name = "NES",
  #   show_colnames = T,
  #   angle_col = "45",
  #   show_rownames = T,
  #   legend = T,
  #   cluster_rows = F,
  #   cluster_cols = F,
  #   annotation_row = row_annotation,
  #   annotation_colors = list(Species = species_pal),
  #   annotation_legend = T,
  #   annotation_names_row = F,
  #   heatmap_legend_param = list(legend_height = unit(1.5, "cm"),
  #                               grid_width = unit(0.25, "cm"),
  #                               title_gp = gpar(fontsize = 8,
  #                                               fontface = "bold"),
  #                               labels_gp = gpar(fontsize = 8),
  #                               border = T),
  #   fontsize_row = 6,
  #   fontsize_col = 5,
  #   fontsize = 6,
  #   color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(25))


  colfun <- generate_inferno_colfun(-3, 3)

  plt <- Heatmap(
    plt_sig_nes_hm,
    width = ncol(plt_sig_nes_hm)*unit(6, "pt"),
    height = nrow(plt_sig_nes_hm)*unit(6, "pt"),
    column_names_rot = -45,
    show_row_names = T,
    cluster_rows = F,
    cluster_columns = F,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 5),
    left_annotation = rowAnnotation(
      df = row_annotation,
      col = list(Species = species_pal),
      show_annotation_name = F,
      simple_anno_size = unit(0.25, "cm"),
      annotation_legend_param = list(
        at = c("Hs", "Mm"),
        labels = c("Human", "Mouse"),
        title_gp = gpar(fontsize = 8,
                        fontface = "bold"),
        grid_height = unit(0.25, "cm"),
        grid_width = unit(0.25, "cm"),
        direction = "horizontal",
        nrow = 2,
        title_position = "lefttop",
        labels_gp = gpar(fontsize = 6),
        border = T),
      border = T),
    heatmap_legend_param = list(
      title = "NES",
      legend_width = unit(1.4, "cm"),
      grid_height = unit(0.25, "cm"),
      title_gp = gpar(fontsize = 8,
                      fontface = "bold"),
      labels_gp = gpar(fontsize = 6),
      direction = "horizontal",
      title_position = "lefttop",
      border = T),
    col = colfun,
    border = T) %>%
    draw(heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)
  
  return (plt)
  
}

#####
# F #
#####

build_fig_3f <- function() {
  goseq.results.df <- readRDS("data/processed/data_fig_3f.rds") %>%
    mutate(species = ifelse(str_remove(study, "^.*_") == "Hs",
                            "Homo sapiens",
                            "Mus musculus")) %>%
    filter(category != "background")

  plt <- ggplot(data = goseq.results.df, aes(x = category, y = over_represented_p_score)) +
    theme_rgb() +
    geom_boxplot(lwd = 0.2) +
    geom_point(aes(fill = frac_de), pch = 21, size = 1.5, alpha = 0.7) +
    scale_fill_gradient2(low = "white", high = "darkred",
                         limits = c(0,1),
                         guide = guide_colorbar(
                           frame.colour = "black",
                           frame.linewidth = 1
                         )) +
    geom_line(aes(group = study), lwd = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.8) +
    # scale_fill_manual(values = colorRampPalette("redcols")) +
    # scale_fill_continuous(guide = guide_legend(override.aes = list(alpha = 0) ) ) +
    labs(fill = "Fraction DE") +
    xlab("Neutrotime category") +
    ylab("*P*-Score") +
    scale_x_discrete(labels = c("Early", "Late")) +
    theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(),
          legend.title = element_text(face = "bold"),
          legend.key.size = unit(0.3, "cm"),
          text = element_text(size = 8),
          strip.text.x = element_text(size = 8),
          axis.text.x = element_markdown(size = 8),
          axis.text.y = element_markdown(size = 8)) +
    facet_wrap(~species, ncol = 2)
  
  
  return(plt)
}



build_fig_3f_1 <- function() {
  goseq.results.df <- readRDS("data/processed/data_fig_3f.rds")
  ind_Hs <- ifelse(str_sub(goseq.results.df$study,start=-2) == "Hs", TRUE, FALSE)
  plt <- ggplot(data = goseq.results.df[goseq.results.df$category != "background" & ind_Hs,], aes(x = category, y = over_represented_p_score)) +
    theme_rgb() +
    geom_boxplot(lwd = 0.2) +
    geom_point(aes(fill = frac_de), pch = 21, size = 1.5, alpha = 0.7) +
    scale_fill_gradient2(low = "white", high = "darkred",
                         limits = c(0,1), guide = guide_legend(override.aes = list(alpha = 0))) +
    geom_line(aes(group = study), lwd = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.8) +
    # scale_fill_manual(values = colorRampPalette("redcols")) +
    # scale_fill_continuous(guide = guide_legend(override.aes = list(alpha = 0) ) ) +
    labs(fill = "Fraction DE") +
    xlab("Neutrotime category") +
    ylab("*P*-Score") +
    scale_x_discrete(labels = c("Early", "Late")) +
    theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(),
          legend.title = element_text(color = "transparent", face = "bold", size = 8),
          legend.text = element_text(color = "transparent"),
          text = element_text(size = 8),
          axis.text.x = element_markdown(size = 8),
          axis.text.y = element_markdown(size = 8),
          legend.background=element_rect(fill = alpha("white", 0)),
          legend.key=element_rect(fill = alpha("white", 0), color = alpha("white", 0))) +
    guides(color = guide_legend(override.aes = list(fill = NA)))
  
  
  return(plt)
}

build_fig_3f_2 <- function() {
  goseq.results.df <- readRDS("data/processed/data_fig_3f.rds")
  ind_Hs <- ifelse(str_sub(goseq.results.df$study,start=-2) == "Hs", TRUE, FALSE)
  plt <- ggplot(data = goseq.results.df[goseq.results.df$category != "background" & !ind_Hs,], aes(x = category, y = over_represented_p_score)) +
    geom_boxplot(lwd = 0.2) +
    geom_point(aes(fill = frac_de), pch = 21, size = 1.5, alpha = 0.7) +
    scale_fill_gradient2(low = "white", high = "darkred",
                         limits = c(0,1)) +
    geom_line(aes(group = study), lwd = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.8) +
    # scale_fill_manual(values = colorRampPalette("redcols")) +
    labs(fill = "Fraction DE") +
    xlab("Neutrotime category") +
    ylab("*P*-Score") +
    scale_x_discrete(labels = c("Early", "Late")) +
    theme_rgb() +
    theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.title = element_text(face = "bold", size = 8),
          legend.key.size = unit(0.3, "cm"),
          text = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))
  return(plt)
}

###########
# G and H #
###########
filter_counts_by_list <- function(gene_list){
  cross_se <- readRDS("data/processed/cross_se.rds")
  cross_se_inflammatory <- cross_se[, cross_se$inflammatory_analysis]
  annotation <- readRDS("data/processed/data_2_fig_4d_1.rds")
  # series_annotation <- readxl::read_xlsx("data/metadata/archs4_sample_annotation.xlsx")
  series_annotation <- readxl::read_xlsx("data/metadata/inflamm_salmon_metadata_full.xlsx")
  
  counts_df <- log2(assay(cross_se_inflammatory, assay = "abundance")+1)
  counts_df <- data.frame(counts_df)
  counts_df$gene <- rownames(counts_df)
  
  counts_df <- counts_df %>% filter(gene %in% gene_list)
  counts_df <- melt(counts_df, id = "gene")
  counts_df <- merge(counts_df, annotation, all.x, by.x = "variable", by.y = "sample_id")
  counts_df <- merge(counts_df, series_annotation[,c("SRX_accession", "series")], all.x=T, by.x = "variable", by.y = "SRX_accession")
  counts_df$condition <- ifelse(counts_df$species == "Hs", paste0(counts_df$condition, "_Hs"), paste0(counts_df$condition, "_Mm"))
  counts_df$condition <- factor(counts_df$condition, levels = c("HC_Hs", "INFL_Hs", "HC_Mm", "INFL_Mm"))
  
  return(counts_df)
}

build_fig_4g <- function() {
  counts_df <- filter_counts_by_list(c("TRAF1", "JUNB", "NFKBIZ", "IL1A", "IL1B", "IL1RN"))
  counts_df <- filter_counts_by_list(c("TRAF1", "CXCL10", "NFKBIZ", "IL1A", "IL1B", "IL1RN"))
  
  comparisons <- list(
    c('HC_Hs', 'INFL_Hs'),
    c('HC_Mm', 'INFL_Mm')
  )
  
  plt <- ggplot(data = counts_df, aes(x = condition, y = value)) +
    geom_jitter(size = 0.05) +
    geom_violin(fill = NA, lwd = 0.2) +
    theme_rgb() +
    xlab("Condition and species") +
    ylab("Log<sub>2</sub>(TPM+1)") +
    stat_compare_means(method = "wilcox.test",
                       comparisons = comparisons,
                       exact = F,
                       label = 'p.signif',
                       size = 1.5,
                       bracket.size = 0.25,
                       vjust = 0.4
    ) +
    theme(axis.title.x = element_text(),
          axis.title.y = element_markdown(),
          text = element_text(size = 8),
          axis.text.x = element_text(size = 7, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 8)) +
    facet_wrap(~gene, scales = "free_y")
  
  return(plt)
}

build_fig_4h <- function() {
  counts_df <- filter_counts_by_list(c("CDK5R1", "CDKN1B", "FADD", "CXCR4", "KLF3"))
  
  comparisons <- list(
    c('HC_Hs', 'INFL_Hs'),
    c('HC_Mm', 'INFL_Mm')
  )
  
  plt <- ggplot(data = counts_df, aes(x = condition, y = value)) +
    geom_jitter(size = 0.05) +
    geom_violin(fill = NA, lwd = 0.2) +
    theme_rgb() +
    xlab("Condition and species") +
    ylab("Log<sub>2</sub>(TPM+1)") +
    stat_compare_means(method = "wilcox.test",
                       comparisons = comparisons,
                       exact = F,
                       label = 'p.signif',
                       size = 1.5,
                       bracket.size = 0.25,
                       vjust = 0.4
    ) +
    theme(axis.title.x = element_text(),
          axis.title.y = element_markdown(),
          text = element_text(size = 8),
          axis.text.x = element_text(size = 7, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 8)) +
    facet_wrap(~gene, scales = "free_y")
  
  return(plt)
}



############
# FIGURE 4 #
############
############

#####
# A #
#####
build_fig_4a <- function(return_labels = FALSE) {
  
  chea_res_dn_df_all <- readRDS("data/processed/data_fig_5a_2.rds")
  
  highlight <- c(
    "FOXO3",
    "FOXO1",
    "TFEB",
    "SPI1",
    "RARA",
    "STAT5B",
    "ATF7"
  )
  
  # Define genes to highlight
  top10 <- chea_res_dn_df_all %>%
    transmute(qname = `Query Name`,
              species = str_replace(qname, "^.*_", ""),
              TF,
              logscore = -log10(as.numeric(Score))) %>%
    group_by(species, TF) %>%
    summarise(mean_logscore = mean(logscore)) %>%
    group_by(TF) %>%
    summarise(mean_logscore_TF = mean(mean_logscore)) %>%
    arrange(desc(mean_logscore_TF)) %>%
    mutate(rank = 1:length(TF),
           label = ifelse(rank <= 10, TF, NA)) %>%
    dplyr::select(label) %>% unlist %>% na.omit %>% unique
  
  print(paste0("Labeled TFs: (1) Top 10 by mean logscore. (2) Manual: ",   paste(highlight[!highlight %in% top10], collapse = ", ")))
  highlight_all <- unique(c(highlight, top10))
  
  
  plt_df <- chea_res_dn_df_all %>%
    transmute(qname = `Query Name`,
              species = str_replace(qname, "^.*_", ""),
              TF,
              logscore = -log10(as.numeric(Score))) %>%
    group_by(species, TF) %>%
    summarise(mean_logscore = mean(logscore)) %>%
    ungroup() %>%
    pivot_wider(names_from = "species",
                values_from = "mean_logscore") %>%
    arrange(TF %in% highlight_all)
  
  max_score <- max(c(plt_df$Hs, plt_df$Mm))
  min_score <- min(c(plt_df$Hs, plt_df$Mm))
  plt <- ggplot (plt_df, aes(x = Hs,
                y = Mm,
                label = ifelse(TF %in% highlight_all, TF, ""),
                #label = TF,
                fill = TF %in% highlight_all)) +
    geom_point(shape = 21, color = "black", size = 1) +
    stat_cor(aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~"),color = NULL, fill = NULL),
             size = 2,
             method = "pearson",
             cor.coef.name = "r",
             label.x.npc = "middle",
             label.y.npc = "bottom",
             ) +
    # coord_fixed() +
    # xlab(expression("Log"["10"](average ChEA3 score^-1,), Homo sapiens)) +
    # ylab(expression(paste("Log(average ChEA3 score"^-1,"), Mus musculus"))) +
    xlab(bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), human")) +
    ylab(bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), mouse")) +
    xlim(c(min_score, max_score)) +
    ylim(c(min_score, max_score)) +
    scale_fill_manual(values = c("darkgreen",
                                 "darkred")) +
    geom_label_repel(#force = 50,
      color = "darkred", alpha = 0.9,
      fill = "white",
      max.overlaps = Inf,
      size = 1.5,
      seed = 42,
      label.padding = unit(0.1, "lines"),
      box.padding = unit(0.1, "lines"),
      min.segment.length = 0,
      segment.size = 0.25) +
    theme_rgb() +
    theme(#aspect.ratio = 1,
          legend.position = "none",
          axis.title.x = element_text(),
          axis.title.y = element_text())
  
  if(return_labels == TRUE){
    return(highlight_all)
  }else{
    return(plt)
  }
  
  
}

#####
# B #
#####
build_fig_4b <- function(return_labels = FALSE) {
  
  chea_res_up_df_all <- readRDS("data/processed/data_fig_5a_1.rds")
  
  highlight <- c(
    "PLSCR1",
    "NFIL3",
    "CSRNP1",
    "FOSB",
    "FOS",
    "NFKB2",
    "CEBPB",
    "JUNB"
  )
  
  top10 <- chea_res_up_df_all %>%
    transmute(qname = `Query Name`,
              species = str_replace(qname, "^.*_", ""),
              TF,
              logscore = -log10(as.numeric(Score))) %>%
    group_by(species, TF) %>%
    summarise(mean_logscore = mean(logscore)) %>%
    group_by(TF) %>%
    summarise(mean_logscore_TF = mean(mean_logscore)) %>%
    arrange(desc(mean_logscore_TF)) %>%
    mutate(rank = 1:length(TF),
           label = ifelse(rank <= 10, TF, NA)) %>%
    dplyr::select(label) %>% unlist %>% na.omit %>% unique
  

  print(paste0("Labeled TFs: (1) Top 10 by mean logscore. (2) Manual: ",   paste(highlight[!highlight %in% top10], collapse = ", ")))
  highlight_all <- unique(c(highlight, top10))
  
  plt_df <- chea_res_up_df_all %>%
    transmute(qname = `Query Name`,
              species = str_replace(qname, "^.*_", ""),
              TF,
              logscore = -log10(as.numeric(Score))) %>%
    group_by(species, TF) %>%
    summarise(mean_logscore = mean(logscore)) %>%
    ungroup() %>%
    pivot_wider(names_from = "species",
                values_from = "mean_logscore") %>%
    arrange(TF %in% highlight_all)
  
  max_score <- max(c(plt_df$Hs, plt_df$Mm))
  min_score <- min(c(plt_df$Hs, plt_df$Mm))
  plt <- ggplot (plt_df,
                 aes(x = Hs,
                     y = Mm,
                     label = ifelse(TF %in% highlight_all, TF, ""),
                     #label = TF,
                     fill = TF %in% highlight_all)) +
    geom_point(shape = 21, color = "black", size = 1) +
    stat_cor(aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~"),color = NULL, fill = NULL),
             size = 2,
             method = "pearson",
             cor.coef.name = "r",
             label.x.npc = "middle",
             label.y.npc = "bottom") +
    # coord_fixed() +
    # xlab(expression(paste("Log(average ChEA3 score"^-1,"), Homo sapiens"))) +
    # ylab(expression(paste("Log(average ChEA3 score"^-1,"), Mus musculus"))) +
    xlab(bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), human")) +
    ylab(bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), mouse")) +
    xlim(c(min_score, max_score)) +
    ylim(c(min_score, max_score)) +
    scale_fill_manual(values = c("darkgreen",
                                 "darkred")) +
    geom_label_repel(#force = 50,
                     color = "darkred", alpha = 0.9,
                     fill = "white",
                     max.overlaps = Inf,
                     size = 1.5,
                     #nudge_y = -0.11,
                     seed = 42,
                     label.padding = unit(0.1, "lines"),
                     box.padding = unit(0.1, "lines"),
                     min.segment.length = 0,
                     segment.size = 0.25) +
    theme_rgb() +
    theme(#aspect.ratio = 1,
          legend.position = "none",
          axis.title.x = element_text(),
          axis.title.y = element_text())
  
  if(return_labels == TRUE){
    return(highlight_all)
  }else{
    return(plt)
  }  
  
}

#####
# C #
#####

# Loaded globally to enable usage in 4c.
# mark_tfs_sub_hm <- c("ATF3",
#                      "STAT3",
#                      "JUNB",
#                      "JUN",
#                      "NFKB1",
#                      "RELA",
#                      "MYC",
#                      "BATF3",
#                      "MAFF",
#                      "ATF7",
#                      "NANOGNB",
#                      "DUXA",
#                      "DMRT1",
#                      "HOMEZ",
#                      "JRK",
#                      "TCF20",
#                      "FOXO3",
#                      "FOXO1",
#                      "TFEB",
#                      "SPI1",
#                      "STAT5B",
#                      "PLSCR1",
#                      "NFIL3",
#                      "CSRNP1",
#                      "FOSB",
#                      "FOS",
#                      "NFKB2")

# mark_tfs_sub_hm <- c(
#   "FOXO3",
#   "FOXO1",
#   "TFEB",
#   "SPI1",
#   "RARA",
#   "STAT5B",
#   "PLSCR1",
#   "NFIL3",
#   "CSRNP1",
#   "FOSB",
#   "FOS",
#   "NFKB2"
# )

build_fig_4c <- function(label_manual = c("CEBPB", "JUNB")) {
  # Read the TF data
  chea_res_dn_df_all <- readRDS("data/processed/data_fig_5a_2.rds")
  chea_res_up_df_all <- readRDS("data/processed/data_fig_5a_1.rds")
  p_lfc_rank_stats <- readRDS("data/processed/p_lfc_rank_stats.rds")
  
  # Summarize the TF data
  chea_res_dn_df_all %<>%
    transmute(qname = `Query Name`,
              species = str_replace(qname, "^.*_", ""),
              TF,
              logscore = -log10(as.numeric(Score))) %>%
    group_by(species, TF) %>%
    summarise(mean_logscore = mean(logscore)) %>%
    ungroup() %>%
    pivot_wider(names_from = "species",
                values_from = "mean_logscore") 
  
  chea_res_up_df_all %<>%
    transmute(qname = `Query Name`,
              species = str_replace(qname, "^.*_", ""),
              TF,
              logscore = -log10(as.numeric(Score))) %>%
    group_by(species, TF) %>%
    summarise(mean_logscore = mean(logscore)) %>%
    ungroup() %>%
    pivot_wider(names_from = "species",
                values_from = "mean_logscore") 
  
  volcano_stats <- readRDS("data/processed/data_fig_4b.rds")
  
  volcano_stats_up <- volcano_stats[volcano_stats$lfc >= 0,]
  
  volcano_stats_up <- merge(volcano_stats_up, chea_res_up_df_all, by.x = "symbol", by.y = "TF", all.x = T)
  
  volcano_stats_dn <- volcano_stats[volcano_stats$lfc < 0,]
  
  volcano_stats_dn <- merge(volcano_stats_dn, chea_res_dn_df_all, by.x = "symbol", by.y = "TF", all.x = T)
  
  volcano_stats_tf <- rbind(volcano_stats_up, volcano_stats_dn)
  
  if(nrow(volcano_stats_tf) != nrow(volcano_stats)){
    stop("Merge problem.")
  }
  
  # volcano_stats_tf$mean_tf <- (volcano_stats_tf$Hs + volcano_stats_tf$Mm)/2
  
  fisher_logp_cutoff <- p_lfc_rank_stats[p_lfc_rank_stats$fisher_rank == 500, "p_fisher_log"]
  
  volcano_stats_tf %<>%
    filter(!is.na(pval)) %>%
    mutate(mean_tf = (Hs+Mm)/2) %>%
    filter(!is.na(mean_tf)) %>%
    mutate(pval = case_when(
      pval == 0 ~ min(pval*NA^(pval <= 0), na.rm=T),
      TRUE ~ pval)
    ) %>%
    mutate(logp = -log10(pval)) %>%
    mutate(color = ifelse(logp >= fisher_logp_cutoff, logp, 0)) %>%
    arrange(desc(lfc)) %>%
    mutate(rank_lfc = 1:nrow(.)) %>%
    arrange(pval) %>%
    mutate(rank_p = 1:nrow(.)) %>%
    arrange(desc(pval))
  
  print(length(unique(volcano_stats_tf$symbol)))
  
  # ifelse(p_values <= 0.05, -log10(p_values), NA), col = pvalue_col_fun, na_col = "red")
  
  volcano_stats_tf$label <- ifelse(volcano_stats_tf$rank_lfc <= 10 & volcano_stats_tf$pval <= 0.05 | 
                                     volcano_stats_tf$rank_lfc >= max(volcano_stats_tf$rank_lfc)-9 & volcano_stats_tf$pval <= 0.05 | 
                                     volcano_stats_tf$rank_p <= 5 |
                                     volcano_stats_tf$symbol %in% label_manual,
                                   volcano_stats_tf$symbol, NA)
  
  
  
    ggplot(volcano_stats_tf, aes(x = lfc,
                                 y = mean_tf,
                                 size = color,
                                 fill = color,
                                 label = label)) +
    scale_fill_gradientn(colors = c("grey70",
                                    brewer.pal(9, "OrRd")[4:9]),
                         values = rescale(c(0,-log10(0.05),-log10(0.0001))),
                         breaks = c(0,
                                    100,
                                    150,
                                    200,
                                    250),
                         labels = c("ns",
                                    "1e-100",
                                    "1e-150",
                                    "1e-200",
                                    "1e-250")) +
    scale_size_continuous(breaks = c(0,
                                     100,
                                     150,
                                     200,
                                     250),
                          labels = c("ns",
                                     "1e-100",
                                     "1e-150",
                                     "1e-200",
                                     "1e-250")) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               alpha = 0.8) +
    geom_point(color = "black", pch = 21, alpha = 0.9) +
    geom_label_repel(show.legend = F,
                     aes(size = 2),
                     color = "darkred",
                     fill = "white",
                     label.padding = unit(0.1, "lines"),
                     max.overlaps = Inf,
                     seed = 42,
                     fontface = "italic",
                     min.segment.length = 0,
                     segment.size = 0.25) +
    # geom_text() +
    xlab("Log<sub>2</sub> fold change") +
    ylab(bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), across species")) +
    #xlim(-3,3) +
    labs(fill = expression(italic(P)*"-value"),
         size = expression(italic(P)*"-value")) +
    guides(fill=guide_legend(),
           size = guide_legend()) +
    theme_rgb() +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_text(),
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 6),
          legend.position = "right",
          #aspect.ratio = 0.4
          )
}


build_fig_4d <- function() {
  
  plt_data <- readRDS("data/processed/data_1_fig_5c.rds")
  col_annotation <- readRDS("data/processed/data_2_fig_5c.rds")
  row_annotation <- readRDS("data/processed/data_3_fig_5c.rds")
  
  plt_data_scaled <- plt_data %>%
    t() %>%
    scale() %>%
    t()
  
  # viridis inferno
  
  lo <- quantile(plt_data_scaled, 0.01)
  hi <- quantile(plt_data_scaled, 0.99)
  col_fun2 <- generate_inferno_colfun(lo, hi)
  
  # anno_col <- HeatmapAnnotation(`Gene sets` = ifelse(col_annotation$set == "up", "Upregulated", "Downregulated"),
  #                               col = list(
  #                                 `Gene sets` = c(Upregulated = "Orange", Downregulated = "Lightblue")
  #                               ),
  #                               show_legend = FALSE)
  
  anno_col <- HeatmapAnnotation(`Gene sets` = anno_block(
    gp = gpar(
      fill = c(Upregulated = "Orange", Downregulated = "Lightblue")
    ),
    labels = c("Downregulated", "Upregulated"),
    labels_gp = gpar(fontsize = 8)
  ))
  
  # Annotate p values
  # old: calculate
  # p_values <- plt_data_scaled %>% 
  #   data.frame %>% 
  #   rownames_to_column("symbol") %>% 
  #   melt %>%
  #   mutate(group = case_when(
  #     grepl("_up", variable) ~ "UP",
  #     grepl("_dn", variable) ~ "DN"
  #   )) %>%
  #   group_by(symbol) %>%
  #   wilcox_test(value ~ group) %>%
  #   adjust_pvalue(method = "fdr") %>%
  #   add_significance %>%
  #   arrange(match(symbol, rownames(plt_data_scaled))) %>%
  #   dplyr::select(p.adj) %>%
  #   unlist()
  
  # new: take from row_annotation
  p_values <- row_annotation$padj
  
  labels_auto <- c(
    row_annotation %>%
      slice_min(n = 10, order_by = group_score, with_ties = F) %>%
      rownames_to_column("TF") %>%
      pull(TF),
    row_annotation %>%
      slice_max(n = 10, order_by = group_score, with_ties = F) %>%
      rownames_to_column("TF") %>%
      pull(TF)
  )
  labels_4ab <- unique(c(build_fig_4a(return_labels = TRUE), build_fig_4b(return_labels = TRUE)))
  labels_all <- unique(
    c(
      labels_4ab,
      labels_auto
    )
  )
  
  annotate_genes <- unique(c(labels_all , rownames(row_annotation[row_annotation$tf_group == "b_eq",])))
  anno_right <- rowAnnotation(
    Gene = anno_mark(
      at = match(annotate_genes, rownames(plt_data)),
      labels = annotate_genes,
      labels_gp = gpar(fontsize = 6))
    )  
  
  pvalue_col_fun <- generate_or_rd_colfun(min(p_values), max(p_values))
  anno_p <- rowAnnotation(
    Sig. = anno_simple(
      ifelse(p_values <= 0.05, p_values, NA),
      col = pvalue_col_fun, na_col = "grey70",
      border = TRUE),
    annotation_name_side = "top",
    show_annotation_name = T,
    annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
    )
  
  lgd_pvalue = Legend(title = expression(bold(paste("Adjusted"~bolditalic(P), "-value", sep=""))),
                      col_fun = pvalue_col_fun,
                      at = c(0.05,
                             0.025,
                             0.001),
                      labels = c("0.05",
                                 "0.025",
                                 expression(""<="0.001")),
                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                      title_position = "lefttop-rot",
                      labels_gp = gpar(fontsize = 6),
                      border = TRUE
                      )
  
  # lgd_notsig <- Legend(labels = c("Not significant"),
  #                      legend_gp = gpar(fill = "grey70"),
  #                      title = "")
  
  # hm_legend <- packLegend(lgd_pvalue, lgd_notsig)
  hm_legend <- lgd_pvalue
  
  plt <- Heatmap(plt_data_scaled,
                 name = "Relative TF activity",
                 cluster_rows = F,
                 cluster_columns = F,
                 show_row_names = F,
                 show_column_names = F,
                 col = col_fun2,
                 row_split = row_annotation$tf_group,
                 row_title = c(paste0("TFs enriched\nin up DEGs, N = ", sum(row_annotation$tf_group == "a_up")), paste0("TFs enriched\nin down DEGs, N = ", sum(row_annotation$tf_group == "c_dn"))),
                 # row_title = c(paste0("TFs enriched\nin up DEGs, N = ", sum(row_annotation$tf_group == "a_up")), paste0("TFs enriched\nin both groups, N = ", sum(row_annotation$tf_group == "b_eq")), paste0("TFs enriched\nin down DEGs, N = ", sum(row_annotation$tf_group == "c_dn"))),                 row_title = c(paste0("TFs enriched\nin up DEGs, N = ", sum(row_annotation$tf_group == "a_up")), paste0("TFs enriched\nin both groups, N = ", sum(row_annotation$tf_group == "b_eq")), paste0("TFs enriched\nin down DEGs, N = ", sum(row_annotation$tf_group == "c_dn"))),
                 # row_title = c("TFs enriched\nin up DEGs", "TFs enriched\nin down DEGs"),
                 row_title_gp = gpar(fontsize = 8),
                 row_gap = unit(2.5, "mm"),
                 column_split = col_annotation$set,
                 column_title = " ",
                 top_annotation = anno_col,
                 right_annotation = c(anno_p, anno_right),
                 heatmap_legend_param = list(
                   labels_gp = gpar(fontsize = 6),
                   title_gp = gpar(fontsize = 8, fontface = "bold"),
                   title_position = "lefttop-rot",
                   border = TRUE
                 ),
                 border = TRUE)
  #draw(plt, annotation_legend_list = list(hm_legend), merge_legend = TRUE)
  
  return (draw(plt,
               annotation_legend_list = list(hm_legend),
               merge_legend = TRUE))
  
}




############
# FIGURE 5 #
############
############

# Initialize function to make to volcano plot for subfigures

make_volcano <- function(res, annotate_genes, annotate_n, title = "untitled", ratio = 1){
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  res_df <- data.frame(res) %>% rownames_to_column("gene") %>% na.omit
  
  volcano_label <- data.frame(gene = res_df$gene, label = ifelse(res_df$gene %in% annotate_genes, res_df$gene, NA ))
  volcano_label$top <- ifelse(res_df$gene %in% annotate_genes[1:annotate_n], res_df$gene, NA)
  volcano_label$italics <- ifelse(is.na(volcano_label$top), NA, paste0( "italic('", volcano_label$top, "')"))
  #subset df for highlighting
  volcano_highlight <- res_df %>%
    # data.frame %>%
    # rownames_to_column %>%
    filter(gene %in% volcano_label$label, !is.na(padj)) %>%
    column_to_rownames("gene")
  volcano_highlight[volcano_highlight$padj == 0,"padj"] <- min(data.frame(res_df)[data.frame(res_df)$padj != 0,"padj"], na.rm=T)*10^-1
  
  # Order to plot highlights above rest
  res_df %<>%
    mutate(to_highlight = gene %in% rownames(volcano_highlight)) %>%
    arrange(to_highlight)
  
  volcano_label %<>%
    arrange(match(gene, res_df$gene))
  
  print(paste0("Percentage of core_up genes sig up in this comparison: ", nrow(res_df[res_df$log2FoldChange >= 1 & res_df$padj <= 0.05 & res_df$gene %in% up_fisher_mouse,])/length(up_fisher_mouse)))
  print(paste0("Percentage of core_up genes sig up in this comparison: ", nrow(res_df[res_df$log2FoldChange <= -1 & res_df$padj <= 0.05 & res_df$gene %in% up_fisher_mouse,])/length(up_fisher_mouse)))

  # Volcano
  
  keyvals <- ifelse(res_df$gene %in% rownames(volcano_highlight), "red",
    ifelse(abs(res_df$log2FoldChange) >= 1 & res_df$padj <= 0.05, "grey70",
                    ifelse(abs(res_df$log2FoldChange) <= 1 & res_df$padj <= 0.05, "grey70", 
                           ifelse(abs(res_df$log2FoldChange) >= 1 & res_df$padj > 0.05, "grey70",
                                  "grey70"))))
  
  # keyvals <- ifelse(res_df$padj <= 0.05, "blue", "red")
  
  # names(keyvals)[keyvals == 'yellow'] <- 'Core inflammation genes'
  # names(keyvals)[keyvals == 'red'] <- 'P-value and log2(FC)'
  # names(keyvals)[keyvals == 'blue'] <- 'P-value'
  # names(keyvals)[keyvals == 'green'] <- 'Log2(FC)'
  # names(keyvals)[keyvals == 'grey'] <- 'NS'
  
  # All grey except core genes
  names(keyvals)[keyvals == 'red'] <- 'Core inflammation genes'
  names(keyvals)[keyvals == 'grey70'] <- 'P-value and log2(FC)'
  names(keyvals)[keyvals == 'grey70'] <- 'P-value'
  names(keyvals)[keyvals == 'grey70'] <- 'Log2(FC)'
  names(keyvals)[keyvals == 'grey70'] <- 'NS'
  
  sizevals <- ifelse(keyvals == "red", 0.3, 0.1)
  alphavals <- ifelse(keyvals == "red", 1, 0.5)
  
# Custom enhancedVolcano function to use new vars "labAlpha" and "labPad"
  plt <- enhancedVolcanoFAR(
    res_df,
    x = "log2FoldChange",
    y = "padj",
    lab = volcano_label$italics,
    boxedLabels = T,
    parseLabels = T,
    drawConnectors = T,
    max.overlaps = Inf,
    pointSize = sizevals,
    colCustom = keyvals,
    colAlpha = alphavals,
    # shape = 21,
    labSize = 2,
    labAlpha = 0.9,
    labPad = 0.1,
    title = title,
    titleLabSize = 0,
    legendLabSize = 7,
    axisLabSize = 5,
    subtitle = "",
    pCutoff = 0.05,
    FCcutoff = 1,
    FCcutoff_label = 1,
    subtitleLabSize = 0,
    caption = "",
    captionLabSize = 0,
    legendPosition = "none"
  ) +
    theme_classic() +
    theme(legend.position="none", plot.title = element_markdown(size = 8, margin = margin(0,0,-10,0)),
          axis.title = element_text(size = 8), aspect.ratio = ratio)
  return(plt)
}

# Initialize ggvenn helper function for 5a. 
run_ggvenn <- function(X, title = "unnamed"){
  ggvenn(X, fill_color = c("red", "yellow", "blue"), set_name_size = 2, text_size = 2, show_percentage = FALSE, set_name_color = "transparent",
         stroke_size = 0.2) + theme(plot.title = element_text(hjust = 0.5, size = 6, margin = margin(0,0,-5,0))) + ggtitle(title)
}

build_fig_5a_stat <- function(return_lgd = FALSE) {
  mat <- readRDS("data/processed/data_fig_5b.rds")
  mat %<>% data.frame() %>%
    dplyr::rename("1" = "BL_BM",
                  "2" = "MEM_BM", 
                  "3" = "AP_BM", 
                  "4" = "MEM_BL") %>%
    dplyr::mutate(rows = c(2,3,4,5)) %>%
    remove_rownames() %>%
    column_to_rownames("rows") %>% as.matrix()
  
  colfun <- generate_rd_bu_colfun(0.05, 0)
  
  pvalue_col_fun = colorRamp2(c(0.05, 0), c("white", "palegreen")) 
  
  lgd = Legend(labels = c(expression(""<="0.05"), expression("">"0.05")), title = "Significance", legend_gp = gpar(fill = c("green", "white")), border="black")
  
  plt <- Heatmap(mat,
                 # col = pvalue_col_fun,
                 cluster_columns = F,
                 cluster_rows = F, 
                 rect_gp = gpar(type = "none"),
                 show_heatmap_legend = F,
                 show_column_names = T,
                 show_row_names = T,
                 column_names_gp = gpar(fontsize = 4),
                 row_names_gp = gpar(fontsize = 4),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.rect(x = x, y = y, width = width, height = height, 
                             gp = gpar(col = NA, fill = NA))
                   if(i == j) {
                     grid.circle(x = x, y = y, r = abs(0.9)/2 * min(unit.c(width, height)), 
                                 gp = gpar(lwd = 0.5, fill = pvalue_col_fun(mat[i, j]), col = "black"))
                     grid.text(signif(mat[i, j], 2), x, y, gp = gpar(fontsize = 3))
                     
                   } else if(i > j) {
                     grid.circle(x = x, y = y, r = abs(0.9)/2 * min(unit.c(width, height)), 
                                 gp = gpar(lwd = 0.5, fill = pvalue_col_fun(mat[i, j]), col = "black"))
                     grid.text(signif(mat[i, j], 2), x, y, gp = gpar(fontsize = 3))
                   } else {
                     grid.rect(x = x, y = y, width = width, height = height, 
                               gp = gpar(col = NA, fill = NA))
                   }
                 }) #%>% draw(annotation_legend_list = list(lgd))
  
  if(return_lgd == FALSE){
    return(plt)
  }else{
    return(lgd)
  }
}

build_fig_5a <- function(return_lgd = FALSE) {
  list_comparisons <- readRDS("data/processed/data_fig_5a.rds")
  fig_5a_stat <- build_fig_5a_stat(return_lgd = FALSE)
  fig_5a_stat_grob <- grid.grabExpr(draw(fig_5a_stat))
  
  plt1 <- arrangeGrob(run_ggvenn(list_comparisons[[5]],
                                 names(list_comparisons)[5]),
                      fig_5a_stat_grob,
                      ncol = 2)
  
  layout_matrix <- rbind(
    c(1,2),
    c(3,4)
  )
  plt2 <- arrangeGrob(run_ggvenn(list_comparisons[[1]], names(list_comparisons)[1]),
                      run_ggvenn(list_comparisons[[2]], names(list_comparisons)[2]),
                      run_ggvenn(list_comparisons[[3]], names(list_comparisons)[3]),
                      run_ggvenn(list_comparisons[[4]], names(list_comparisons)[4]),
                      layout_matrix = layout_matrix)
  
  plt <- arrangeGrob(
    plt2,
    plt1,
    nrow = 2,
    heights = c(2,1)
    )
  
  if(return_lgd == FALSE){
    return(plt)
  }else{
    lgd <- ggplot(data=data.frame(x = c(1,2,3,4,5), y = c(1,2,3,4,5), col = c("A", "B", "C", "D", "E")), aes (x=x, y=y, fill = col)) +
      geom_point(pch=21) +
      scale_fill_manual(labels = c("\u2264 0.05", "> 0.05", "Core inflammation genes", "Genes with increased accessibility", "Genes with decreased accessibility"), 
                        values = c("green", "white", "red", "yellow", "blue")) +
      labs(fill = NULL) +
      guides(fill = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 0.5))) +
      theme(legend.text = element_text(size = 6),
            legend.background=element_rect(fill = alpha("white", 0)),
            legend.title = element_blank(),
            legend.margin=margin(c(0,0,0,0)),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.key.height = unit(0,"mm"))
    return(get_legend(lgd))
  }
}

# Function for Figure 5b and Sx
plot_motif_enrichment <- function(comparison, ylab, y_axis, title = "untitled"){
  
  if(y_axis == "lfc"){
    current_label <- comparison %>%
      arrange(desc(lfc)) %>%
      mutate(lfc_rank = 1:nrow(comparison)) %>%
      arrange(desc(log_p_value_score)) %>%
      mutate(log_p_rank = 1:nrow(comparison)) %>%
      filter(log_p_rank %in% 1:10 | log_p_rank %in% (nrow(comparison)-9):nrow(comparison))
  } else if(y_axis == "mean_tf") {
    current_label <- comparison %>%
      arrange(desc(mean_tf)) %>%
      mutate(tf_rank = 1:nrow(comparison)) %>%
      arrange(desc(log_p_value_score)) %>%
      mutate(log_p_rank = 1:nrow(comparison)) %>%
      filter(log_p_rank %in% 1:10 | log_p_rank %in% (nrow(comparison)-9):nrow(comparison))
  } else{
    stop("Unknown y-axis")
  }
  
  comparison$col <- ifelse(comparison$fdr <= 0.05, "darkred", "grey70")
  
  plt <- comparison %>%
    ggplot(aes(x = log_p_value_score,
               y = get(y_axis),
               label = symbol)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dotted") +
    geom_point(pch = 16, color = comparison$col, size = 1) +
    geom_text_repel(
      data = current_label,
      size = 2,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.1
      ) +
    xlab(bquote("Directional log<sub>2</sub>(*P*)-score, motif enrichment")) +
    ylab(ylab) +
    stat_cor(
      cor.coef.name = "rho",
      method = "spearman", size = 2,
      aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~"),color = NULL)) +
    ggtitle(title) +
    theme_rgb() +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          legend.position = "bottom",
          aspect.ratio = 1)
  return(plt)
}

build_fig_5b <- function() {
  load("data/processed/data_fig_5b.rda")
  set.seed(42)
  plt <- plot_motif_enrichment(comparison = TF_merged_AP_VS_BL, 
                               ylab = bquote("Log<sub>2</sub>(FC), cross-species"),
                               y_axis = "lfc",
                               title = "Air pouch vs. blood")
  return(plt)
  
}

build_fig_S6a <- function(full_list = FALSE, height, width) {
  p_lfc_df <- readRDS("data/processed/fisher_p_lfc_df.rds")
  ortho_mouse_to_human <- readRDS("data/processed/ortho_mouse_to_human.rds")
  
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  list_comparisons <- readRDS("data/processed/data_fig_5a.rds")
  get_atac_vector <- function(n) {
    return(ifelse(up_fisher_mouse %in% list_comparisons[[n]][[2]] & up_fisher_mouse %in% list_comparisons[[n]][[3]], "Both",
                  ifelse(up_fisher_mouse %in% list_comparisons[[n]][[2]], "Increased",
                         ifelse(up_fisher_mouse %in% list_comparisons[[n]][[3]], "Decreased", "None"))))
  }
  
  df <- data.frame(
    gene = up_fisher_mouse
  ) %>% 
    mutate("1" = get_atac_vector(1),
           "2" = get_atac_vector(2),
           "3" = get_atac_vector(3),
           "4" = get_atac_vector(4),
           "5" = get_atac_vector(5)) %>%
    left_join(ortho_mouse_to_human[,c("external_gene_name","hsapiens_homolog_associated_gene_name")],
              by = c("gene" = "external_gene_name")) %>%
    left_join(p_lfc_df, by = c("hsapiens_homolog_associated_gene_name" = "symbol")) %>%
    arrange(fisher_adjusted) %>%
    dplyr::select(-hsapiens_homolog_associated_gene_name, -mean_lfc, -fisher_adjusted) %>%
    column_to_rownames("gene") 
    # as.matrix()
  
  collect.counts <- data.frame()  
  for(i in 1:nrow(df)){
    count_increased <- sum(df[i,] == "Increased")
    count_decreased <- sum(df[i,] == "Decreased")
    count_both <- sum(df[i,] == "Both")
    count_none <- sum(df[i,] == "None")
    current.counts <- c(count_increased, count_decreased, count_both, count_none)
    collect.counts <- rbind(collect.counts, current.counts)
  }
  colnames(collect.counts) <- c("count_increased", "count_decreased", "count_both", "count_none")
  
  df %<>%
    cbind(collect.counts) %>%
    arrange(desc(count_increased), desc(count_both), desc(count_none), desc(count_decreased)) %>%
    dplyr::select(-c(count_increased, count_both, count_none, count_decreased)) %>%
    as.matrix()
  
  rdbu <- brewer.pal(9, "RdBu")
  
  colors = structure(
    c(
      rdbu[1],
      rdbu[9],
      rdbu[3],
      rdbu[5]
      ),
    names = c(
      "Increased",
      "Decreased",
      "Both",
      "None"
      )
    )
  
  plt <- grid.grabExpr(
    draw(
      Heatmap(
        df,
        name = "Accessibility",
        col = colors,
        show_row_names = ifelse(full_list == FALSE, F, T),
        column_names_gp = gpar(fontsize = 6),
        border_gp = gpar(col = "black", lty = 1),
        row_title = "Core genes",
        row_title_gp = gpar(
          fontface = "bold",
          fontsize = 8
          ),
        cluster_rows = F,
        heatmap_legend_param = list(
          title_gp = gpar(
            fontsize = 8,
            fontface = "bold"
            ),
          grid_height = unit(0.25, "cm"),
          grid_width = unit(0.25, "cm"),
          direction = "horizontal",
          nrow = 2,
          title_position = "lefttop",
          labels_gp = gpar(fontsize = 6),
          border = T
          )
        ),
      heatmap_legend_side = "bottom"),
    height = height,
    width = width
    )
  
  # Make table
  t <- apply(df, 2, table)
  order_ind <- match(c("Increased", "Decreased", "Both", "None"), rownames(t))
  t <- t[order_ind,]
  t
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.5)),
    colhead = list(fg_params=list(cex = 0.5)),
    rowhead = list(fg_params=list(cex = 0.5)))
  plt_table <- tableGrob(t, theme = mytheme)
  return(list(plt, plt_table))
}

build_fig_S9_2 <- function() {
  load("data/processed/data_fig_5b.rda")
  plt1 <- plot_motif_enrichment(comparison = TF_merged_BL_VS_BM, 
                               ylab = bquote("-Log<sub>2</sub>(fold change), cross-species"),
                               y_axis = "lfc",
                               title = "Blood vs. Bone marrow")
  plt2 <- plot_motif_enrichment(comparison = TF_merged_MEM_VS_BM, 
                                ylab = bquote("-Log<sub>2</sub>(fold change), cross-species"),
                                y_axis = "lfc",
                                title = "Membrane vs. Bone marrow")
  plt3 <- plot_motif_enrichment(comparison = TF_merged_AP_VS_BM, 
                                ylab = bquote("-Log<sub>2</sub>(fold change), cross-species"),
                                y_axis = "lfc",
                                title = "Air pouch vs. Bone marrow")
  plt4 <- plot_motif_enrichment(comparison = TF_merged_MEM_VS_BL, 
                                ylab = bquote("-Log<sub>2</sub>(fold change), cross-species"),
                                y_axis = "lfc",
                                title = "Membrane vs. Blood")
  
  plt5 <- plot_motif_enrichment(comparison = TF_merged_BL_VS_BM, 
                                ylab = bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), across species"),
                                y_axis = "mean_tf",
                                title = "Blood vs. Bone marrow")
  
  plt6 <- plot_motif_enrichment(comparison = TF_merged_MEM_VS_BM, 
                                ylab = bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), across species"),
                                y_axis = "mean_tf",
                                title = "Membrane vs. Bone marrow")
  
  plt7 <- plot_motif_enrichment(comparison = TF_merged_AP_VS_BM, 
                                ylab = bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), across species"),
                                y_axis = "mean_tf",
                                title = "Air pouch vs. Bone marrow")
  
  plt8 <- plot_motif_enrichment(comparison = TF_merged_MEM_VS_BL, 
                                ylab = bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), across species"),
                                y_axis = "mean_tf",
                                title = "Membrane vs. Blood")
  
  plt9 <- plot_motif_enrichment(comparison = TF_merged_AP_VS_BL, 
                                ylab = bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), across species"),
                                y_axis = "mean_tf",
                                title = "Air pouch vs. Blood")
  
  plt_list <- list(plt1,plt2,plt3,plt4,plt5,plt6,plt7,plt8,plt9)
  
  return(plt_list)
  
}

build_fig_5d_by_row <- function(row = 1) {
  load("data/processed/data_fig_5cd.rda")
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  # dn_fisher_mouse <- readRDS("data/processed/dn_fisher_mouse.rds")
  
  if (row == 1) {
    toptable <- res_WT
    
    comp_string <- "WT: Zymosan vs control"
    p_down <- signif(goseq_WT_up$under_represented_pvalue, 2)
    p_up <- signif(goseq_WT_up$over_represented_pvalue, 2)
    
  } else if (row == 2) {
    toptable <- res_JUNB
    
    comp_string <- "Resting: *JunB*<sup>-/-</sup> vs WT"
    p_down <- signif(goseq_JUNB_dn$over_represented_pvalue, 2)
    p_up <- signif(goseq_JUNB_dn$under_represented_pvalue, 2)
    
  } else if (row == 3) {
    toptable <- res_CEBPB
    
    comp_string <- "Resting: *Cebp&#946;*<sup>-/-</sup> vs WT"
    p_down <- signif(goseq_CEBPB_dn$over_represented_pvalue, 2)
    p_up <- signif(goseq_CEBPB_dn$under_represented_pvalue, 2)
    
  } else if (row == 4) {
    toptable <- res_JUNB_activated
    
    comp_string <- "Zymosan: *JunB*<sup>-/-</sup> vs WT"
    p_down <- signif(goseq_JUNB_activated_dn$over_represented_pvalue, 2)
    p_up <- signif(goseq_JUNB_activated_dn$under_represented_pvalue, 2)
    
  } else if (row == 5) {
    toptable <- res_CEBPB_activated
    
    comp_string <- "Zymosan: *Cebp&#946;*<sup>-/-</sup> vs WT"
    p_down <- signif(goseq_CEBPB_activated_dn$over_represented_pvalue, 2)
    p_up <- signif(goseq_CEBPB_activated_dn$under_represented_pvalue, 2)
    
  } else {
    stop("Invalid row definition")
  }
  
  toptable %<>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(
      p_fixed = ifelse(
        padj == 0,
        min(padj[padj != 0], na.rm = T),
        padj
      ),
      scaled_p = scale(-log10(p_fixed)),
      scaled_lfc = scale(log2FoldChange),
      rank_stat = ifelse(
        gene %in% up_fisher_mouse,
        abs(scaled_p) * abs(scaled_lfc),
        -1
        ),
      toprank_fisher = rank_stat >= sort(rank_stat, decreasing = T)[12],
      label = ifelse(
        toprank_fisher,
        gene,
        NA
      ),
      color = ifelse(
        gene %in% up_fisher_mouse,
        "in_up",
        "not_in_up"
        )
      ) %>%
    arrange(desc(color))
  
  n_up_in_fisher_up <- toptable %>%
    filter(
      padj <= 0.05,
      log2FoldChange >= 1,
      gene %in% up_fisher_mouse
    ) %>%
    nrow()
  
  n_down_in_fisher_up <- toptable %>%
    filter(
      padj <= 0.05,
      log2FoldChange <= -1,
      gene %in% up_fisher_mouse
    ) %>%
    nrow()
  set.seed(42)
  plt <- ggplot(toptable,
                aes(
                  x = log2FoldChange,
                  y = -log10(p_fixed),
                  color = color,
                  label = label
                  ),
                ) +
    geom_point(size = 0.25) +
    scale_color_manual(
      breaks = c("in_up", "not_in_up"),
      values = c("darkred", "grey70")
    ) +
    geom_text_repel(
      aes(point.size = 0.25),
      color = "black",
      size = 2,
      fontface = "italic",
      min.segment.length = 0,
      segment.size = 0.1,
      max.overlaps = Inf
      ) +
    theme_rgb() +
    ylim(c(0, 350)) +
    ggtitle("") +
    xlab("Log<sub>2</sub> fold change") +
    ylab("-log<sub>10</sub>(<i>P</i><sub>adjusted</sub>)") +
    theme(
      legend.position = "none",
      title = element_blank(),
      axis.title.x = element_markdown(),
      axis.title.y = element_markdown()
    )
    # theme(
    #   axis.text = element_text(size = 6),
    #   axis.title = element_text(size = 6)
    # )
  
  
  title_string <- paste0(
    comp_string,
    "<br>",
    "<i>P<sub>downregulated</sub></i> = ",
    p_down,
    " <i>P<sub>upregulated</sub></i> = ",
    p_up,
    "<br>",
    n_down_in_fisher_up,
    " downregulated, ",
    n_up_in_fisher_up,
    " upregulated")
  
  title_grob <- richtext_grob(
    title_string,
    gp = gpar(fontsize = 8),
    halign = 0,
    valign = 0
  )
  
  
  return (list(title_grob, plt))
  
}



build_fig_5c <- function (return_lgd = FALSE, ratio = 1) {
  load("data/processed/data_fig_5cd.rda")
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  
  n_up <- nrow(
    res_WT %>%
    data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange >= 1,
             gene %in% up_fisher_mouse
           )
    )
  n_dn <- nrow(
    res_WT %>%
      data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange <= -1,
             gene %in% up_fisher_mouse
      )
  )
  
  

  plt <- make_volcano(res = res_WT, annotate_genes = up_fisher_mouse, annotate_n = 20,
               title = paste0("WT: Zymosan vs control<br>*P<sub>downregulated</sub>* = ", signif(goseq_WT_up$under_represented_pvalue, 2),
                              ", *P<sub>upregulated</sub>* = ", signif(goseq_WT_up$over_represented_pvalue, 2),
                              "<br>", n_dn, " downregulated, ", n_up, " upregulated"), ratio = ratio)
  
  
  if(return_lgd == FALSE){
    return(plt)
  }else{
  lgd <- ggplot(data=data.frame(x = c(1,2), y = c(1,2), col = c("A", "B")), aes (x=x, y=y, fill = col)) +
    geom_point(pch=21) +
    # scale_fill_manual(labels = c("NS", "Log<sub>2</sub>(FC)", "*P*-value", "*P*-value and log<sub>2</sub>(FC)", "Core inflammation genes"), 
    #                   values = c("grey30", "forestgreen", "royalblue", "red2", "yellow")) +
    scale_fill_manual(labels = c("Other", "Core inflammation genes"), 
                      values = c("grey70", "red")) +
    labs(fill = NULL) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 0.5, color = c("grey70", "red")))) +
    theme(legend.text = element_markdown(size = 6),
          legend.background=element_rect(fill = alpha("white", 0)),
          legend.title = element_blank(),
          legend.margin=margin(c(0,0,0,0)),
          legend.spacing.x = unit(0, "mm"),
          legend.spacing.y = unit(0, "mm"),
          legend.key.height = unit(0,"mm"))
  return(get_legend(lgd))
  }
}

build_fig_s6c_1 <- function (ratio = 1) {
  load("data/processed/data_fig_5cd.rda")
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  
  n_up <- nrow(
    res_JUNB %>%
      data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange >= 1,
             gene %in% up_fisher_mouse
      )
  )
  n_dn <- nrow(
    res_JUNB %>%
      data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange <= -1,
             gene %in% up_fisher_mouse
      )
  )
  
  plt <- make_volcano(res = res_JUNB, annotate_genes = up_fisher_mouse, annotate_n = 20,
                      title = paste0("Resting: *JunB*<sup>-/-</sup> vs WT<br>*P<sub>downregulated</sub>* = ", signif(goseq_JUNB_dn$over_represented_pvalue,2),
                                     " ,*P<sub>upregulated</sub>* = ", signif(goseq_JUNB_dn$under_represented_pvalue,2),
                                     "<br>", n_dn, " downregulated, ", n_up, " upregulated"), ratio = ratio)
  return(plt)
}

build_fig_s6c_2 <- function (ratio = 1) {
  load("data/processed/data_fig_5cd.rda")
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  
  n_up <- nrow(
    res_CEBPB %>%
      data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange >= 1,
             gene %in% up_fisher_mouse
      )
  )
  n_dn <- nrow(
    res_CEBPB %>%
      data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange <= -1,
             gene %in% up_fisher_mouse
      )
  )
  
  plt <- make_volcano(res = res_CEBPB, annotate_genes = up_fisher_mouse, annotate_n = 20,
                      title = paste0("Resting: *Cebp&#946;*<sup>-/-</sup> vs WT<br>*P<sub>downregulated</sub>* = ", signif(goseq_CEBPB_dn$over_represented_pvalue,2),
                                     " ,*P<sub>upregulated</sub>* = ", signif(goseq_CEBPB_dn$under_represented_pvalue,2),
                                     "<br>", n_dn, " downregulated, ", n_up, " upregulated"), ratio = ratio)
  return(plt)
}

build_fig_s6c_3 <- function (ratio = 1) {
  load("data/processed/data_fig_5cd.rda")
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  
  n_up <- nrow(
    res_JUNB_activated %>%
      data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange >= 1,
             gene %in% up_fisher_mouse
      )
  )
  n_dn <- nrow(
    res_JUNB_activated %>%
      data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange <= -1,
             gene %in% up_fisher_mouse
      )
  )
  
  plt <- make_volcano(res = res_JUNB_activated, annotate_genes = up_fisher_mouse, annotate_n = 20,
                      title = paste0("Zymosan: *JunB*<sup>-/-</sup> vs WT<br>*P<sub>downregulated</sub>* = ", signif(goseq_JUNB_activated_dn$over_represented_pvalue,2),
                                     " ,*P<sub>upregulated</sub>* = ", signif(goseq_JUNB_activated_dn$under_represented_pvalue,2),
                                     "<br>", n_dn, " downregulated, ", n_up, " upregulated"), ratio = ratio)
  return(plt)
}

build_fig_s6c_4 <- function (ratio = 1) {
  load("data/processed/data_fig_5cd.rda")
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  
  n_up <- nrow(
    res_CEBPB_activated %>%
      data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange >= 1,
             gene %in% up_fisher_mouse
      )
  )
  n_dn <- nrow(
    res_CEBPB_activated %>%
      data.frame %>%
      rownames_to_column("gene") %>%
      filter(padj <= 0.05,
             log2FoldChange <= -1,
             gene %in% up_fisher_mouse
      )
  )
  
  plt <- make_volcano(res = res_CEBPB_activated, annotate_genes = up_fisher_mouse, annotate_n = 20,
                      title = paste0("Zymosan: *Cebp&#946;*<sup>-/-</sup> vs WT<br>*P<sub>downregulated</sub>* = ", signif(goseq_CEBPB_activated_dn$over_represented_pvalue,2),
                                     " ,*P<sub>upregulated</sub>* = ", signif(goseq_CEBPB_activated_dn$under_represented_pvalue,2),
                                     "<br>", n_dn, " downregulated, ", n_up, " upregulated"), ratio = ratio)
  return(plt)
}

build_fig_s6b <- function() {
  load("data/processed/data_1_fig_5d.rda")
  load("data/processed/data_fig_5cd.rda")
  
  
  plot_histo <- function(res_p, current.goseq, title, original_goseq, direction) {
    current.p <- signif(current.goseq$over_represented_pvalue, 2)
    plt <- res_p$hitsPerc %>% data.frame() %>%
      dplyr::rename("Percentage" = ".") %>%
      ggplot(aes(x = Percentage)) +
      # geom_vline(xintercept = 0.05) +
      geom_histogram(binwidth = 1, fill = "black", color = "transparent") +
      xlim(0,100) +
      # ggtitle(title) +
      xlab(paste0("Percentage ", direction)) +
      ylab("Count") +
      theme_rgb() +
      theme(
        axis.title.x = element_markdown(),
        axis.title.y = element_text(),
        # plot.title = element_markdown()
        title = element_blank()
        )
    plt2 <- plt +
      geom_segment(aes(x = original_goseq$hitsPerc, y = max(ggplot_build(plt)$data[[1]]$count)/10, xend = original_goseq$hitsPerc, yend = 0),
                   arrow = arrow(length = unit(0.1, "cm")), color = "darkred") +
      annotate(hjust = 0.05, geom = "text", size = 2, parse = T, label = deparse(bquote(italic(P) * "=" * .(current.p))),  x = original_goseq$hitsPerc, y = (max(ggplot_build(plt)$data[[1]]$count)/10)+max(ggplot_build(plt)$data[[1]]$count)/15) +
      annotate(angle = 90, geom = "text", size = 1.5, parse = F, label = paste0(round(min(res_p$hitsPerc), digits = 1), "%"), x = min(res_p$hitsPerc), y = (max(ggplot_build(plt)$data[[1]]$count)*0.9)) +
      annotate(angle = 90, geom = "text", size = 1.5, parse = F, label = paste0(round(max(res_p$hitsPerc), digits = 1), "%"), x = max(res_p$hitsPerc), y = (max(ggplot_build(plt)$data[[1]]$count)*0.9))
    return(plt2)
  }
  fig_s6b_1 <- plot_histo(res_WT_p, goseq_WT_up, "WT: Zymosan vs resting", goseq_WT_up, direction = "upregulated")
  fig_s6b_2 <- plot_histo(res_JUNB_p, goseq_JUNB_dn, "Resting: *JunB -/-* vs WT", goseq_JUNB_dn, direction = "downregulated")
  fig_s6b_3 <- plot_histo(res_CEBPB_p, goseq_CEBPB_dn,"Resting: *Cebp&#946;*<sup>-/-</sup> vs WT", goseq_CEBPB_dn, direction = "downregulated")
  fig_s6b_4 <- plot_histo(res_JUNB_activated_p, goseq_JUNB_activated_dn,"Zymosan: *JunB -/-* vs WT", goseq_JUNB_activated_dn, direction = "downregulated")
  fig_s6b_5 <- plot_histo(res_CEBPB_activated_p, goseq_CEBPB_dn,"Zymosan: *Cebp&#946;*<sup>-/-</sup> vs WT", goseq_CEBPB_activated_dn, direction = "downregulated")
  return(list(fig_s6b_1, fig_s6b_2, fig_s6b_3, fig_s6b_4, fig_s6b_5))
}



build_fig_5d <- function() {
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  dds <- readRDS("data/processed/data_fig_5e.rds")
  vsd <- vst(dds)
  expr <- as.data.frame(assay(vsd))
  expr <- expr[rownames(expr) %in% up_fisher_mouse,] %>%
    t() %>%
    scale(scale = FALSE) %>%
    t()
  
  expr_df <- expr %>%
    as.data.frame() %>%
    t() %>%
    scale() %>%
    t()
  
  annotation <- dds %>%
    colData() %>%
    as.data.frame() %>%
    transmute(id,
              treatment
    ) %>%
    filter(id %in% colnames(expr_df)) %>%
    arrange(match(id, colnames(expr_df)))
  
  expr_df_sorted <- expr_df %>%
    as.data.frame() %>%
    rownames_to_column("symbol") %>%
    pivot_longer(cols = -"symbol",
                 names_to = "id",
                 values_to = "count") %>%
    left_join(annotation, by = "id") %>%
    arrange(treatment)
  
  expr_df_sorted_rank <- expr_df_sorted %>%
    filter(treatment == "zymosan") %>%
    group_by(treatment, symbol) %>%
    summarise(mean_treatment = mean (count)) %>%
    arrange(desc(mean_treatment))
  
  expr_df_sorted_wide <- expr_df_sorted %>%
    pivot_wider(id_cols = c("symbol"),
                names_from = "id",
                values_from = "count") %>%
    column_to_rownames("symbol")
  
  annotation_sorted <- annotation %>%
    arrange(treatment)
  
  expr_df_sorted_wide <- expr_df_sorted_wide[
    match(expr_df_sorted_rank$symbol, rownames(expr_df_sorted_wide)),
    match(annotation_sorted %>%
            pull(id), colnames(expr_df_sorted_wide))]
  
  expr_df_sorted_wide %<>% t
  
  colfun <- generate_rd_bu_colfun(-2, 2)
  
  anno_col <- rowAnnotation(block = anno_block(labels = c("Control", "Zymosan"), labels_gp = gpar(fontsize = 6)))
  
  
  annotate_all <- up_fisher_mouse[1:20]
  
  anno_right <- columnAnnotation(Gene = anno_mark(
    at = match(annotate_all, colnames(expr_df_sorted_wide)),
    labels = annotate_all,
    labels_gp = gpar(fontsize = 6, fontface = "italic"),
    link_gp = gpar(lwd = 0.5),
    padding = 2))
  
  
  plt <- draw(Heatmap(expr_df_sorted_wide,
                      name = "Relative expression",
                      left_annotation = anno_col,
                      top_annotation =  anno_right,
                      show_row_names = F,
                      show_column_names = F,
                      row_title = c(),
                      cluster_columns = F,
                      cluster_rows = F,
                      col = colfun,
                      show_heatmap_legend = T,
                      row_split = annotation_sorted$treatment,
                      heatmap_legend_param = list(
                        direction = "horizontal",
                        title_gp = gpar(fontsize = 8),
                        labels_gp = gpar(fontsize = 8),
                        border = F
                        ),
                      border = TRUE),
                 heatmap_legend_side = "bottom")
  
  return(plt)
}

build_fig_5e <- function() {
  load("data/processed/data_1_fig_5e.rda")
  up_fisher_mouse <- readRDS("data/processed/up_fisher_mouse.rds")
  

  make_scatter <- function(x, y, xtitle, ytitle) {
    
    # color by condition: red = DE in both, black = DE in a, green = DE in b, grey = DE in none
    
    df <- y %>%
      left_join(x, by = "symbol") %>%
      filter(complete.cases(.[,2:5])) 
    
    df_core <- df %>%
      filter(symbol %in% up_fisher_mouse) %>%
    mutate(
      condition = ifelse(
        .[,2] >= 1 & .[,3] <= 0.05 & .[,4] >= 1 & .[,5] <= 0.05, "UP_both",
        ifelse(
          .[,2] <= -1 & .[,3] <= 0.05 & .[,4] <= -1 & .[,5] <= 0.05, "DN_both",
          ifelse(
          .[,2] >= 1 & .[,3] <= 0.05, "UP_x",
          ifelse(
            .[,2] <= -1 & .[,3] <= 0.05, "DN_x",
            ifelse(
              .[,4] >= 1 & .[,5] <= 0.05, "UP_y",
              ifelse(
                .[,4] <= -1 & .[,5] <= 0.05, "DN_y",
                "None"
              )
            )
          )
        )
      )
      ), color = ifelse(
        condition %in% c("UP_both", "DN_both"), "red",
        ifelse(
          condition %in% c("UP_x", "DN_x"), "black",
          ifelse(
            condition %in% c("UP_y", "DN_y"), "green",
            "grey70"
          )
        )
      ))
    
    df_label <- df %>%
      arrange(across(names(.)[2])) %>%
      mutate(rank_1 = 1:nrow(.)) %>%
      arrange(across(names(.)[3])) %>%
      mutate(rank_2 = 1:nrow(.)) %>%
      mutate(label = ifelse(rank_1 <= 10 | 
                              rank_1 >= (nrow(.)-9) | 
                              rank_2 <= 10 | 
                              rank_2 >= (nrow(.)-9), symbol, NA))
    
    
  
    
    plt <- df %>%
      ggplot(aes(x = get(colnames(x)[2]), y = get(colnames(y)[2]))) +
      geom_point(size = 0.2, color = "grey70") +
      geom_point(data = df_core, color = df_core$color, size = 0.6, alpha = 1, pch = 16) +
      scale_x_continuous(limits = symmetric_limits) +
      scale_y_continuous(limits = symmetric_limits) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_text_repel(data = df_label, aes(label = label), size = 1.5, fontface = "italic") +
      geom_text_npc(aes(npcx = x, npcy = y, label=label), 
                    data = data.frame(x = 0.05, y = 0.95, 
                                      label=paste(
                                        paste0("Upregulated in both: ", sum(df_core$condition == "UP_both")), 
                                                  paste0("Downregulated in both: ", sum(df_core$condition == "DN_both")),
                                                  sep = "\n")
                                      ), size = 2, color = "red") +
      geom_text_npc(aes(npcx = x, npcy = y, label=label), 
                    data = data.frame(x = 0.05, y = 0.83, 
                                      label=paste(
                                                  paste0("Upregulated in x: ", sum(df_core$condition == "UP_x")),
                                                  paste0("Downregulated in x: ", sum(df_core$condition == "DN_x")),
                                                  sep = "\n")
                    ), size = 2, color = "green") +
      geom_text_npc(aes(npcx = x, npcy = y, label=label), 
                    data = data.frame(x = 0.05, y = 0.71, 
                                      label=paste(
                                                  paste0("Upregulated in y: ", sum(df_core$condition == "UP_y")),
                                                  paste0("Downregulated in y: ", sum(df_core$condition == "DN_y")),
                                                  sep = "\n")
                    ), size = 2, color = "black") +
      geom_text_npc(aes(npcx = x, npcy = y, label=label), 
                    data = data.frame(x = 0.05, y = 0.59, 
                                      label=paste(
                                                  paste0("Not DE: ", sum(df_core$condition == "None")),
                                                  sep = "\n")
                    ), size = 2, color = "grey70") +
      xlab(xtitle) +
      ylab(ytitle) +
      theme_rgb() +
      theme(axis.title.x = element_markdown(),
            axis.title.y = element_markdown())
    return(plt)
    
  }
  
  fig_5e1 <- make_scatter(x = res_JUNB_activated_df, 
               y = res_JUNB_df, 
               xtitle = paste0("Zymosan, *JunB*<sup>-/-</sup> vs WT"),
               ytitle = paste0("Resting, *JunB*<sup>-/-</sup> vs WT"))
  fig_5e2 <- make_scatter(x = res_CEBPB_activated_df, 
               y = res_CEBPB_df, 
               xtitle = paste0("Zymosan, *Cebp&#946;*<sup>-/-</sup> vs WT"),
               ytitle = paste0("Resting, *Cebp&#946;*<sup>-/-</sup> vs WT"))
  

  fig_5e3 <- make_scatter(x = res_JUNB_zymosan_df, 
               y = res_WT_df, 
               xtitle = paste0("*JunB*<sup>-/-</sup>, zymosan vs resting"),
               ytitle = paste0("WT, zymosan vs resting"))
  
  fig_5e4 <- make_scatter(x = res_CEBPB_zymosan_df, 
               y = res_WT_df, 
               xtitle = paste0("*Cebp&#946;*<sup>-/-</sup>, zymosan vs resting"),
               ytitle = paste0("WT, zymosan vs. resting"))
  
  return(list(fig_5e1, fig_5e2, fig_5e3, fig_5e4))
}


  


############
# FIGURE 6 #
############
############

#####
# A #
#####
build_fig_6a <- function () {
  plt_surfaceome <- readRDS("data/processed/data_fig_6a.rds")
  
  up_genes <- readRDS("data/processed/fisher_up_genes.rds")
  dn_genes <- readRDS("data/processed/fisher_dn_genes.rds")
  
  annotation_col <- data.frame(group = ifelse(colnames(plt_surfaceome) %in% up_genes, "up",
                                              ifelse(colnames(plt_surfaceome) %in% dn_genes,
                                                     "dn", NA)))
  
  annotation_top <- columnAnnotation(
    block = anno_block(
      labels = c(paste0("Down\n (N = ", sum(annotation_col$group == "dn"), ")"), 
                 paste0("Up\n (N = ", sum(annotation_col$group == "up"), ")")),
      labels_gp = gpar(fontsize = 6),
      height = unit(6, "mm")
    )
  )
  
  
  annotation_row <- rowAnnotation(Species = str_remove(rownames(plt_surfaceome), "^.*_"),
                                  col = list (Species = species_pal),
                                  annotation_name_gp = gpar(fontsize = 10),
                                  annotation_legend_param = list(title = expression(bold("Species")),
                                                                 labels = c("Human", "Mouse"),
                                                                 at = c("Hs", "Mm"),
                                                                 title_gp = gpar(fontsize = 8),
                                                                 title_position = "lefttop",
                                                                 labels_gp = gpar(fontsize = 6),
                                                                 grid_height = unit(0.25, "cm"),
                                                                 grid_width = unit(0.25, "cm"),
                                                                 nrow = 1,
                                                                 border = TRUE),
                                  # show_legend = F,
                                  show_annotation_name = F,
                                  border = TRUE
                                  )
  
  rownames(plt_surfaceome) <- str_remove(rownames(plt_surfaceome), "_(Hs|Mm)$")
  rownames(plt_surfaceome) <- str_remove(rownames(plt_surfaceome), "condition_")
  rownames(plt_surfaceome) <- str_replace_all(rownames(plt_surfaceome), "_", " ")
  
  surfaceome_mean <- data.frame(protein = colnames(plt_surfaceome), meanlogFC = colMeans(plt_surfaceome))
  sum(surfaceome_mean$meanlogFC >= 0)
  sum(surfaceome_mean$meanlogFC < 0)
  
  # l <- length(colnames(plt_surfaceome))
  # label_genes <- colnames(plt_surfaceome)[c(1:5, 10, 15, 20,
  #                                           l-(0:4), l-9, l-14, l-19)]
  
  label_genes <- c("CD14",
                   "CD40",
                   "CD69",
                   "CD274",
                   "IL4R")
  
  # gene_annotation <- columnAnnotation(genes = anno_mark(at = which(colnames(plt_surfaceome) %in% label_genes),
  #                                                       labels = label_genes,
  #                                                       labels_gp = gpar(fontsize = 5),
  #                                                       link_gp = gpar(lwd = 0.5),
  #                                                       side = "bottom", padding = 2))
  
  col_ind <- data.frame(gene = colnames(plt_surfaceome), color = ifelse(colnames(plt_surfaceome) %in% label_genes, "red", "black"))
    
  gene_annotation <- columnAnnotation(genes = anno_mark(at = which(colnames(plt_surfaceome) %in% colnames(plt_surfaceome)),
                                                        labels = colnames(plt_surfaceome),
                                                        labels_gp = gpar(fontsize = 6,
                                                                         col = col_ind$color),
                                                        link_gp = gpar(lwd = 0.3),
                                                        side = "bottom", padding = 0.2))
  
  
  
  # RdBu Brewer
  colfun <- generate_rd_bu_colfun(-4, 4)
  
  Heatmap(plt_surfaceome,
          # name = "Log2FC",
          heatmap_legend_param = list(title = expression(bold("Log"["2"](FC))),
                                      direction = "horizontal",
                                      title_gp = gpar(fontsize = 8),
                                      title_position = "lefttop",
                                      grid_height = unit(0.25, "cm"),
                                      labels_gp = gpar(fontsize = 6),
                                      border = TRUE),
          # show_heatmap_legend = F,
          left_annotation = annotation_row,
          bottom_annotation = gene_annotation,
          top_annotation = annotation_top,
          show_column_names = FALSE,
          column_split = annotation_col$group,
          border_gp = gpar(col = "black", lty = 1),
          column_title = c(),
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          use_raster = FALSE,
          col = colfun,
          row_names_gp = gpar(fontsize = 6),
          row_names_max_width = unit(100, "cm"),
          border = TRUE) %>% draw(
            heatmap_legend_side = "bottom"
            )
}

#####
# B #
#####
build_fig_6b_1 <- function() {
  mfi_human <- readRDS("data/processed/data_fig_6c_1.rds")

  mfi_human <- mfi_human %>%
    filter (day %in% c(2),
            marker %in% c("CD14", "CD40", "CD69", "IL4R", "PDL1")) %>%
    mutate(group  = factor(paste(day, condition, sep = "_"), levels = c("2_unstim",
                                                                        "2_GL",
                                                                        "2_GI")))
  
  comparisons <- list(
    c('2_unstim', '2_GI'),
    c('2_unstim', '2_GL')
  )
  
  
  
  mfi_human %>% group_by(group) %>% shapiro_test(mfi) #sig = non-normality
  
  maxmfi <- mfi_human %>%
    group_by(marker) %>%
    summarise(y.position = max(mfi))
  
  kruskal.res <- mfi_human %>%
    filter(day == 2) %>%
    group_by(marker) %>%
    kruskal_test(mfi ~ group)
  
  stopifnot(!any(kruskal.res$p > 0.05))
  
  stat.test <- mfi_human %>%
    filter(day == 2) %>%
    group_by(marker) %>%
    dunn_test(mfi ~ group, p.adjust.method = "holm") %>%
    add_significance() %>% add_x_position(x = "group") %>%
    left_join(maxmfi) %>%
    filter(group1 == "2_unstim")
  
  
  stat.test[stat.test$group2 == "2_GL","y.position"] <- stat.test[stat.test$group2 == "2_GL","y.position"] * 1.2
  stat.test[stat.test$group2 == "2_GI","y.position"] <- stat.test[stat.test$group2 == "2_GI","y.position"] * 1.4
  set.seed(42)
  mfi_plt_human <- ggplot(mfi_human,
                          mapping = aes(x = group,
                                        y = mfi#,
                                        #color = as.character(day)
                          )
  ) +
    geom_violin(scale = "width", trim = F, fill = "lightblue") +
    geom_vline(xintercept = 1.5, linetype = "dashed", size = 0.3) +
    geom_jitter(size = 1, shape = 21, fill = "white", stroke = 0.25) +
    # stat_compare_means(method = "wilcox.test",
    #                    comparisons = comparisons,
    #                    exact = FALSE,
    #                    label = 'p.signif',
    #                    vjust = 0.5,
    #                    tip.length = 0, size = 1.5
    # ) +
    stat_pvalue_manual(stat.test,
                       vjust = 0.5,
                       tip.length = 0, size = 1.5) +
    theme_rgb() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          strip.background = element_blank(),
          axis.title.x = element_text(size = 8),
          strip.text = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(angle = 45,
                                     vjust = 0,
                                     hjust = 1,
                                     size = 6)) +
    facet_wrap(~marker,
               nrow = 1,
               scales = "free_y") +
    scale_x_discrete(
      breaks = c(
        "2_unstim",
        "2_GI",
        "2_GL"
      ),
      labels = c(
        "U",
        "GI",
        "GL"
      )
    ) +
    xlab("Conditions, human, day 2 after stimulation") +
    ylab("MFI, human")
  
  return(mfi_plt_human)
}

build_fig_6b_2 <- function () {
  mfi_mouse <- readRDS("data/processed/data_fig_6c_2.rds")
  
  mfi_mouse <- mfi_mouse %>%
    filter (day %in% c(2),
            marker %in% c("CD14", "CD40", "CD69", "IL4R", "PDL1")) %>%
  mutate(group  = factor(paste(day, condition, sep = "_"), levels = c("2_unstim",
                                                                      "2_GL",
                                                                      "2_GI")))
  
  comparisons <- list(
    c('2_unstim', '2_GI'),
    c('2_unstim', '2_GL')
  )
  
  mfi_mouse %>% group_by(group) %>% shapiro_test(mfi) #sig = non-normality
  
  maxmfi <- mfi_mouse %>%
    group_by(marker) %>%
    summarise(y.position = max(mfi))

  kruskal.res <- mfi_mouse %>%
    filter(day == 2) %>%
    group_by(marker) %>%
    kruskal_test(mfi ~ group)
  
  stopifnot(!any(kruskal.res$p > 0.05))
  
  stat.test <- mfi_mouse %>%
    filter(day == 2) %>%
    group_by(marker) %>%
    dunn_test(mfi ~ group, p.adjust.method = "holm") %>%
    add_significance() %>% add_x_position(x = "group") %>%
    left_join(maxmfi) %>%
    filter(group1 == "2_unstim")
  
  
  stat.test[stat.test$group2 == "2_GL","y.position"] <- stat.test[stat.test$group2 == "2_GL","y.position"] * 1.2
  stat.test[stat.test$group2 == "2_GI","y.position"] <- stat.test[stat.test$group2 == "2_GI","y.position"] * 1.4
  
  set.seed(42)
  mfi_plt_mouse <- ggplot(mfi_mouse,
                          mapping = aes(x = group,
                                        y = mfi#,
                                        #color = as.character(day)
                          )
  ) +
    geom_violin(scale = "width", trim = F, fill = "lightblue") +
    geom_vline(xintercept = 1.5, linetype = "dashed", size = 0.3) +
    geom_jitter(size = 1, shape = 21, fill = "white", stroke = 0.25) +
    # stat_compare_means(method = "wilcox.test",
    #                    comparisons = comparisons,
    #                    exact = FALSE,
    #                    label = 'p.signif',
    #                    vjust = 0.5,
    #                    tip.length = 0, size = 1.5) +
    stat_pvalue_manual(stat.test,
                       vjust = 0.5,
                       tip.length = 0, size = 1.5) +
    #stat_pvalue_manual(stat.test, tip.length = 0) +
    theme_rgb() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6,
                                     angle = 45,
                                     vjust = 0,
                                     hjust = 1)) +
    facet_wrap(~marker,
               nrow = 1,
               scales = "free_y") +
  scale_x_discrete(
    breaks = c(
      "2_unstim",
      "2_GI",
      "2_GL"
    ),
    labels = c(
      "U",
      "GI",
      "GL"
    )
  ) +
    xlab("Conditions, mouse, day 2 after stimulation") +
    ylab("MFI, mouse")
  
  return (mfi_plt_mouse)
  
}

build_fig_s8_1 <- function() {
  mfi_human <- readRDS("data/processed/data_fig_6c_1.rds")
  
  mfi_human <- mfi_human %>%
    filter (day %in% c(0, 2),
            condition == "unstim",
            marker %in% c("CD14", "CD40", "CD69", "IL4R", "PDL1")) %>%
    mutate(group  = factor(paste(day, condition, sep = "_"), levels = c("0_unstim",
                                                                        "2_unstim")))
  
  comparisons <- list(
    c('0_unstim', '2_unstim')
  )
  
  
  
  mfi_human %>% group_by(group) %>% shapiro_test(mfi) #sig = non-normality
  
  maxmfi <- mfi_human %>%
    group_by(marker) %>%
    summarise(y.position = max(mfi))
  
  kruskal.res <- mfi_human %>%
    # filter(day == 2) %>%
    group_by(marker) %>%
    kruskal_test(mfi ~ group)
  
  # stopifnot(!any(kruskal.res$p > 0.05))
  
  stat.test <- mfi_human %>%
    # filter(day == 2) %>%
    group_by(marker) %>%
    dunn_test(mfi ~ group, p.adjust.method = "holm") %>%
    add_significance() %>% add_x_position(x = "group") %>%
    left_join(maxmfi) %>%
    filter(group1 == "0_unstim")
  
  
  stat.test[stat.test$group2 == "2_GL","y.position"] <- stat.test[stat.test$group2 == "2_GL","y.position"] * 1.2

  set.seed(42)
  mfi_plt_human <- ggplot(mfi_human,
                          mapping = aes(x = group,
                                        y = mfi#,
                                        #color = as.character(day)
                          )
  ) +
    geom_violin(scale = "width", trim = F, fill = "lightblue") +
    geom_vline(xintercept = 1.5, linetype = "dashed", size = 0.3) +
    geom_jitter(size = 1, shape = 21, fill = "white", stroke = 0.25) +
    # stat_compare_means(method = "wilcox.test",
    #                    comparisons = comparisons,
    #                    exact = FALSE,
    #                    label = 'p.signif',
    #                    vjust = 0.5,
    #                    tip.length = 0, size = 1.5
    # ) +
    stat_pvalue_manual(stat.test,
                       vjust = 0.2,
                       tip.length = 0, size = 1.5) +
    theme_rgb() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          strip.background = element_blank(),
          axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(angle = 45, vjust = 0, hjust = 1)) +
    facet_wrap(~marker,
               nrow = 1,
               scales = "free_y") +
    xlab("Condition") +
    ylab("MFI") +
    ggtitle("Human")
  
  return(mfi_plt_human)
}

build_fig_s8_2 <- function() {
  mfi_mouse <- readRDS("data/processed/data_fig_6c_2.rds")
  
  mfi_mouse <- mfi_mouse %>%
    filter (day %in% c(0,2),
            condition %in% "unstim",
            marker %in% c("CD14", "CD40", "CD69", "IL4R", "PDL1")) %>%
    mutate(group  = factor(paste(day, condition, sep = "_"), levels = c("0_unstim",
                                                                        "2_unstim")))
  
  comparisons <- list(
    c('0_unstim', '2_unstim')
  )
  
  mfi_mouse %>% group_by(group) %>% shapiro_test(mfi) #sig = non-normality
  
  maxmfi <- mfi_mouse %>%
    group_by(marker) %>%
    summarise(y.position = max(mfi))
  
  kruskal.res <- mfi_mouse %>%
    filter(day %in% c(0,2)) %>%
    group_by(marker) %>%
    kruskal_test(mfi ~ group)
  
  # stopifnot(!any(kruskal.res$p > 0.05))
  
  stat.test <- mfi_mouse %>%
    filter(day %in% c(0,2)) %>%
    group_by(marker) %>%
    dunn_test(mfi ~ group, p.adjust.method = "holm") %>%
    add_significance() %>% add_x_position(x = "group") %>%
    left_join(maxmfi) %>%
    filter(group1 == "0_unstim")
  
  
  stat.test[stat.test$group2 == "2_unstim","y.position"] <- stat.test[stat.test$group2 == "2_unstim","y.position"] * 1.2
  # stat.test[stat.test$group2 == "2_GI","y.position"] <- stat.test[stat.test$group2 == "2_GI","y.position"] * 1.4
  # 
  
  set.seed(42)
  mfi_plt_mouse <- ggplot(mfi_mouse,
                          mapping = aes(x = group,
                                        y = mfi#,
                                        #color = as.character(day)
                          )
  ) +
    geom_violin(scale = "width", trim = F, fill = "lightblue") +
    geom_vline(xintercept = 1.5, linetype = "dashed", size = 0.3) +
    geom_jitter(size = 1, shape = 21, fill = "white", stroke = 0.25) +
    # stat_compare_means(method = "wilcox.test",
    #                    comparisons = comparisons,
    #                    exact = FALSE,
    #                    label = 'p.signif',
    #                    vjust = 0.5,
    #                    tip.length = 0, size = 1.5) +
    stat_pvalue_manual(stat.test,
                       vjust = 0.2,
                       tip.length = 0, size = 1.5) +
    #stat_pvalue_manual(stat.test, tip.length = 0) +
    theme_rgb() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          strip.background = element_blank(),
          axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(angle = 45, vjust = 0, hjust = 1)) +
    facet_wrap(~marker,
               nrow = 1,
               scales = "free_y") +
    xlab("Condition") +
    ylab("MFI") +
    ggtitle("Mouse")
  
  return (mfi_plt_mouse)
}

#####
# C #
#####

build_fig_6c_conditions <- function () {
  
  df_conditions <- readRDS("data/processed/data_fig_6d.rds") %>%
    mutate(species = case_when(
      group == "human" ~ "Hs",
      TRUE ~ "Mm"
    ))
  
  # plot_colors <- c("red", "grey70")
  # plot_colors <- c("#023e8a", "#caf0f8")
  # 
  # names(plot_colors) <- c("human", "mouse")
  # human_overlay <- ggplot(df_conditions, aes(x = x, y = y, colour = factor(df_conditions[,"group"]))) + 
  #   rasterize(geom_point(shape = 16, stroke = 0, size = 0.3), dpi = 300) + 
  #   scale_colour_manual(values = plot_colors) + 
  #   scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  #   scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  #   theme_classic() +
  #   labs(x = "DMAP 1",
  #        y = "DMAP. 2",
  #        color = "Condition",
  #        title = "Human") +
  #   theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  #         axis.line = element_blank(),
  #         legend.position = "none",
  #         plot.title = element_text(size = 8),
  #         axis.title = element_text(size = 8),
  #         axis.text = element_text(size = 6),
  #         aspect.ratio=1) +
  #   guides(color = guide_legend(override.aes = list(size = 3)))
  set.seed(42)
  p_conditions <- ggplot(df_conditions, aes(x = x, y = y, color = condition)) + 
    rasterize(geom_point(shape = 16, stroke = 0, size = 0.4), dpi = 300) +
    scale_color_manual(
      values = exp_cond_pal,
      breaks = names(exp_cond_pal),
      labels = c(
        "unstimulated",
        "GM-CSF + IFNy",
        "GM-CSF + LPS"
        )
    ) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
    theme_classic() +
    labs(x = "DMAP 1",
         y = "DMAP 2",
         color = "Experimental condition") +
    theme(panel.border = element_rect(colour = "black",
                                      fill=NA,
                                      size=0.5),
          axis.line = element_blank(),
          #legend.position = c(.75, .18),
          plot.title = element_text(size = 8),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          legend.title = element_text(face = "bold", size = 8),
          legend.text = element_text(size = 6),
          legend.background=element_rect(fill = alpha("white", 0)),
          legend.key=element_rect(fill = alpha("white", 0), color = alpha("white", 0)),
          legend.key.size = unit(0.2, "cm"),
          aspect.ratio=1) +
    guides(color = guide_legend(
      override.aes = list(size = 3,
                          fill = NA),
      title.position = "left"
      ))
  
  return (p_conditions)
}



build_fig_6c_species <- function () {
  
  df_conditions <- readRDS("data/processed/data_fig_6d.rds") %>%
    mutate(species = case_when(
      group == "human" ~ "Hs",
      TRUE ~ "Mm"
    ))
  
  # plot_colors <- c("red", "grey70")
  # plot_colors <- c("#023e8a", "#caf0f8")
  # 
  # names(plot_colors) <- c("human", "mouse")
  # human_overlay <- ggplot(df_conditions, aes(x = x, y = y, colour = factor(df_conditions[,"group"]))) + 
  #   rasterize(geom_point(shape = 16, stroke = 0, size = 0.3), dpi = 300) + 
  #   scale_colour_manual(values = plot_colors) + 
  #   scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  #   scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  #   theme_classic() +
  #   labs(x = "DMAP 1",
  #        y = "DMAP. 2",
  #        color = "Condition",
  #        title = "Human") +
  #   theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  #         axis.line = element_blank(),
  #         legend.position = "none",
  #         plot.title = element_text(size = 8),
  #         axis.title = element_text(size = 8),
  #         axis.text = element_text(size = 6),
  #         aspect.ratio=1) +
  #   guides(color = guide_legend(override.aes = list(size = 3)))
  
  # names(plot_colors) <- c("mouse", "human")
  set.seed(42)
  p_species <- ggplot(df_conditions,
                      aes(x = x,
                          y = y,
                          colour = species)) + 
    rasterize(geom_point(shape = 16, stroke = 0, size = 0.3), dpi = 300) + 
    scale_colour_manual(
      values = species_pal,
      breaks = names(species_pal),
      labels = c("Human", "Mouse")
    ) + 
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
    theme_classic() +
    labs(x = "DMAP 1",
         y = "DMAP 2",
         color = "Species") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          # legend.position = "none",
          plot.title = element_text(size = 8),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          legend.title = element_text(face = "bold", size = 8),
          legend.text = element_text(size = 6),
          aspect.ratio=1,
          legend.key.size = unit(0.2, "cm")) +
    guides(color = guide_legend(
      override.aes = aes(size = 2.5),
      title.position = "left"
      ))
  
  return (p_species)
}


# build_fig_6d_2 <- function () {
#   ### CONDITION UMAP ###
#   df_conditions_mouse <- readRDS("data/processed/data_fig_6d_2.rds")
#   
#   p_conditions <- ggplot(df_conditions_mouse, aes(x = x, y = y, color = condition)) + 
#     rasterize(geom_point(shape = 16, stroke = 0, size = 0.1), dpi = 300) +
#     #geom_density_2d() +
#     scale_color_brewer(palette = "Set1",
#                        labels = c("GM-CSF + IFN", "GM-CSF + LPS", "unstimulated")) +
#     theme_classic() +
#     labs(x = "UMAP1",
#          y = "UMAP2",
#          color = "Condition",
#          title = "Mus musculus, no CD40") +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#           axis.line = element_blank(),
#           legend.position = "left") +
#     guides(color = guide_legend(override.aes = list(size = 3)))
#   
#   # done in builder/notebooks, as ggExtra causes trouble with ggplot
#   # p_conditions <- ggExtra::ggMarginal(p_conditions, groupColour = TRUE, groupFill = TRUE) %>%
#   #   as_ggplot()
#   
#   return (p_conditions)
# }

#####
# D #
#####
build_fig_6d <- function() {
  #### MARKER UMAP ####
  df_markers <- readRDS("data/processed/data_fig_6e.rds")
  
  p_markers <- ggplot(df_markers %>%
           dplyr::filter(variable != "CD40") %>%
           mutate(group = c(human = "Human",
                            mouse = "Mouse")[group]), aes(x = x, y = y, color = value)) + 
    rasterize(geom_point(shape = 16, stroke = 0, size = 0.1), dpi = 300) +
    # scale_color_viridis(option = "inferno",
    #                     guide = guide_colorbar(frame.colour = "black",
    #                                            frame.linewidth = 0.5,
    #                                            draw.ulim = F,
    #                                            draw.llim = F)) +
    # scale_color_gradientn(colors = ice_blue_pal,
    #                       guide = guide_colorbar(
    #                         direction = "horizontal",
    #                         frame.colour = "black",
    #                         ticks.colour = "black",
    #                         frame.linewidth = 1,
    #                         title.hjust = -1
    #                       ))+
    scale_color_gradientn(colors = viridis(n = 25),
                          guide = guide_colorbar(
                            direction = "horizontal",
                            frame.colour = "black",
                            ticks.colour = "black",
                            frame.linewidth = 1,
                            title.hjust = -1
                          ))+
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
    theme_classic() +
    labs(x = "DMAP 1",
         y = "DMAP 2",
         color = "Scaled Expression") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          plot.title = element_text(size = 8),
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          # legend.position = "right",
          legend.title = element_text(face = "bold", size = 8),
          legend.text = element_text(size = 6),
          legend.key.size = unit(4, "mm"),
          aspect.ratio = 1) +
    facet_grid(group~variable)
  
  return (p_markers)
}

# build_fig_6e_2 <- function() {
#   #### MARKER UMAP ####
#   df_markers_mouse <- readRDS("data/processed/data_fig_6e_2.rds")
#   
#   p_markers_mouse <- ggplot(df_markers_mouse %>% dplyr::filter(variable != "CD40"), aes(x = x, y = y, color = value)) + 
#     rasterize(geom_point(shape = 16, stroke = 0, size = 0.1), dpi = 300) +
#     scale_color_viridis() +
#     theme_classic() +
#     labs(x = "UMAP1",
#          y = "UMAP2",
#          color = "Expr",
#          title = "Mus musculus, no CD40") +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#           axis.line = element_blank(),
#           legend.position = "right") +
#     facet_wrap(~variable)
#   
#   return (p_markers_mouse)
# }


############
# FIGURE 7 #
############
############

build_fig_7a <- function(){
  
  df <- read.csv("data/processed/MFI_scaled_pca.csv")
  
  scale_color_option = "viridis"
  
  organ <- ggplot(data = df, aes(x = PC1, y = PC2, colour = organ)) +
    rasterize(geom_point(shape = 16, stroke = 0, size = 1), dpi = 300) +
    scale_color_manual(
      values = c(PB = "#E41A1C", BM = "#377EB8", SPL = "#4DAF4A"),
      labels = c(
        BM = "Bone Marrow",
        PB = "Peripheral Blood",
        SPL = "Spleen"
      ),
      name = ""
    ) +
    theme_classic() +
    labs(x = "PC1",
         y = "PC2",
         title = "Organ") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.text=element_text(size = 6),
          axis.title=element_text(size = 6, face="bold"),
          plot.margin = unit(c(0,0,0,0), "cm"),
          aspect.ratio = 1,
          legend.text=element_text(size = 6),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold")) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  condition <- ggplot(data = df, aes(x = PC1, y = PC2, colour = condition)) +
    rasterize(geom_point(shape = 16, stroke = 0, size = 1), dpi = 300) +
    scale_color_manual(
      values = exp_cond_long_pal,
      breaks = names(exp_cond_long_pal),
      labels = c(
        "Unstimulated",
        "GM-CSF + IFNy",
        "GM-CSF + LPS"
      ),
      name = ""
    ) +    
    theme_classic() +
    labs(x = "PC1",
         y = "PC2",
         title = "Experimental condition") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          axis.text=element_text(size = 6),
          axis.title=element_text(size = 6, face="bold"),
          aspect.ratio = 1,
          legend.text=element_text(size = 6),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold")) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  PDL1 <- ggplot(data = df, aes(x = PC1, y = PC2, colour = PDL1)) +
    rasterize(geom_point(shape = 16, stroke = 0, size = 1), dpi = 300) +
    scale_color_viridis(option = scale_color_option) +
    theme_classic() +
    labs(x = "PC1",
         y = "PC2",
         title = "PDL1") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.text=element_text(size = 6),
          axis.title=element_text(size = 6, face="bold"),
          plot.margin = unit(c(0,0,0,0), "cm"),
          aspect.ratio = 1,
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold"))
  
  CD40 <- ggplot(data = df, aes(x = PC1, y = PC2, colour = CD40)) +
    rasterize(geom_point(shape = 16, stroke = 0, size = 1), dpi = 300) +
    scale_color_viridis(option = scale_color_option) +
    theme_classic() +
    labs(x = "PC1",
         y = "PC2",
         title = "CD40") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.text=element_text(size = 6),
          axis.title=element_text(size = 6, face="bold"),
          plot.margin = unit(c(0,0,0,0), "cm"),
          aspect.ratio = 1,
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold"))
  
  CD69 <- ggplot(data = df, aes(x = PC1, y = PC2, colour = CD69)) +
    rasterize(geom_point(shape = 16, stroke = 0, size = 1), dpi = 300) +
    scale_color_viridis(option = scale_color_option) +
    theme_classic() +
    labs(x = "PC1",
         y = "PC2",
         title = "CD69") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.text=element_text(size = 6),
          axis.title=element_text(size = 6, face="bold"),
          plot.margin = unit(c(0,0,0,0), "cm"),
          aspect.ratio = 1,
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold"))
  
  IL4R <- ggplot(data = df, aes(x = PC1, y = PC2, colour = IL4R)) +
    rasterize(geom_point(shape = 16, stroke = 0, size = 1), dpi = 300) +
    scale_color_viridis(option = scale_color_option,
                        name = "Scaled\nExpression",
                        labels = c(0, 0.25, 0.5, 0.75, 1)) +
    theme_classic() +
    labs(x = "PC1",
         y = "PC2",
         title = "IL4R") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.text=element_text(size = 6),
          axis.title=element_text(size = 6, face="bold"),
          plot.margin = unit(c(0,0,0,0), "cm"),
          aspect.ratio = 1,
          #legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
          legend.box.margin = margin(0,0,0,0),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6)) + 
    guides(color = guide_colorbar(barwidth = 0.5,
                                  barheight = 3,
                                  label = TRUE))
  
  CD14 <- ggplot(data = df, aes(x = PC1, y = PC2, colour = CD14)) +
    rasterize(geom_point(shape = 16, stroke = 0, size = 1), dpi = 300) +
    scale_color_viridis(option = scale_color_option) +
    theme_classic() +
    labs(x = "PC1",
         y = "PC2",
         title = "CD14") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.text=element_text(size = 6),
          axis.title=element_text(size = 6, face="bold"),
          plot.margin = unit(c(0,0,0,0), "cm"),
          aspect.ratio = 1,
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold"))
  
  legend <- get_legend(IL4R)
  
  markers <- cowplot::plot_grid(
    CD69,
    CD40 + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank()
    ),
    PDL1 + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank()
    ),
    CD14 + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank()
    ),
    IL4R + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 legend.position = "none"
    ),
    
    scale = c(1,1,1,1,1),
    nrow = 1,
    align = "v",
    axis = "tb"
  )
  
  organ_condition <- cowplot::plot_grid(organ, condition, nrow = 1, align = "v")
  
  fig_7a <- grid.arrange(organ_condition,
                         markers,
                         legend,
                         heights = c(1, 1.2),
                         widths = c(1, 0.1),
                         ncol = 2,
                         nrow = 2,
                         layout_matrix = rbind(c(1, 1),
                                               c(2, 3))
  )
  
  
  
  return(fig_7a)
}


build_fig_7b <- function(){
  
  bm_df <- readRDS("data/processed/data_fig_7b_bm.rds")
  spl_df <- readRDS("data/processed/data_fig_7b_spl.rds")
  pb_df <- readRDS("data/processed/data_fig_7b_pb.rds")
  full_data <- rbind(bm_df, pb_df, spl_df)
  
  organ_annotations <- c(BM = "Bone Marrow",
                         PB = "Peripheral Blood",
                         SPL = "Spleen")
  set.seed(42)
  fig_7b <- ggplot(full_data, aes(x = x, y = y, colour = condition)) + 
    rasterize(geom_point(shape = 16, stroke = 0, size = 0.4), dpi = 300) +  
    scale_color_manual(
      values = exp_cond_long_pal,
      breaks = names(exp_cond_long_pal),
      labels = c(
        "Unstimulated",
        "GM-CSF + IFNy",
        "GM-CSF + LPS"
      )
    ) +    
    theme_classic() +
    labs(x = "Diff. comp. 1",
         y = "Diff. comp. 2",
         color = NULL) +
    scale_x_continuous(breaks = seq(0, 1, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.text=element_text(size = 6),
          axis.title=element_text(size = 6, face="bold"),
          plot.margin = unit(c(0,0,0,0), "cm"),
          aspect.ratio = 1,
          legend.text=element_text(size = 6),
          legend.position = c(.10, .25),
          legend.key = element_blank(),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0,0,0,0)) +
    guides(color = guide_legend(override.aes = list(size = 3)),
           legend.spacing.y = unit(0, 'cm'),
           legend.margin = margin(0,0,0,0),
           legend.box.margin = margin(0,0,0,0)) +
    facet_grid(~organ,
               scales = "free",
               labeller = labeller(organ = organ_annotations))
  
  
  return(fig_7b)
}

build_fig_7c <- function(){
  
  
  bm_data <- readRDS("data/processed/data_fig_7c_bm.rds")
  spl_data <- readRDS("data/processed/data_fig_7c_spl.rds")
  pb_data <- readRDS("data/processed/data_fig_7c_pb.rds")
  
  full_data <- rbind(bm_data, pb_data, spl_data)
  
  organ_annotations <- c(BM = "Bone Marrow",
                         PB = "Peripheral Blood",
                         SPL = "Spleen")
  
  fig_7c <- ggplot(data = full_data, aes(x = x, y = y, color = value)) + 
    rasterize(geom_point(shape = 16, stroke = 0, size = 0.05), dpi = 300) +
    scale_color_viridis() +
    theme_classic() +
    labs(x = "Diff. comp. 1",
         y = "Diff. comp. 2",
         color = "Scaled Expr.") +
    scale_x_continuous(breaks = seq(0, 1, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          axis.text=element_text(size = 6),
          axis.title=element_text(size = 6, face="bold"),
          legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "cm"),
          aspect.ratio = 1,
          panel.spacing.x = unit(4, "mm")) +
    facet_grid(organ ~ variable,
               scales = "fixed",
               labeller = labeller(organ = organ_annotations))
  return (fig_7c)
}


###########################
###########################
## SUPPLEMENTARY FIGURES ##
###########################
###########################


#############
# FIGURE S1 #
#############
#############

build_fig_s1_1 <- function() {
  enrichr_go_bp <- readRDS("data/processed/addition_fig3c.rds")
  
  enrichr_go_bp_long <- Reduce(function(df1, df2) rbind(df1, df2), enrichr_go_bp)
  
  plt_data <- enrichr_go_bp_long %>%
    group_by(term_name) %>%
    # filter(length(term_name) >= 4) %>%
    mutate(mean_rank = mean(rank)) %>%
    group_by(query_description) %>%
    slice_min(order_by = mean_rank, n = 10, with_ties = F) %>%
    ungroup() %>%
    filter(adj_p_value < 0.05) %>%
    mutate(overlap = as.numeric(lapply(overlapping_genes, length))) %>%
    group_by(term_name) %>%
    mutate(mean_term_overlap = mean(overlap)) %>%
    ungroup()
  
  plt <- ggplot(plt_data, aes(x = query_description,
                              y = reorder(term_name, mean_term_overlap),
                              color = -log10(adj_p_value),
                              size = overlap)) +
    geom_point() +
    xlab("Gene set") +
    ylab("GO term") +
    # illegal use of whitespace incoming
    scale_size(name = "                             Overlap",
               guide = guide_legend(title.position = "left",
                                    title.hjust = 1)) +
    scale_color_continuous(name = expression("-Log"[10]*"(adjusted"~italic(P)*"-value)"),
                           guide = guide_colorbar(frame.colour = "black",
                                                  frame.linewidth = 1,
                                                  title.position = "left",
                                                  title.hjust = 1)) +
    theme_minimal() +
    theme(axis.title = element_text(size = 8),
          axis.text.x = element_text(size = 6,
                                     angle = 45,
                                     vjust = 1,
                                     hjust = 1),
          axis.text.y = element_text(size = 6),
          panel.border = element_rect(fill = NA),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8, angle = 90),
          legend.justification = "left",
          legend.direction = "vertical",
          legend.box.just = "right",
          legend.box = "horizontal",
          legend.position = "bottom"
    )
  return (plt)
}

build_fig_s1_2 <- function() {
  enrichr_go_bp <- readRDS("data/processed/addition_fig3c.rds")
  
  enrichr_go_bp_long <- Reduce(function(df1, df2) rbind(df1, df2), enrichr_go_bp)
  
  plt_data <- enrichr_go_bp_long %>%
    filter(adj_p_value <= 0.01) %>%
    mutate(overlap = as.numeric(lapply(overlapping_genes, length))) %>%
    arrange(rank)
  
  plt <- ggplot(plt_data, aes(x = -log10(adj_p_value),
                              y = rank,
                              color = -log10(adj_p_value),
                              size = overlap)) +
    geom_point() +
    xlab("Overlap") +
    ylab("GO term") +
    scale_color_continuous(name = "-log(adjusted P Value)") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA)) +
    facet_wrap(~query_description,
               scales = "free_y",
               ncol = 1) +
    geom_text_repel(aes(label = ifelse(rank <= 3, term_name, ""),
                        size = NULL,
                        color = NULL),
                    # box.padding = 0.3,
                    max.iter = 1000000,
                    max.time = 50)
  
  return (plt)
}

build_fig_s1_3 <- function () {
  enrichr_go_bp <- readRDS("data/processed/addition_fig3c.rds")
  
  enrichr_go_bp_long <- Reduce(function(df1, df2) rbind(df1, df2), enrichr_go_bp)
  
  plt_data <- enrichr_go_bp_long %>%
    filter(adj_p_value <= 0.05, rank <= 10) %>%
    arrange(query_description, rank) %>%
    mutate(overlap = as.numeric(lapply(overlapping_genes, length)),
           overall_rank = 1:length(term_name)) %>%
    dplyr::select(term_name, query_description, rank, overall_rank, adj_p_value, overlap)
  
  custom_y <- plt_data %>%
    pull(overall_rank, term_name)
  
  plt <- ggplot(plt_data, aes(x = -log10(adj_p_value),
                              y = as.factor(overall_rank),
                              color = -log10(adj_p_value),
                              size = overlap)) +
    geom_point() +
    xlab("-log(adjusted P Value)") +
    ylab("GO term") +
    scale_y_discrete(breaks = custom_y,
                     # labels = str_remove(names(custom_y), " \\(.*$")) +
                     labels = names(custom_y)) +
    scale_color_continuous(name = "-log(adjusted P Value)") +
    scale_size(name = "Overlap") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill = NA),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 6),
          strip.text = element_text(size = 8),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.justification = "right",
          legend.direction = "vertical",
          legend.box.just = "right",
          legend.box = "horizontal",
          legend.position = "bottom",
          aspect.ratio = 0.9
    ) +
    facet_wrap(~query_description,
               scales = "free_y",
               ncol = 1)
  plt
  return (plt)
}


######
# S? #
######

# build_fig_s00_1 <- function(){
#   goseq.results.df <- readRDS("data/processed/data_fig_4f.rds")
#   ind_Hs <- ifelse(str_sub(goseq.results.df$study,start=-2) == "Hs", TRUE, FALSE)
#   plt <- ggplot(data = goseq.results.df[goseq.results.df$category != "background" & ind_Hs,], aes(x = category, y = over_represented_p_score)) +
#     geom_bar(aes(fill = frac_de), stat="identity", color = "black") +
#     facet_wrap(~study, scale = "fixed", nrow=5) +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.8) +
#     scale_fill_gradient2(low = "white", high = "darkred",
#                          limits = c(0,1)) +
#     coord_flip() +
#     theme_rgb() +
#     labs(fill = "Fraction DE") +
#     theme(strip.text = element_text(size = 5.5), legend.title = element_text(face = "bold"))
#   return(plt)
# }
# 
# build_fig_s00_2 <- function(){
#   goseq.results.df <- readRDS("data/processed/data_fig_4f.rds")
#   ind_Hs <- ifelse(str_sub(goseq.results.df$study,start=-2) == "Hs", TRUE, FALSE)
#   plt <- ggplot(data = goseq.results.df[goseq.results.df$category != "background" & !ind_Hs,], aes(x = category, y = over_represented_p_score)) +
#     geom_bar(aes(fill = frac_de), stat="identity", color = "black") +
#     facet_wrap(~study, scale = "fixed", nrow=5) +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.8) +
#     scale_fill_gradient2(low = "white", high = "darkred",
#                          limits = c(0,1)) +
#     coord_flip() +
#     theme_rgb() +
#     labs(fill = "Fraction DE") +
#     theme(strip.text = element_text(size = 7), legend.title = element_text(face = "bold"))
#   return(plt)
# }


#############
# Figure S2 #
#############
#############

# with DE
build_fig_s2a <- function() {
  df <- readRDS("data/processed/data_fig_s2.rds")
  df <- df %>%
    arrange(desc(n_fisher)) %>%
    mutate(comparison = factor(comparison, levels = unique(comparison)))
  
  median_de <- median(df$n_de)
  median_up <- median(df$n_up)
  median_dn <- median(df$n_dn)
  
  df$n <- NULL
  df$n_minus_de <- NULL
  df$n_up <- NULL
  df$n_dn <- NULL
  df_melted <- melt(df)
  df_melted$variable <- factor(df_melted$variable, levels = c("n_de", "n_fisher"))
  
  df_melted_sum <- df_melted %>%
    group_by(comparison) %>%
    summarise(sum_value = sum(value))
  
  df_melted <- merge(df_melted, df_melted_sum, by = "comparison", all.x=T)
  
  
  label_median_de <- paste0(expression(tilde("Median: ")), median_de)
  
  plt <- ggplot(df_melted, aes(y = comparison, fill = variable)) +
    # geom_bar(aes(x = value), stat = "identity", position = "fill") +
    geom_bar(aes(x = value), stat = "identity") +
    geom_text(data = df_melted[df_melted$variable == "n_fisher",], aes(label=value, x = -30, angle = 270), vjust=0.5, color = "red", size = 1.5) +
    geom_text(data = df_melted[df_melted$variable == "n_de",], aes(label=value, x = sum_value+30, angle = 270), vjust=0.5, color = "grey70", size = 1.5) +
    ggplot2::annotate(size = 1.5, "text", hjust = 0, x = 4000, y = last(levels(df_melted$comparison)), label = dplyr::expr(paste(tilde("N"), ""[DE], " = ", !!median_de,sep=""))) +
    ggplot2::annotate(size = 1.5, "text", hjust = 0, x = 4000, y = levels(df_melted$comparison)[length(levels(df_melted$comparison))-1], label = expr(paste(tilde("N"), ""[UP], " = ", !!median_up,sep=""))) +
    ggplot2::annotate(size = 1.5, "text", hjust = 0, x = 4000, y = levels(df_melted$comparison)[length(levels(df_melted$comparison))-2], label = expr(paste(tilde("N"), ""[DOWN], " = ", !!median_dn,sep=""))) +
    scale_fill_manual("Subgroups of genes", values = c("grey70", "red"), labels = c("DE, not conserved", "DE + conserved")) +
    # scale_x_break(c(1750, 2500))+
    theme_rgb() +
    theme(axis.title.x = element_text(),
          axis.title.y = element_text()) +
    xlab("Number of genes") +
    ylab("Comparison")
  return(plt)
}

build_fig_s2b <- function() {
  lfc_df <- readRDS("data/processed/fisher_lfc_df.rds")
  p_lfc_df <- readRDS("data/processed/fisher_p_lfc_df.rds")
  up_fisher_human <- readRDS("data/processed/up_fisher_human.rds")
  
  validated_genes = c(
    "CD14",
    "CD40",
    "CD69",
    "CD274",
    "IL4R"
  )
  
  df <- lfc_df %>%
    pivot_longer(cols = ends_with(c("Hs", "Mm"))) %>%
    mutate(species = gsub("^.*_", "", name)) %>%
    group_by(species, symbol) %>%
    summarise(value = mean(value)) %>%
    pivot_wider(names_from = "species", values_from = "value") %>%
    left_join(p_lfc_df[,c("symbol", "fisher_adjusted")]) %>%
    arrange(desc(Hs+Mm)) %>%
    mutate(rank = 1:nrow(.)) %>%
    mutate(label = ifelse(rank <= 20 | rank >= (nrow(.)-19) | symbol %in% validated_genes, symbol, NA))
  
  df_core <- df %>%
    filter(symbol %in% up_fisher_human)
  
  set.seed(42)
  plt <- df %>%
    ggplot(aes(x = Hs, y = Mm)) +
    geom_hline(yintercept = 0, lwd = 0.25) +
    geom_vline(xintercept = 0, lwd = 0.25) +
    geom_point(pch = 19, size = 0.2, color = "grey70") +
    geom_point(data = df_core, color = "red", size = 0.2) +
    xlab("Mean log<sub>2</sub>(FC), Human") +
    ylab("Mean log<sub>2</sub>(FC), Mouse") +
    geom_text_repel(aes(label = label),
                    color = "black",
                    size = 1.5,
                    fontface = "italic",
                    min.segment.length = 0,
                    segment.size = 0.1,
                    max.overlaps = Inf) +
    scale_x_continuous(limits = symmetric_limits) +
    scale_y_continuous(limits = symmetric_limits) +
    theme_rgb() +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          legend.position = "none",
          aspect.ratio = 1)
  return(plt)
}

build_fig_s2c <- function() {
  study_stats <- readRDS("data/processed/data_fig_4a.rds")
  lfc_df <- readRDS("data/processed/fisher_lfc_df.rds")
  heatmap_matrix <- readRDS("data/processed/data_fig_4c_1.rds")
  
  
  merge_anno <- rownames(heatmap_matrix) %>%
    str_split(., pattern = "_") %>%
    sapply(., "[[",1) %>%
    data.frame() %>%
    dplyr::rename("study_id" = ".") %>%
    left_join(study_stats[,c("study_id", "model", "species")]) %>%
    mutate(orig = rownames(heatmap_matrix))
  
  merge_anno_sum <- merge_anno %>%
    filter(!duplicated(study_id)) %>%
    dplyr::select(-orig)
  
  df <- lfc_df %>%
    pivot_longer(cols = ends_with(c("Hs", "Mm"))) %>%
    mutate(study_id = stringr::str_split_fixed(name, '_', n=Inf)[, 3]) %>%
    left_join(merge_anno_sum) %>%
    filter(symbol %in% c("CD101", "CXCR4"))
  
    stat.test <- df %>%
      group_by(symbol) %>%
      wilcox_test(value ~ model) %>%
      adjust_pvalue(method = "holm") %>%
      add_significance() %>% add_xy_position(x = "model")
  
  


    set.seed(42)
  plt <- ggplot(df, aes(x = model, y = value)) +
    geom_boxplot(outlier.shape = NA, fill = "lightblue") +
    geom_jitter(shape = 21, alpha = 0.9, size = 1, fill = "white", stroke = 0.25) +
    facet_wrap(~symbol, scales = "free_y") +
        stat_pvalue_manual(stat.test, size = 2) +
    theme_rgb() +
    xlab("Model") +
    ylab("Mean log<sub>2</sub>(FC)") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          strip.text = element_text(face = "italic"),
          aspect.ratio = 1)
    
  return(plt)
}

# with all genes
build_fig_s2_2 <- function() {
  df <- readRDS("data/processed/data_fig_s2.rds")
  df <- df %>%
    arrange(desc(n_fisher)) %>%
    mutate(comparison = factor(comparison, levels = unique(comparison)))
  df$n <- NULL
  df_melted <- melt(df)
  df_melted$variable <- factor(df_melted$variable, levels = c("n_minus_de", "n_de", "n_fisher"))
  plt <- ggplot(df_melted, aes(y = comparison, fill = variable)) +
    # geom_bar(aes(x = value), stat = "identity", position = "fill") +
    geom_bar(aes(x = value), stat = "identity") +
    scale_fill_manual("Subgroups of genes", values = c("lightblue", "grey", "red"), labels = c("Not DE", "DE", "DE and conserved")) +
    theme_rgb() +
    theme(axis.title.x = element_text(),
          axis.title.y = element_text()) +
    xlab("Number of genes") +
    ylab("Comparison")
  return(plt)
}


#############
# FIGURE S3 #
#############
#############

build_fig_s3a <- function () {
  
  load("data/processed/LMM_res_composite.rda")
  
  fisher_core_genes <- readRDS("data/processed/up_fisher_human.rds")
  
  validated_genes <- c("CD14",
                       "CD40",
                       "CD69",
                       "CD274",
                       "IL4R")
  
  nrow(model_res_df[model_res_df$Beta >= 1 & model_res_df$adj_Pvalue <= 0.05,])
  nrow(model_res_df[model_res_df$Beta <= -1 & model_res_df$adj_Pvalue <= 0.05,])
  
  nrow(model_res_df[model_res_df$Beta >= 1 & model_res_df$adj_Pvalue <= 0.05 & model_res_df$symbol %in% fisher_core_genes,])

  
  
  plt_df <- model_res_df %>%
    mutate(color = ifelse(symbol %in% fisher_core_genes, "fisher", color),
           colfac = factor(color, levels = c("fisher", "both", "lfc", "p", "none")),
           label = case_when(
             -log10(adj_Pvalue)*Beta >=
               sort(-log10(adj_Pvalue)*Beta, decreasing = T)[20] ~ symbol,
             -log10(adj_Pvalue)*-Beta >=
               sort(-log10(adj_Pvalue)*-Beta, decreasing = T)[15] ~ symbol,
             symbol %in% validated_genes ~ symbol,
             TRUE ~ ""
           )) %>%
    arrange(desc(colfac))
  
  set.seed(42)
  ggplot(plt_df, aes(x = Beta, y = -log10(adj_Pvalue), color = color)) +
    geom_point(size = 0.75) +
    scale_color_manual(name = "Group",
                       breaks = c("fisher", "both", "lfc", "p", "none"),
                       labels = c("Core inflammation<br>program", "Log<sub>2</sub>(FC) and *P*", "Only log<sub>2</sub>(FC)", "Only *P*", "NS"),
                       values = c("orange", "red", "darkgreen", "blue", "black")) +
    geom_text_repel(aes(label = label, color = NULL, alpha = NULL),
                    #box.padding = 0.5,
                    max.overlaps = Inf,
                    show.legend = F,
                    size = 2,
                    seed = 42,
                    # nudge_x = 0.2,
                    min.segment.length = 0.1,
                    force_pull = 4,
                    fontface = "italic") +
    xlab(~beta) +
    ylab("-Log<sub>10</sub>(adjusted *P*-value)") +
    theme_rgb() +
    theme(legend.position = "right",
          legend.title = element_text(size = 8),
          legend.text = element_markdown(size = 6),
          axis.title.y = element_markdown(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text = element_text(size = 6))
  
 }

build_fig_s3b <- function () {
  lmm_individual_reproducibility <- readRDS("data/processed/lmm_individual_reproducibility.rds")
  
  plt_data <- lapply(lmm_individual_reproducibility, function(res){
    return(data.frame(comparison = res$comp, pi1 = res$pi1))
  }) %>% Reduce(function(df1, df2) rbind(df1, df2), .)
  
  set.seed(42)
  ggplot(plt_data, aes(x = "DESeq2 hits", y = pi1)) +
    geom_boxplot() +
    geom_jitter() +
    xlab("Gene set") +
    ylab(~pi[1]~"-statistic") +
    theme_rgb() +
    theme(axis.title = element_text(size = 8),
          axis.text =  element_text(size = 6))
  
}


build_fig_s3c <- function () {
  
  # sig_fisher_events, lmm_fisher_q, lmm_fisher_pi1
  load("data/processed/lmm_fisher_reproducibility.rda")
  
  ggplot(sig_fisher_events, aes(x = fisher)) +
    geom_histogram() +
    xlab("Fisher *P*-value") +
    ylab("Event count") +#bquote(.(labNames[1]) ~ x^2)
    ggtitle(substitute(paste(pi[1], spaceholder),
                       list(spaceholder=paste0("-statistic = ",
                                               lmm_fisher_pi1,
                                               ", N = ",
                                               nrow(sig_fisher_events))))) +
    theme_rgb() +
    theme(axis.title.y = element_text(size = 8),
          axis.title.x = element_markdown(size = 8),
          axis.text =  element_text(size = 6),
          plot.title = element_text(size = 8))
  
}

#############
# FIGURE S4 #
#############
#############

build_fig_s4a <- function () {
  load("data/processed/DE_testing_replicability.rda")
  
  df <- replication_res %>%
    filter(events_from != hits_from) %>%
    mutate(events_spec = str_remove(events_from, "^.*_"),
           hits_spec = str_remove(hits_from, "^.*_"))
  
  df_median <- df %>%
    group_by(events_from) %>%
    mutate(median_pi1 = median(pi1))
  
  print(paste0("Median probabilities ranked from ", min(df_median$median_pi1), " to ", max(df_median$median_pi1)))
  
  
  # argh
  rename_spec <- function(variable, value) {
    return(c(
      "Hs" = "Human",
      "Mm" = "Mouse"
    )[value])
  }
  
  
  df %>%
    ggplot(aes(x = events_from, y = pi1, color = hits_spec)) +
    geom_boxplot(aes(color = NULL), outlier.shape = NA, show.legend = F) +
    geom_jitter(size = 1) +
    coord_flip() +
    ylim(0,1) +
    facet_grid(rows = vars(events_spec),
               scales = "free_y",
               space = "free_y",
               labeller = rename_spec) +
    xlab("Comparison for gene set selection") +
    ylab(~pi[1]-statistic) +
    scale_color_discrete(
      name = "Species",
      breaks = c("Hs", "Mm"),
      labels = c("Human", "Mouse")
      ) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6))
}

build_fig_s4b <- function() {
  res_df <- readRDS("data/processed/study_fisher_NES.rds")
  
  plt_data <-res_df %>%
    rownames_to_column("pathway") %>%
    pivot_longer(cols = -pathway,
                 names_to = "ranked_study",
                 values_to = "NES") %>%
    mutate(pathway_group = str_remove(pathway, "__.*$"),
           pathway_study = str_remove(pathway, "^.*__"),
           rank_species = str_extract(ranked_study, "(?<=HC_).*(?=_NES)"))
  
  ggplot(plt_data, aes(x = NES, y = pathway_study, color = rank_species)) +
    geom_boxplot(aes(color = NULL),
                 outlier.shape = NA, show.legend = F) +
    geom_jitter(size = 1) +
    facet_wrap(~ pathway_group, ncol = 2) +
    xlab("NES in ranked comparisons") +
    ylab("Comparison for gene set selection") +
    scale_color_discrete(
      name = "Species",
      breaks = c("Hs", "Mm"),
      labels = c("Human", "Mouse")
      ) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6))
  
}

build_fig_s4c <- function() {
  load("data/processed/DE_testing_replicability.rda")
  
  ggplot(lfc_res %>%
           filter(name_outer != name_inner) %>%
           mutate(padj = ifelse(padj < 0.05, padj, NA)),
         aes(x = name_outer, y = cor, color = -log10(padj), size = rel_overlap)) +
    geom_boxplot(aes(color = NULL, size = NULL),
                 outlier.shape = NA,
                 show.legend = F) +
    geom_jitter() +
    coord_flip() +
    xlab("Comparison for gene set selection") +
    ylab("Pearson of"~log[2]~"(FC)") +
    scale_size_continuous(name = "Size of\noverlap",
                          range = c(0, 3)) +
    scale_color_gradientn(colors = c(brewer.pal(9, "OrRd")[4:9]),
                          na.value = "grey70",
                          name = ~-log[10](italic(P)),
                          guide = guide_colorbar(
                            frame.colour = "black",
                            frame.linewidth = 1
                          )) +
    # scale_color_continuous(name = ~-log[10](P)) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.key.height = unit(3, "mm"),
          legend.key.width = unit(3, "mm"))
}

#############
# Figure S5 #
#############
#############

build_fig_s5a <- function() {
  lfc_means_tf <- readRDS("data/processed/data_fig_4e.rds")
  p_lfc_rank_stats <- readRDS("data/processed/p_lfc_rank_stats.rds")
  
  lfc_means_tf_stats <- lfc_means_tf %>%
    left_join(p_lfc_rank_stats) %>%
    # remove TFs not tested in meta
    filter(!is.na(fisher_rank))
  
  lfc_means_tf_stats %<>%
    arrange(desc(mean_Hs+mean_Mm)) %>%
    mutate(enrichment_rank = 1:nrow(.)) %>%
    mutate(label = ifelse(enrichment_rank <= 10 | enrichment_rank >= (nrow(.) -9), symbol, NA),
           logp_plot = ifelse(fisher_rank <= 500, p_fisher_log_fixed, 0)) %>%
    arrange(desc(fisher_rank))
  
  set.seed(42)
  plt <- lfc_means_tf_stats %>%
    ggplot(
      aes(
        x = mean_Hs,
        y = mean_Mm,
        fill = logp_plot,
        size = logp_plot,
        label = label
        )
      ) + 
    geom_hline(yintercept = 0, lwd = 0.1) +
    geom_vline(xintercept = 0, lwd = 0.1) +
    # geom_point(size = 0.5) +
    geom_point(color = "black", pch = 21, alpha = 0.9) +
    scale_fill_gradientn(colors = c("grey70",
                                    brewer.pal(9, "OrRd")[4:9]),
                         # values = rescale(c(0,-log10(0.05),-log10(0.0001))),
                         breaks = c(0,
                                    100,
                                    150,
                                    200,
                                    250),
                         labels = c("ns",
                                    "1e-100",
                                    "1e-150",
                                    "1e-200",
                                    "1e-250")) +
    scale_size_continuous(breaks = c(0,
                                     100,
                                     150,
                                     200,
                                     250),
                          labels = c("ns",
                                     "1e-100",
                                     "1e-150",
                                     "1e-200",
                                     "1e-250")) +
    # geom_text_repel(aes(label = label), color = "black", size = 2,
    #                 fontface = "italic") +
    geom_label_repel(show.legend = F,
                     aes(size = 2),
                     color = "darkred",
                     fill = "white",
                     label.padding = unit(0.1, "lines"),
                     max.overlaps = Inf,
                     seed = 42,
                     fontface = "italic",
                     min.segment.length = 0,
                     segment.size = 0.25) +
    stat_cor(
      cor.coef.name = "r",
      label.x.npc = "left",
      label.y.npc = "top",
      method = "pearson",
      aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~"),color = NULL),
      size = 3) +
    xlim(-3,3) +
    ylim(-3,3) +
    xlab("TF-encoding gene, mean log<sub>2</sub>(FC), human") +
    ylab("TF-encoding gene, mean log<sub>2</sub>(FC), mouse") +
    theme_rgb() +
    labs(fill = expression(italic(P)*"-value"),
         size = expression(italic(P)*"-value")) +
    guides(fill=guide_legend(),
           size = guide_legend()) +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          legend.position = "right",
          aspect.ratio = 1)
  return(plt)
}

build_fig_s5b <- function() {
  act_merged <- readRDS("data/processed/data_fig_4f.rds")
  p_lfc_rank_stats <- readRDS("data/processed/p_lfc_rank_stats.rds") %>%
    mutate(logp_plot = ifelse(fisher_rank <= 500, p_fisher_log_fixed, 0))
  
  df <- act_merged %>%
    group_by(species, symbol_hs) %>%
    summarise(mean_score = mean(score),
              mean_p = mean(p_value)) %>%
    ungroup() %>%
    pivot_wider(names_from = "species",
                values_from = c("mean_score", "mean_p")) %>%
    arrange(desc(mean_score_Hs+mean_score_Mm)) %>%
    mutate(rank = 1:nrow(.)) %>%
    mutate(label = ifelse(rank <= 10 | rank >= (nrow(.) -9), symbol_hs, NA)) %>%
    mutate(mean_p = (mean_p_Hs+mean_p_Mm)/2) %>%
    left_join(p_lfc_rank_stats %>%
                transmute(
                  symbol_hs = symbol,
                  logp_plot
                )) %>%
    filter(!is.na(logp_plot)) %>%
    arrange(logp_plot)
  
  plt <- df %>%
    ggplot(
      aes(
        x = mean_score_Hs,
        y = mean_score_Mm,
        fill = logp_plot,
        size = logp_plot,
        label = label
      )
      ) +
    geom_vline(xintercept = 0, lwd = 0.1) +
    geom_hline(yintercept = 0, lwd = 0.1) +
    # geom_point(size = 0.5) + 
    geom_point(color = "black", pch = 21, alpha = 0.9) +
    scale_fill_gradientn(colors = c("grey70",
                                    brewer.pal(9, "OrRd")[4:9]),
                         # values = rescale(c(0,-log10(0.05),-log10(0.0001))),
                         breaks = c(0,
                                    100,
                                    150,
                                    200,
                                    250),
                         labels = c("ns",
                                    "1e-100",
                                    "1e-150",
                                    "1e-200",
                                    "1e-250")) +
    scale_size_continuous(breaks = c(0,
                                     100,
                                     150,
                                     200,
                                     250),
                          labels = c("ns",
                                     "1e-100",
                                     "1e-150",
                                     "1e-200",
                                     "1e-250")) +
    # geom_text_repel(aes(label = label), color = "black", size = 2,
    #                 fontface = "italic") +
    geom_label_repel(
      show.legend = F,
      aes(size = 2),
      color = "darkred",
      fill = "white",
      label.padding = unit(0.1, "lines"),
      max.overlaps = Inf,
      seed = 42,
      fontface = "italic",
      min.segment.length = 0,
      segment.size = 0.25
      ) +
    stat_cor(
      cor.coef.name = "r",
             method = "pearson",
             label.x.npc = "left",
             label.y.npc = "top",
             aes(label = paste(..r.label.., gsub("p", "P", ..p.label..), sep = "~`,`~"),color = NULL),
             size = 3) +
    xlim(-2,5) +
    ylim(-2,5) +
    xlab("TF enrichment score, human") + 
    ylab("TF enrichment score, mouse") +
    theme_rgb() +
    labs(fill = expression(italic(P)*"-value"),
         size = expression(italic(P)*"-value")) +
    guides(fill=guide_legend(),
           size = guide_legend()) +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          legend.position = "right",
          aspect.ratio = 1)
  return(plt)
}

#############
# Figure S7 #
#############
#############



#############
# Figure S7 #
#############
#############

build_fig_s7a <- function() {
  
  fisher_core <- readRDS("data/processed/up_fisher_human.rds")
  set.seed(42)
  plt_meta <- readRDS("data/processed/data_1_fig_s6a.rds") %>%
    as_tibble() %>%
    transmute(
      SRX_accession,
      species = ifelse(species == "Hs", "Human", "Mouse"),
      condition = ifelse(condition == "HC", "HC", "INFL")
      ) %>%
    group_by(species, condition) %>%
    mutate(SRX_shuffled = sample(SRX_accession))
  plt <- readRDS("data/processed/data_2_fig_s6a.rds")
  gene_df <- readRDS("data/processed/data_3_fig_s6a.rds") %>%
    mutate(
      annotation_color = ifelse(MEs %in% c("ME33", "ME10"), "black", "darkgrey"),
      MEs = factor(
        x = str_replace(MEs, "^ME", "Mod "),
        levels = c("Mod 33", "Mod 34", "Mod 10", "Mod 21")
      ))
  
  plt <- plt[, plt_meta$SRX_shuffled]
  
  stopifnot(all.equal(colnames(plt), plt_meta$SRX_shuffled))
  stopifnot(all.equal(rownames(plt), gene_df$gene))
  
  topAnno <- columnAnnotation(
    Condition = plt_meta$condition,
    Species = plt_meta$species,
    col = list(
      Condition = hc_infl_pal,
      Species = species_pal_en
    ),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(border = TRUE),
    border = T
  )
  rightAnno <- rowAnnotation(
    # Gene = anno_mark(
    #   at = c(match(fisher_core, rownames(plt))),
    #   labels = fisher_core,
    #   labels_gp = gpar(fontsize = 6)
    # ),
    infisher = anno_simple(
      as.character(rownames(plt) %in% fisher_core),
      col = c(`TRUE` = "black", `FALSE` = "white"),
      border = TRUE
    ),
    genenames = anno_mark(
      at = which(rownames(plt) %in% fisher_core),
      labels = rownames(plt)[rownames(plt) %in% fisher_core],
      padding = 0.1,
      labels_gp = gpar(
        col = gene_df$annotation_color[rownames(plt) %in% fisher_core],
        fontface = "italic",
        fontsize = 3
        ),
      link_gp = gpar(
        col = gene_df$annotation_color[rownames(plt) %in% fisher_core],
        lwd = 0.2
        )
    ),
    show_annotation_name = c(infisher = FALSE, genenames = FALSE)
  )
  
  colfun <- generate_rd_bu_colfun(min = quantile(plt, 0.01), quantile(plt, 0.99))
  
  Heatmap(plt,
          top_annotation = topAnno,
          right_annotation = rightAnno,
          show_row_names = F,
          show_column_names = F,
          cluster_columns = F,
          show_row_dend = F,
          col = colfun,
          column_split = plt_meta$condition,
          column_title = c(),
          # name = "z-score",
          row_split = gene_df$MEs,
          cluster_row_slices = F,
          row_title_gp = gpar(fontsize = 8),
          border = TRUE,
          heatmap_legend_param = list(
            # labels_gp = gpar(fontsize = 6),
            # title_gp = gpar(fontsize = 8, fontface = "bold"),
            title_position = "lefttop-rot",
            title = expression(bold(paste(bolditalic(z), "-score", sep=""))),
            border = TRUE
          )
  )
}


build_fig_s7b <- function () {
  
  plt_mod_betas <- readRDS("data/processed/data_1_fig_s6b.rds") %>%
    dplyr::select(-`(Intercept)`)
  module_stats <- readRDS("data/processed/data_2_fig_s6b.rds") %>%
    column_to_rownames("module_name")
  module_stats <- module_stats[rownames(plt_mod_betas), ] %>%
    mutate(padj = ifelse(padj < 0.05, padj, NA))
  
  stopifnot(all(rownames(plt_mod_betas) == rownames(module_stats)))
  
  rightAnno <-rowAnnotation(
    P = module_stats$padj,
    Size = module_stats$setsize,
    col = list(
      P = generate_or_rd_colfun(min(module_stats$padj, na.rm = T), 0.05),
      Size = generate_greens_colfun(min(module_stats$setsize), max(module_stats$setsize))
      ),
    na_col = "grey70",
    border = T,
    # annotation_name_gp = gpar(fontsize = 8, fontface = "italic"),
    annotation_legend_param = list(
      P = list(
        title = expression(bolditalic(P)),
        border = T
      ),
      Size = list(
        title = "Size",
        border = T
      )),
    show_annotation_name = c(F, F)
  )
  
  maxb <- max(
    abs(
      c(
        quantile(as.matrix(plt_mod_betas), 0.01),
        quantile(as.matrix(plt_mod_betas), 0.99)
        )
      )
  )
  
  colfun = generate_rd_bu_colfun(
    -maxb,
    maxb
  )
  
  Heatmap(plt_mod_betas,
          right_annotation = rightAnno,
          column_split = factor(str_extract_all(colnames(plt_mod_betas),
                                         "condition|species", simplify = T),
                                levels = c("species", "condition")),
          column_title = c("Species", "Condition"),
          column_labels = str_remove_all(colnames(plt_mod_betas), "species|condition_all"),
          row_labels = str_replace(rownames(plt_mod_betas), "^ME", ""),
          row_title = "Modules",
          row_title_side = "left",
          row_names_side = "left",
          col = colfun,
          # name = "beta",
          heatmap_legend_param = list(
            title = "Beta"
          ),
          clustering_method_columns = "ward.D2",
          clustering_method_rows = "ward.D2",
          border = T,
          column_title_gp = gpar(fontsize = 8),
          row_title_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 8),
          column_names_rot = 45
          )
}

# New S7
build_fig_S7_1 <- function() {
  mfi_human <- readRDS("data/processed/data_fig_6c_1.rds")
  mfi_human <- mfi_human %>%
    filter (day %in% c(0,2), condition == "unstim",
            marker %in% c("CD14", "CD40", "CD69", "IL4R", "PDL1")) %>%
    mutate(group  = factor(paste(day, condition, sep = "_"), levels = c("0_unstim",
                                                                        "2_unstim"
                                                                        )))
  
  comparisons <- list(
    c('0_unstim', '2_unstim')
  )
  
  
  
  mfi_human %>% group_by(group) %>% shapiro_test(mfi) #sig = non-normality
  
  maxmfi <- mfi_human %>%
    group_by(marker) %>%
    summarise(y.position = max(mfi))
  
  kruskal.res <- mfi_human %>%
    filter(day == 2) %>%
    group_by(marker) %>%
    kruskal_test(mfi ~ group)
  
  stopifnot(!any(kruskal.res$p > 0.05))
  
  stat.test <- mfi_human %>%
    filter(day == 2) %>%
    group_by(marker) %>%
    dunn_test(mfi ~ group, p.adjust.method = "holm") %>%
    add_significance() %>% add_x_position(x = "group") %>%
    left_join(maxmfi) %>%
    filter(group1 == "0_unstim")
  
  
  stat.test[stat.test$group2 == "2_unstim","y.position"] <- stat.test[stat.test$group2 == "2_unstim","y.position"] * 1.2
  stat.test[stat.test$group2 == "2_unstim","y.position"] <- stat.test[stat.test$group2 == "2_unstim","y.position"] * 1.4
  
  mfi_plt_human <- ggplot(mfi_human,
                          mapping = aes(x = group,
                                        y = mfi#,
                                        #color = as.character(day)
                          )
  ) +
    geom_violin(scale = "width", trim = F, fill = "lightblue") +
    geom_jitter(size = 1, shape = 21, fill = "white", stroke = 0.25) +
    # stat_compare_means(method = "wilcox.test",
    #                    comparisons = comparisons,
    #                    exact = FALSE,
    #                    label = 'p.signif',
    #                    vjust = 0.5,
    #                    tip.length = 0, size = 1.5
    # ) +
    stat_pvalue_manual(stat.test,
                       vjust = 0.5,
                       tip.length = 0, size = 1.5) +
    theme_rgb() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_blank(),
          strip.background = element_blank(),
          axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(angle = 45, vjust = 1, hjust = 0)) +
    facet_wrap(~marker,
               scales = "free_y") +
    xlab("Condition") +
    ylab("Median fluorescence intensity") +
    ggtitle("Homo sapiens")
  
  return(mfi_plt_human)
}

#############
# Figure S6 #
#############
#############

build_fig_s8a <- function() {
  
  p_lfc_df <- readRDS("data/processed/fisher_p_lfc_df.rds")
  fisher_up_genes <- readRDS("data/processed/fisher_up_genes.rds")
  fisher_dn_genes <- readRDS("data/processed/fisher_dn_genes.rds")
  validated_genes <- c(
    "CD14",
    "CD40",
    "CD69",
    "CD274",
    "IL4R"
  )
  
  plt_df <- p_lfc_df %>%
    arrange(fisher_adjusted) %>%
    mutate(rank = 1:length(fisher_adjusted),
           logfish = -log10(fisher_adjusted),
           group = case_when(
             symbol %in% validated_genes ~ "a_validated",
             rank <= 500 & !(symbol %in% validated_genes) ~ "p",
             TRUE ~ "none"
           ),
           label = case_when(
             symbol %in% validated_genes ~ symbol
           )) %>%
    arrange(desc(group))
  # blubb
  last_fish <- plt_df[plt_df$rank == 500, "fisher_adjusted"]
  
  ggplot (plt_df, aes (x = rank,
                       y = logfish,
                       color = group,
                       label = label)) +
    geom_vline(xintercept = 500, lwd = 0.25,
               linetype = "dashed") +
    geom_hline(yintercept = plt_df[plt_df$rank == 500, "logfish"],
               lwd = 0.25,
               linetype = "dashed") +
    annotate("text",
             x = 550,
             y = 300,
             label = "Rank = 500",
             size = 2,
             hjust = 0,
             vjust = 1) +
    annotate("text",
             x = max(plt_df$rank),
             y = plt_df[plt_df$rank == 500, "logfish"] + 5,
             label = list(bquote(italic(P)==.(last_fish))),
             size = 2,
             hjust = 1,
             vjust = 0,
             parse = T) +
    geom_point(size = 0.5) +
    scale_color_manual(name = "Group",
                       values = c(
                         "grey70",
                         brewer.pal(9, "OrRd")[3],
                         brewer.pal(9, "OrRd")[9]
                       ),
                       breaks = c(
                         "none",
                         "p",
                         "a_validated"
                       ),
                       labels = c(
                         "Below cutoff",
                         "Above cutoff",
                         "Validated"
                       )
                       ) +
    geom_label_repel(show.legend = F,
                     fontface = "italic",
                     nudge_x = 400,
                     min.segment.length = 0,
                     size = 2,
                     label.padding = unit(0.1, "lines"),
                     box.padding = unit(0.1, "lines")) +
    theme_rgb() +
    theme(axis.title = element_text(),
          legend.position = "bottom") +
    xlab("Rank") +
    ylab(expression("-Log"["10"]("adjusted Fisher "*italic(P)*"-value")))
  
}

build_fig_s8b <- function() {
  
  p_lfc_df <- readRDS("data/processed/fisher_p_lfc_df.rds")
  fisher_up_genes <- readRDS("data/processed/fisher_up_genes.rds")
  fisher_dn_genes <- readRDS("data/processed/fisher_dn_genes.rds")
  validated_genes <- c(
    "CD14",
    "CD40",
    "CD69",
    "CD274",
    "IL4R"
  )
  
  plt_df <- p_lfc_df %>%
    arrange(fisher_adjusted) %>%
    mutate(rank = 1:length(fisher_adjusted),
           logfish = -log10(fisher_adjusted),
           group = case_when(
             symbol %in% fisher_up_genes ~ "up",
             symbol %in% fisher_dn_genes ~ "dn",
             rank <= 500 & !(symbol %in% c(fisher_up_genes, fisher_dn_genes)) ~ "p",
             TRUE ~ "none"
           ),
           label = case_when(
             symbol %in% validated_genes ~ symbol
           ),
           masked_lfc = case_when(
             abs(mean_lfc) >= 0.5 ~ mean_lfc,
             TRUE ~ 0
           )) %>%
    arrange(desc(fisher_adjusted)) %>%
    filter(rank <= 500) %>%
    arrange(abs(mean_lfc))
  
  mean_lfc_vec <- plt_df$mean_lfc
  
  ggplot (plt_df, aes (x = rank,
                       y = logfish,
                       color = mean_lfc,
                       label = label)) +
    geom_point(size = 0.5) +
    geom_label_repel(show.legend = F,
                     min.segment.length = 0,
                     size = 2,
                     label.padding = unit(0.1, "lines"),
                     box.padding = unit(0.1, "lines"),
                     nudge_x = 10,
                     nudge_y = 10) +
    theme_rgb() +
    scale_color_gradientn(
      name = expression("Mean Log"["2"]*" Fold Change"),
      colors = c(
        brewer.pal(9, "RdBu")[9:7],
        "grey70",
        "grey70",
        brewer.pal(9, "RdBu")[3:1]
      ),
      values = c(
        rescale(c(
          seq(min(mean_lfc_vec), -0.5, abs(min(mean_lfc_vec)+0.5)/(3-1)),
          -0.4999,
          0.4999,
          seq(0.5, max(mean_lfc_vec), (max(mean_lfc_vec)-0.5)/(3-1))
          )
        )
      )
      ) +
    theme(legend.position = "bottom",
          axis.title = element_text()) +
    xlab("Rank") +
    ylab(expression("-log"["10"]("adjusted Fisher "*italic(P)*" value")))
    
  
}



build_fig_s2e <- function() {
  
  p_lfc_df <- readRDS("data/processed/fisher_p_lfc_df.rds")
  fisher_up_genes <- readRDS("data/processed/fisher_up_genes.rds")
  fisher_dn_genes <- readRDS("data/processed/fisher_dn_genes.rds")
  validated_genes <- c(
    "CD14",
    "CD40",
    "CD69",
    "CD274",
    "IL4R"
  )
  
  plt_df <- p_lfc_df %>%
    arrange(fisher_adjusted) %>%
    mutate(rank = 1:length(fisher_adjusted),
           logfish = -log10(fisher_adjusted),
           group = case_when(
             symbol %in% fisher_up_genes ~ "up",
             symbol %in% fisher_dn_genes ~ "dn",
             rank <= 500 & !(symbol %in% c(fisher_up_genes, fisher_dn_genes)) ~ "p",
             TRUE ~ "none"
           ),
           label = case_when(
             symbol %in% validated_genes ~ symbol
           ),
           masked_lfc = case_when(
             abs(mean_lfc) >= 0.5 ~ mean_lfc,
             TRUE ~ 0
           ),
           shape = ifelse(masked_lfc > 0, "16", NA)) %>%
    arrange(desc(fisher_adjusted)) %>%
    filter(rank <= 500) %>%
    arrange(abs(mean_lfc))
  
  mean_lfc_vec <- plt_df$mean_lfc

  ggplot (plt_df, aes (x = rank,
                       y = logfish,
                       color = mean_lfc,
                       label = label,
                       shape = shape)) +
    geom_point(size = 0.5) +
    scale_shape(guide = "none") +
    geom_label_repel(show.legend = F,
                     min.segment.length = 0,
                     fontface = "italic",
                     size = 2,
                     label.padding = unit(0.1, "lines"),
                     box.padding = unit(0.1, "lines"),
                     nudge_x = 10,
                     nudge_y = 10) +
    theme_rgb() +
    scale_color_gradientn(
      name = expression(bold("Mean log"["2"]*" (FC)")),
      colors = c(
        brewer.pal(9, "RdBu")[9:7],
        "grey70",
        "grey70",
        brewer.pal(9, "RdBu")[3:1]
      ),
      values = c(
        rescale(c(
          seq(min(mean_lfc_vec), -0.5, abs(min(mean_lfc_vec)+0.5)/(3-1)),
          -0.4999,
          0.4999,
          seq(0.5, max(mean_lfc_vec), (max(mean_lfc_vec)-0.5)/(3-1))
        )
        )
      ),
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        frame.linewidth = 1
      )
    ) +
    theme(legend.position = "bottom",
          axis.title = element_text()) +
    xlab("Rank") +
    ylab(expression("-Log"["10"]("adjusted Fisher "*italic(P)*" value"))) +
    ggtitle(paste0("Upregulated subset in top 500 genes, N = ", sum(plt_df$shape == "16", na.rm=T)))
  
  
}


build_fig_s2f <- function() {
  
  p_lfc_df <- readRDS("data/processed/fisher_p_lfc_df.rds")
  fisher_up_genes <- readRDS("data/processed/fisher_up_genes.rds")
  fisher_dn_genes <- readRDS("data/processed/fisher_dn_genes.rds")
  label_genes <- c(
    "CXCR4",
    "CD101"
  )
  
  plt_df <- p_lfc_df %>%
    arrange(fisher_adjusted) %>%
    mutate(rank = 1:length(fisher_adjusted),
           logfish = -log10(fisher_adjusted),
           group = case_when(
             symbol %in% fisher_up_genes ~ "up",
             symbol %in% fisher_dn_genes ~ "dn",
             rank <= 500 & !(symbol %in% c(fisher_up_genes, fisher_dn_genes)) ~ "p",
             TRUE ~ "none"
           ),
           label = case_when(
             symbol %in% label_genes ~ symbol
           ),
           masked_lfc = case_when(
             abs(mean_lfc) >= 0.5 ~ mean_lfc,
             TRUE ~ 0
           ),
           shape = ifelse(masked_lfc < 0, "16", NA)) %>%
    arrange(desc(fisher_adjusted)) %>%
    filter(rank <= 500) %>%
    arrange(abs(mean_lfc))
  
  mean_lfc_vec <- plt_df$mean_lfc
  
  
  ggplot (plt_df, aes (x = rank,
                       y = logfish,
                       color = mean_lfc,
                       label = label,
                       shape = shape)) +
    geom_point(size = 0.5) +
    scale_shape(guide = "none") +
    geom_label_repel(show.legend = F,
                     min.segment.length = 0,
                     size = 2,
                     fontface = "italic",
                     label.padding = unit(0.1, "lines"),
                     box.padding = unit(0.1, "lines"),
                     nudge_x = 10,
                     nudge_y = 10) +
    theme_rgb() +
    scale_color_gradientn(
      name = expression("Mean Log"["2"]*" Fold Change"),
      colors = c(
        brewer.pal(9, "RdBu")[9:7],
        "grey70",
        "grey70",
        brewer.pal(9, "RdBu")[3:1]
      ),
      values = c(
        rescale(c(
          seq(min(mean_lfc_vec), -0.5, abs(min(mean_lfc_vec)+0.5)/(3-1)),
          -0.4999,
          0.4999,
          seq(0.5, max(mean_lfc_vec), (max(mean_lfc_vec)-0.5)/(3-1))
        )
        )
      )
    ) +
    theme(legend.position = "bottom",
          axis.title = element_text()) +
    xlab("Rank") +
    ylab(expression("-Log"["10"]("adjusted Fisher "*italic(P)*" value"))) +
    ggtitle(paste0("Downregulated subset in top 500 genes, N = ", sum(plt_df$shape == "16", na.rm=T)))
  
  
  
}

##########################
##########################
## SUPPLEMENTARY TABLES ##
##########################
##########################
# Table builders

build_table_S5 <- function(){
  plt_ranking <- readRDS("data/processed/data_1_fig_1c.rds")
  colnames(plt_ranking) <- c("Lineage", "Symbol", "Mean count")
  write_xlsx(plt_ranking, "tables/Table S5.xlsx")
}
