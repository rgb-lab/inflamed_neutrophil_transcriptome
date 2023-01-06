

if (!require(marge)) devtools::install_github('robertamezquita/marge', ref = 'master')

# options('homer_path' = '/Users/fr929/homer')

read_homer_results <- function(comparison){
  current.path <- paste0("data/raw_data/zymosan/atac/", comparison, "/knownResults.txt")
  current.comparison <- read.delim(current.path, row.names = NULL,
                                  colClasses=c("P.value"="character"))
  
  colnames(current.comparison) <- c("motif_name",
                                    "consensus",
                                    "p_value",
                                    "log_p_value_score",
                                    "fdr",
                                    "tgt_num",
                                    "tgt_pct",
                                    "bgd_num",
                                    "bgd_pct")
  current.comparison %<>%
    separate("motif_name", c("motif_name", 
                             "experiment", "database"), "/", extra = "drop") %>%
    separate("motif_name", into = "motif_name", sep = '\\(', extra = "drop") %>%
    separate("p_value", into =c("base","log_p_value"), sep = "e", extra = "drop") %>%
    dplyr::select(-base) %>%
    mutate(log_p_value = abs(as.numeric(log_p_value)),
           log_p_value_score = abs(log_p_value_score),
           tgt_pct = as.numeric(sub("%", "", tgt_pct)),
           bgd_pct = as.numeric(sub("%", "", bgd_pct)))
  
  return(current.comparison)
}

BL_VS_BM_increase <- read_homer_results("BL_VS_BM_increase")
BL_VS_BM_decrease <- read_homer_results("BL_VS_BM_decrease")
MEM_VS_BM_increase <- read_homer_results("MEM_VS_BM_increase")
MEM_VS_BM_decrease <- read_homer_results("MEM_VS_BM_decrease")
AP_VS_BM_increase <- read_homer_results("AP_VS_BM_increase")
AP_VS_BM_decrease <- read_homer_results("AP_VS_BM_decrease")
MEM_VS_BL_increase <- read_homer_results("MEM_VS_BL_increase")
MEM_VS_BL_decrease <- read_homer_results("MEM_VS_BL_decrease")
AP_VS_BL_increase <- read_homer_results("AP_VS_BL_increase")
AP_VS_BL_decrease <- read_homer_results("AP_VS_BL_decrease")

                              

# Create TF - humanize-gene table. Input can be any HOMER result (all query the full motif list of interest)
ortho_mouse_to_human <- readRDS("data/processed/ortho_mouse_to_human.rds")
TF_human <- AP_VS_BL_increase %>%
left_join(ortho_mouse_to_human[,c("external_gene_name", "hsapiens_homolog_associated_gene_name")], by = c("motif_name" = "external_gene_name"))
TF_human <- TF_human[!duplicated(TF_human$motif_name),c("motif_name", "hsapiens_homolog_associated_gene_name")]
TF_human <- apply(TF_human,2,as.character)
write.table(TF_human, "motif_genes.txt", sep = "\t", col.names = NA)
# Read manually annotated TF_human table, merge it with the query of interest

join_TF_human_genes <- function(ME_dn, ME_up) {
  chea_res_dn_df_all <- readRDS("data/processed/data_fig_5a_2.rds")
  chea_res_up_df_all <- readRDS("data/processed/data_fig_5a_1.rds")
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
  
  
  TF_human <- read.delim("motif_genes_edited.txt", sep = "\t")
  
  TF_merged_dn <- merge(ME_dn, TF_human, by = "motif_name")
  TF_merged_up <- merge(ME_up, TF_human, by = "motif_name")
  
  TF_merged_dn %<>% 
    group_by(hsapiens_homolog_associated_gene_name) %>%
    summarise(log_p_value = mean(log_p_value),
              log_p_value_score = mean(log_p_value_score),
              fdr = mean(fdr),
              tgt_num = mean(tgt_num), 
              tgt_pct = mean(tgt_pct), 
              bgd_num = mean(bgd_num), 
              bgd_pct = mean(bgd_pct)) %>%
    arrange(desc(log_p_value)) %>%
    mutate(rank_motif = 1:nrow(.))
  TF_merged_up %<>% 
    group_by(hsapiens_homolog_associated_gene_name) %>%
    summarise(log_p_value = mean(log_p_value),
              log_p_value_score = mean(log_p_value_score),
              fdr = mean(fdr),
              tgt_num = mean(tgt_num), 
              tgt_pct = mean(tgt_pct), 
              bgd_num = mean(bgd_num), 
              bgd_pct = mean(bgd_pct)) %>%
    arrange(desc(log_p_value)) %>%
    mutate(rank_motif = 1:nrow(.))
  
  
  TF_merged_dn <- merge(TF_merged_dn, chea_res_dn_df_all, by.x = "hsapiens_homolog_associated_gene_name", by.y = "TF")
  TF_merged_up <- merge(TF_merged_up, chea_res_up_df_all, by.x = "hsapiens_homolog_associated_gene_name", by.y = "TF")
  
  TF_merged_dn %<>%
    mutate(mean_tf = (Hs+Mm)/2) %>%
    filter(!is.na(mean_tf))
  TF_merged_up %<>%
    mutate(mean_tf = (Hs+Mm)/2) %>%
    filter(!is.na(mean_tf))

  TF_merged_dn %<>%
    group_by(hsapiens_homolog_associated_gene_name) %>%
    summarise(log_p_value = mean(log_p_value),
              log_p_value_score = mean(log_p_value_score),
              fdr = mean(fdr),
              tgt_num = mean(tgt_num), 
              tgt_pct = mean(tgt_pct), 
              bgd_num = mean(bgd_num), 
              bgd_pct = mean(bgd_pct),
              mean_tf = mean(mean_tf),
              rank_motif = mean(rank_motif))
  TF_merged_up %<>%
    group_by(hsapiens_homolog_associated_gene_name) %>%
    summarise(log_p_value = mean(log_p_value),
              log_p_value_score = mean(log_p_value_score),
              fdr = mean(fdr),
              tgt_num = mean(tgt_num), 
              tgt_pct = mean(tgt_pct), 
              bgd_num = mean(bgd_num), 
              bgd_pct = mean(bgd_pct),
              mean_tf = mean(mean_tf),
              rank_motif = mean(rank_motif))
  
  volcano_stats <- readRDS("data/processed/data_fig_4b.rds")

  volcano_stats_up <- volcano_stats[volcano_stats$lfc >= 0,]
  volcano_stats_up <- merge(volcano_stats_up, TF_merged_up, by.x = "symbol", by.y = "hsapiens_homolog_associated_gene_name", all.x = T)
  
  volcano_stats_dn <- volcano_stats[volcano_stats$lfc < 0,]
  volcano_stats_dn <- merge(volcano_stats_dn, TF_merged_dn, by.x = "symbol", by.y = "hsapiens_homolog_associated_gene_name", all.x = T)
  volcano_stats_dn %<>%
    mutate(log_p_value = log_p_value*(-1),
           log_p_value_score = log_p_value_score*(-1))
  
  volcano_stats_tf <- rbind(volcano_stats_up, volcano_stats_dn)
  volcano_stats_tf %<>%
    filter(!is.na(log_p_value))
  
    
  return(volcano_stats_tf)
}

TF_merged_BL_VS_BM <- join_TF_human_genes(ME_dn = BL_VS_BM_decrease, ME_up = BL_VS_BM_increase)
TF_merged_MEM_VS_BM <- join_TF_human_genes(ME_dn = MEM_VS_BM_decrease, ME_up = MEM_VS_BM_increase)
TF_merged_AP_VS_BM <- join_TF_human_genes(ME_dn = AP_VS_BM_decrease, ME_up = AP_VS_BM_increase)
TF_merged_MEM_VS_BL <- join_TF_human_genes(ME_dn = MEM_VS_BL_decrease, ME_up = MEM_VS_BL_increase)
TF_merged_AP_VS_BL <- join_TF_human_genes(ME_dn = AP_VS_BL_decrease, ME_up = AP_VS_BL_increase)

save(TF_merged_BL_VS_BM, 
     TF_merged_MEM_VS_BM,
     TF_merged_AP_VS_BM,
     TF_merged_MEM_VS_BL,
     TF_merged_AP_VS_BL,
     file = "data/processed/data_fig_5b.rda")


TF_merged_label <- TF_merged %>%
  arrange(desc(mean_tf)) %>%
  mutate(Score_rank = 1:nrow(TF_merged)) %>%
  arrange(desc(log_p_value)) %>%
  mutate(log_p_rank = 1:nrow(TF_merged)) %>%
  filter(Score_rank %in% 1:10 | log_p_rank %in% 1:10)
  
TF_merged$col <- ifelse(TF_merged$fdr <= 0.05, "red", NA)

plt <- TF_merged_AP_VS_BL %>%
  ggplot(aes(x = log_p_value_score, y = lfc, label = symbol)) +
  geom_point(pch = 21, fill = TF_merged$col, aes(size = tgt_pct)) +
  geom_text_repel(data = TF_merged_label, size = 2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted") +
  xlab(bquote("Directional log<sub>2</sub>(*P*), motif enrichment")) +
  ylab(bquote("Mean "*log[10]*"(ChEA3 scor"*e^-1*"), across species")) +
  stat_cor(method = "spearman") +
  theme_rgb() +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_text(),
        legend.position = "bottom")


