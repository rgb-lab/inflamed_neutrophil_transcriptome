#####################
# UTILITY FUNCTIONS #
#####################
#
# This file contains utility functions that are called by the analysis scripts.
#
#####################


###########
# IMPORTS #
###########
if (!require(httr)) install.packages("httr")
if (!require(jsonlite)) install.packages("jsonlite")
if (!require(tools)) install.packages("tools")

if (!require(BiocManager)) install.packages("BiocManager")
if (!require(SummarizedExperiment)) BiocManager::install("SummarizedExperiment")
if (!require(rhdf5)) BiocManager::install("rhdf5")
if (!require(GEOquery)) BiocManager::install("GEOquery")


#################
# DATA DOWNLOAD #
#################

#' \code {download_arch_data}
#' 
#' Downloads the ARCHS4 data matrix to local storage.
#' 
#' @param species target species, indicating, whether the "human" or the "mouse"
#' matrix should be downloaded.
#' @param target_directory the target directory to download the data to.
#' @param target_filename the filename of the downloaded file.
#' @param checksum an MD5 checksum that is used to check the integrity of
#' existing/previously downloaded files.
#' @param category one of "matrix"(gene counts), "transcript"(transcript counts)
#' or "tpm"(transcript tpm)
#' 
#' @returns the filename/path of the downloaded or existing file.
download_arch_data <- function (
  species = "human",
  target_directory = "./data/raw_data",
  target_filename = "matrix_download.h5",
  checksum = checksum,
  category = "matrix",
  version = "v11") {
  
  if (!species %in% c("human", "mouse") |
      !category %in% c("matrix", "transcript", "tpm") |
      !version %in% c("v11")) {
    stop ("Pleas supply valid species, category and version args. Refer to function docstring.")
  }
  
  destination_file <- get_destination_file(species, target_directory, target_filename)
  download_url <- paste0("https://s3.amazonaws.com/mssm-seq-matrix/", species, "_", category, "_", version, ".h5")

  # check if gene expression file was already downloaded to target dir,
  # if not in current directory download file form repository
  if(!file.exists(destination_file)){
    message ("Downloading compressed gene expression matrix...")
    download.file(download_url, destination_file, quiet = FALSE)
    message ("Saved file to ", destination_file)
  }
  
  message("Computing checksum of existing or previously downloaded file...")
  
  if(md5sum(destination_file) == checksum){
    message(paste0("Checksum of local ", species, " file correct."))
  } else {
    message(paste0("Checksum of species ,",
                   species,
                   " incorrect due to old or corrupt file."))
    message("Renaming old file with current date suffix and downloading new file...")
    destination_file_renamed <- paste0(sub('.([^.]*)$', '', destination_file), "_", Sys.Date(), ".h5")
    file.rename(destination_file, destination_file_renamed)
    message ("Downloading new compressed gene expression matrix...")
    download.file(download_url, destination_file, quiet = FALSE)
    message ("Saved file to ", destination_file)
  }
  
  return(destination_file)
}


#' \code {extract_arch_h5_to_se}
#' 
#' Extracts downloaded ARCHS4 file to a SummarizedExperiment data structure.
#' 
#' @param species target species.
#' @param destination_file the target file (as obtained from download_arch_data).
#' @param extract_samples a string vector containing the names of the samples
#' to extract from the data matrix.
#' 
#' @returns the SummarizedExperiment that was built from the ARCHS4 data.
extract_arch_h5_to_se <- function (species,
                                   destination_file,
                                   extract_samples = c("GSM996198")) {
  
  if (!species %in% c("Hs", "Mm")) stop ("Please define a valid species.")
  
  if (!file.exists(destination_file)) {
    stop ("ERROR: Destination file does not exist!\n",
          "Please run download_arch_data first and pass its return value as destination_file arg.")
  }
  
  # Retrieve count matrix from compressed data
  message("Retrieving count matrix.")
  sample_names <- h5read(destination_file, "meta/Sample_geo_accession")
  # Identify columns to be extracted
  sample_locations <- which(sample_names %in% extract_samples)
  
  genes <- h5read(destination_file, "meta/genes")
  
  # extract gene expression from compressed data
  counts <- as.data.frame(h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations)))
  
  rownames(counts) <- as.character(genes)
  colnames(counts) <- as.character(sample_names[sample_locations])
  
  # retrieve coldata from data
  message("Retrieving coldata.")
  sample_series <- h5read(destination_file, "meta/Sample_series_id")[sample_locations]
  sample_title <- h5read(destination_file, "meta/Sample_title")[sample_locations]
  sample_source_name <- h5read(destination_file, "meta/Sample_source_name_ch1")[sample_locations]
  sample_characteristics <-h5read(destination_file, "meta/Sample_characteristics_ch1")[sample_locations]
  sample_description <- h5read(destination_file, "meta/Sample_description")[sample_locations]
  H5close()
  coldata <- data.frame("sample" = sample_names[sample_locations], 
                        "series" = sample_series, 
                        "title" = sample_title,
                        "source_name" = sample_source_name,
                        "characteristics" = sample_characteristics,
                        "description" = sample_description,
                        "species" = species)
  
  rownames(coldata) <- as.character(coldata[["sample"]])
  coldata[["sample"]] <- NULL
  
  if(!identical(rownames(coldata), colnames(counts))) {
    stop ("Cannot create SummarizedExperiment from provided data, as coldata is not matching count matrix.")
  }
  
  message("Building SE object.")
  se <- SummarizedExperiment(assays = list(counts = counts),
                             colData = coldata)
  
  return(se)
  
}


#' \code {get_chea}
#' 
#' Takes a gene list and passes it to the ChEA3 API. Returns the ChEA3 results.
#' 
#' @param gene_set list of genes to analyze for enrichment.
#' @param query_name name of the query.
#' 
#' @returns a list of dataframes containing the ChEA3 results.
get_chea <- function (gene_set, query_name) {
  
  url = "https://maayanlab.cloud/chea3/api/enrich/"
  encode = "json"
  
  payload = list(query_name = query_name, gene_set = gene_set)
  
  #POST to ChEA3 server, get response as JSON
  response = POST(url = url, body = payload, encode = encode)
  json = httr::content(response, "text")
  
  #results as list of R dataframes
  results = fromJSON(json)
  
  return (results)
  
}


#' \code {get_enrichr}
#' 
#' Takes a gene list and passes it to the Enrichr API. Returns the Enrichr results.
#' 
#' @param gene_set list of genes to analyze for enrichment.
#' @param query_name name of the query.
#' 
#' @returns a list of dataframes containing the ChEA3 results.
add_enrichr <- function (gene_set, set_description) {
  
  url = "https://maayanlab.cloud/Enrichr/addList"
  encode = "multipart"
  
  gene_set_str <- paste(gene_set, collapse = "\n")
  
  payload = list(
    list = gene_set_str,
    description = set_description)
  
  # POST to Enrichr server, get response as JSON
  response <- POST(url = url, body = payload, encode = encode)
  json <- httr::content(response, "text")
  
  # parse results to list
  results <- fromJSON(json)
  
  # add meta
  results$submitted_genes <- gene_set
  results$submitted_description <- set_description
  
  return (results)
}

view_enrichr <- function (add_result) {
  
  userListId <- add_result$userListId
  
  url = paste0('https://maayanlab.cloud/Enrichr/view?userListId=', userListId)
  
  response <- GET(url = url)
  json <- httr::content(response, "text")
  
  results <- fromJSON(json)
  
  return(results)
}

validate_enrichr <- function (add_result) {
  view_enrichr_res <- view_enrichr(add_result)
  
  # TODO: handle removed mouse-genes (capitalized but removed, ...)
  
  if (all(view_enrichr_res$genes %in% add_result$submitted_genes) &&
      all(add_result$submitted_genes %in% view_enrichr_res$genes)) {
    add_result$passed_validation <- TRUE
    add_result$removed_genes <- c("")
    add_result$server_genes <- view_enrichr_res$genes
    return (add_result)
    
  } else if(all(view_enrichr_res$genes %in% add_result$submitted_genes) &&
            !all(add_result$submitted_genes %in% view_enrichr_res$genes)) {
    
    warning ("EnrichR filtered out these genes: ", paste(add_result$submitted_genes[
      !add_result$submitted_genes %in% view_enrichr_res$genes
    ], collapse = ", ") )
    
    add_result$passed_validation <- TRUE
    add_result$removed_genes <- add_result$submitted_genes[
      !add_result$submitted_genes %in% view_enrichr_res$genes
    ]
    add_result$server_genes <- view_enrichr_res$genes
    return (add_result)
  } 
  else if(all(view_enrichr_res$genes %in% str_to_upper(add_result$submitted_genes)) &&
            all(str_to_upper(add_result$submitted_genes) %in% view_enrichr_res$genes)) {
  
    warning ("EnrichR capitalized genes.")
    add_result$passed_validation <- TRUE
    add_result$removed_genes <- c("")
    add_result$server_genes <- view_enrichr_res$genes
    return (add_result)
    
  } else {
    
    add_result$passed_validation <- FALSE
    return (add_result)
  
  }
}

get_enrichr <- function (add_result, gene_set_library = 'KEGG_2015') {
  url = 'https://maayanlab.cloud/Enrichr/enrich'
  userListId <- add_result$userListId
  
  query_string <- paste0('?userListId=', userListId, '&backgroundType=', gene_set_library)
  
  query_url <- paste0(url, query_string)
  response <- GET(query_url, )
  
  json <- httr::content(response, "text")
  results <- fromJSON(json)
  
  return (results)
}

get_enrichr_df <- function (add_result, gene_set_library = 'KEGG_2015') {
  validated_query <- validate_enrichr(add_result)
  
  stopifnot(validated_query$passed_validation)
  
  get_enrichr_res <- get_enrichr(add_result, gene_set_library)
  
  res_df <- tibble(
    rank = as.numeric(lapply(get_enrichr_res[[1]], `[[`, 1)),
    term_name = as.character(lapply(get_enrichr_res[[1]], `[[`, 2)),
    p_value = as.numeric(lapply(get_enrichr_res[[1]], `[[`, 3)),
    z_score = as.numeric(lapply(get_enrichr_res[[1]], `[[`, 4)),
    combined_score = as.numeric(lapply(get_enrichr_res[[1]], `[[`, 5)),
    overlapping_genes = lapply(get_enrichr_res[[1]], `[[`, 6),
    adj_p_value = as.numeric(lapply(get_enrichr_res[[1]], `[[`, 7)),
    old_p_value = as.numeric(lapply(get_enrichr_res[[1]], `[[`, 8)),
    old_adj_p_value = as.numeric(lapply(get_enrichr_res[[1]], `[[`, 9)),
    submitted_genes = list(validated_query$submitted_genes),
    server_genes = list(validated_query$server_genes),
    query_description = validated_query$submitted_description,
    query_passed_validation = validated_query$passed_validation,
    query_removed_genes = list(validated_query$removed_genes),
    queried_library = gene_set_library
  )
  
  return (res_df)
}

# Function to retrieve goseq results.
get_goseq <- function(res, direction, input){
  
  fdr.threshold <- 0.05
  assayed.genes <- rownames(res)
  if(direction == "up"){
    de.genes <- rownames(res)[which(res$padj <= fdr.threshold & res$log2FoldChange >= 1)]
  }else{
    de.genes <- rownames(res)[which(res$padj <= fdr.threshold & res$log2FoldChange <= -1)]
  }
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  
  pwf <- nullp(gene.vector, "mm10", "geneSymbol")
  
  # cats <- rbind(data.frame(external_gene_name = up_genes, cat = "up"), data.frame(external_gene_name = dn_genes, cat = "dn"))
  
  cats <- rbind(data.frame(external_gene_name = input, cat = "up"))
  
  # cats %<>%
  #   left_join(ortho_human_to_mouse[,c("external_gene_name", "mmusculus_homolog_associated_gene_name")], by = "external_gene_name")
  
  up_mouse <- cats[cats$cat == "up", "external_gene_name"]
  # dn_mouse <-  cats[cats$cat == "dn", "mmusculus_homolog_associated_gene_name"]
  
  core_cats <- data.frame(gene = assayed.genes, cat = NA)
  # core_cats$cat <- ifelse(core_cats$gene %in% up_mouse, "core_up", ifelse(core_cats$gene %in% dn_mouse, "core_dn", NA))
  core_cats$cat <- ifelse(core_cats$gene %in% input, "core_up", NA)
  
  goseq.results <- goseq(pwf, genome = "mm10", gene2cat = core_cats, use_genes_without_cat = TRUE)
  goseq.results$hitsPerc <- goseq.results$numDEInCat*100/goseq.results$numInCat
  return(goseq.results)
}

#' \code {get_data_files}
#' 
#' Function to retrieve RNA-Seq data from the GEO and save it to local storage.
#' This function searches the GEO supplementary files for supplements containing
#' the requested (count- and meta-) data based on user input.
#' 
#' Currently, the data dir is hardcoded into the function.
#' 
#' @param study_name name of the study, used for local folder creation.
#' @param series_name a GEO series identifier (GSE...) to locate the series.
#' @param filter_GEO_filenames a regular expression used for filtering of the 
#' GEO supplementary files (e.g. "raw_counts")
#' @param compression a string indicating the file format of the supplementary
#' files to download
#' @param super_tar a bool indicating whether the file in question is a 
#' tar archive containing multiple e.g. gzipped files.
#' 
#' @returns a list containing the names of the downloaded files as well as
#' a relative path to the compressed and extracted files on local storage.
get_data_files <- function (study_name, series_name, filter_GEO_filenames, compression, super_tar = FALSE) {
  
  data_dir <- paste0("./data/raw_data/", study_name, "/", series_name)
  
  if (!dir.exists(data_dir)) {
    message ("Data dir for study ", study_name, "(", series_name, ")", " does not exist.\nCreating data dir.")
    dir.create(data_dir, recursive = TRUE)
  } else {
    message ("Data dir for study ", study_name, "(", series_name, ")", " already existing.\n")
  }
  
  filenames <- getGEOSuppFiles(series_name, makeDirectory = FALSE, fetch_files = FALSE,
                               baseDir = data_dir, filter_regex = filter_GEO_filenames)[["fname"]]
  
  filepaths <- paste (data_dir, filenames, sep = "/")
  
  
  if (any(!file.exists(filepaths))) {
    message ("Data files for study ", study_name, "(", series_name, ")", " not existing or corrupted.\nDownloading data files.\n")
    getGEOSuppFiles(series_name, makeDirectory = FALSE, fetch_files = TRUE,
                    baseDir = data_dir, filter_regex = filter_GEO_filenames)[["fname"]]
  } else {
    message ("Data files for study ", study_name, "(", series_name, ")", " already existing.\nSkipping data file download.\n")
  }
  
  
  if (super_tar) {
    
    message ("Untarring data archive files for study ", study_name, "(", series_name, ").\n")
    
    filenames <- NULL
    for (i in length (filepaths)) {
      filenames <- append (untar(filepaths[[i]], list = TRUE), filenames)
    }
    filepaths_extracted_tar <- paste (data_dir, filenames, sep = "/")
    
    if (any(!file.exists(filepaths_extracted_tar))) {
      
      for (i in length(filepaths)) {
        # untar and overwrite
        message ("Extracted data archive files for study ", study_name, "(", series_name, ")", " not existing or corrupted.\nExtracting data archive files.\n")
        untar(filepaths[[i]], exdir = data_dir)
      }
    } else {
      message ("Extracted data archive files for study ", study_name, "(", series_name, ")", " already existing.\nSkipping...\n")
    }
    
    filepaths <- filepaths_extracted_tar
    
  }
  
  
  if (compression == "gz") {
    filepaths_extracted <- sub(".gz$", "", filepaths)
    
    if (any(!file.exists(filepaths_extracted))) {
      message ("Extracted data files for study ", study_name, "(", series_name, ")", " not existing or corrupted.\nExtracting data files.\n")
      for (i in 1:length(filenames)) {
        gunzip(filename = filepaths[[i]], destname = filepaths_extracted[[i]], overwrite = TRUE, remove = FALSE)
      }
    } else {
      message ("Extracted data files for study ", study_name, "(", series_name, ")", " already existing.\nSkipping data file extraction\n")
    }
    
    
  } else if (compression %in% c("xlsx", "none")) {
    
    filepaths_extracted <- filepaths
    
    if (any(!file.exists(filepaths_extracted))) {
      stop ("Extracted data files for study ", study_name, "(", series_name, ")", " not existing or corrupted.\n.\n")
    }
    
  } else {
    # TODO: come up with a more elegant way to deal with this case
    stop ("Compression not implemented")
  }
  
  
  return (list(filenames, filepaths, filepaths_extracted))
  
}



##################
# DATA WRANGLING #
##################

#' \code {split_se}
#' 
#' Take a summarized experiment and an annotation category and get a list of 
#' SEobjects split by unique values of the annotation category.
#' 
#' @param SEobject the SEobject to split into smaller objects.
#' @param split_by a string specifying the category. This string must be one of
#' the SEobjects colDatas colnames
#' 
#' @returns a list of SEobjects split by unique values occuring in the specified
#' colData column.
split_se <- function(SEobject, split_by) {
  
  if (!is.character(split_by)) {
    stop ("Please pass a character string as split_by arg.")
  }
  
  if (length(split_by) != 1) {
    stop ("Currently this function only supports splitting objects by one category.")
  }
  
  split_levels <- factor(unique(SEobject[[split_by]]))
  
  split_objects <- list()
  
  for (i in 1:length(split_levels)) {
    
    split_objects[[i]] <- SEobject[ , SEobject[[split_by]] %in% split_levels[[i]]]
    
  }
  
  return(list(split_objects = split_objects, split_levels = split_levels))
  
}


#######################
# "PRIVATE" FUNCTIONS #
#######################

#' \code {get_destination_file}
#' 
#' Helper function used to generate a filepath to the desired archs4 download location.
#' 
#' @param species defines the current species.
#' @param target_directory defines the target directory relative to the
#' current workdir (intended to be R project root).
#' @param target_filename the filename for the downloaded file.
#' 
#' @returns the destination file.
get_destination_file <- function (species, target_directory, target_filename) {
  allowed_species <- c("Hs", "Mm")
  
  # check species
  if (!(length(species) == 1 & species %in% allowed_species)) {
    stop("Please define a valid species ('Hs', 'Mm').")
  }
  
  # check target dir
  if (!dir.exists(target_directory)) {
    stop("Target directory does not exist.")
  }
  
  if(!endsWith(target_filename, ".h5")) {
    stop("Please specify a .h5 download file.")
  }
  
  # checks passed
  # strip tailing slash
  target_directory <- ifelse(endsWith(target_directory, "/"), sub("/$", "", target_directory), target_directory)
  # build destfile
  destination_file <- paste(target_directory, paste(species, target_filename, sep = "_"), sep = "/")
  
  return(destination_file)
}
