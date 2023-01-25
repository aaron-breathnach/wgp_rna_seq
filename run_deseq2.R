library(tidyverse)
library(DESeq2)

##################
## run_deseq2() ##
##################

run_deseq2 <- function(meta, id, group, count_data, out_dir, rerun = FALSE) {
  
  out_dir <- paste0(gsub("/$", "", out_dir), "/")
  
  target <- paste0(out_dir, "de.all_res.tsv")
  
  run_analysis <- any(!file.exists(target) | rerun == TRUE)
  
  if (run_analysis) {
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    col_data <- meta %>%
      select(all_of(c(id, group))) %>%
      setNames(c("sample_id", "group")) %>%
      mutate_if(is.character, as.factor) %>%
      column_to_rownames("sample_id")
    
    cts <- count_data %>%
      mutate_if(is.numeric, as.integer) %>%
      column_to_rownames("ENSEMBL_ID")
    
    dds <- DESeqDataSetFromMatrix(
      countData = cts,
      colData = col_data,
      design = ~ group
    )
    
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    dds <- DESeq(dds)
    res <- results(dds)
    
    temp <- res@elementMetadata %>%
      as.data.frame() %>%
      .[[2, 2]]
    
    coef <- gsub(" ", "_", gsub(".* group", "group", temp))
    
    resLFC <- lfcShrink(dds,
                        coef = coef,
                        type = "apeglm")
    
    norm_data <- fpm(dds) %>%
      as.data.frame() %>%
      rownames_to_column("ENSEMBL") %>%
      setNames(paste0(names(.), "_fpm")) %>%
      dplyr::rename(ENSEMBL = 1)
    
    tab <- resLFC %>%
      as.data.frame() %>%
      rownames_to_column("ENSEMBL") %>%
      replace(is.na(.), 1) %>%
      inner_join(norm_data, by = "ENSEMBL")
    
    write.table(tab,
                file = paste0(out_dir, "de.all_res.tsv"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
    
    sig <- tab %>%
      filter(padj <= 0.05 & abs(log2FoldChange) >= 1)
    
    write.table(sig,
                file = paste0(out_dir, "de.sig_res.tsv"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
  } else {
    
    print("Output detected. Skipping...")
    print(target)
    
  }
}

wrapper_deseq2 <- function(file, out_dir = "output/deseq2") {
  
  tmp <- file %>%
    gsub("\\.txt", "", .) %>%
    gsub(".*\\/", "", .)
  
  if (grepl("^GSE", tmp)) {
    sub_dir <- gsub("_.*", "", tmp)
  } else {
    sub_dir <- tmp
  }
  
  out_dir <- paste0(gsub("/$", "", out_dir), "/", sub_dir, "/")
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  count_data <- read_tsv(file) %>%
    setNames(gsub("\\-", "_", names(.)))
  
  sample_ids <- colnames(count_data)[-1]
  groups <- gsub("..$", "", colnames(count_data)[-1])
  meta <- tibble(sample_id = sample_ids, group = groups)
  
  run_deseq2(meta,
             id = "sample_id",
             group = "group",
             count_data,
             out_dir,
             rerun = FALSE)
  
}
