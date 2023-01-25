library(tidyverse)

make_out_dir <- function(out_dir) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
}

gsea_bar <- function(sig_res, out_dir, gs, title = NULL, rerun = FALSE) {
  
  out_dir <- paste0(gsub("/$", "", out_dir), "/")
  output  <- paste0(out_dir, gs, ".bar.png")
  
  if (!file.exists(output) | rerun == TRUE) {
  
    make_out_dir(out_dir)
    
    legend_title <- "-log<sub>10</sub>(*p*)"
    
    bar_inp <- sig_res %>%
      mutate(log10_pval = -log10(pvalue))
    
    p <- ggplot(bar_inp, aes(x = NES, y = reorder(ID, NES))) +
      geom_bar(aes(fill = log10_pval),
               stat = "identity",
               colour = "black",
               width = 1) +
      theme_bw(base_size = 12.5) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.title = ggtext::element_markdown(),
        axis.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold")
      ) +
      scale_y_discrete(expand = expansion(mult = c(0, 0)),
                       labels = function(x) { str_trunc(x, 30) }) +
      ggtitle(title) +
      labs(x = "Normalised effect size",
           y = "Gene set",
           fill = legend_title) +
      scale_fill_distiller(palette = "Reds", direction = 1)
    
    min_h <- 2
    
    max_h <- 12.5
    
    h <- 0.5 * nlevels(as.factor(bar_inp$ID))
    
    if (h < min_h) {
      h <- min_h
    } else if (h > max_h) {
      h <- max_h
    }
    
    ggsave(p,
           file = output,
           height = h,
           width = 8.75,
           dpi = 300)
    
  } else {
    print(paste0(output, " already exists. Skipping..."))
  }
}

gs_hmap <- function(pwy, gene_set, df, out_dir, rerun = FALSE) {
  
  out_dir <- paste0(gsub("/$", "", out_dir), "/")
  output  <- paste0(out_dir, pwy, ".png")
  
  if (!file.exists(output) | rerun == TRUE) {
    
    make_out_dir(out_dir)
    
    filtered_gs <- gene_set %>%
      filter(gs_name == pwy)
    
    ensembl_to_symbol <- filtered_gs %>%
      select(ensembl_gene, gene_symbol)
    
    gene_ids <- filtered_gs$ensembl_gene
    
    hmap_in <- df %>%
      filter(gene_id %in% gene_ids) %>%
      inner_join(ensembl_to_symbol, by = c("gene_id" = "ensembl_gene")) %>%
      select(ncol(.), 2:(ncol(.) - 1)) %>%
      pivot_longer(!gene_symbol, names_to = "sample_id", values_to = "fpkm") %>%
      group_by(gene_symbol, sample_id) %>%
      summarise(fpkm = sum(fpkm)) %>%
      ungroup() %>%
      pivot_wider(names_from = "sample_id", values_from = "fpkm") %>%
      column_to_rownames("gene_symbol") %>%
      as.matrix()
    
    n <- nrow(hmap_in)
    
    show_row_names <- ifelse(n < 50, TRUE, FALSE)
    
    png(filename = output, width = 5, height = 8, units = "in", res = 300)
    
    p <- ComplexHeatmap::pheatmap(log(hmap_in + 0.001),
                                  fontsize = 7.5,
                                  fontsize_col = 10,
                                  fontsize_row = 10,
                                  main = str_trunc(pwy, 50),
                                  name = "z-score",
                                  scale = "row",
                                  show_rownames = show_row_names,
                                  border_color = NA)
    
    print(p)
    
    dev.off()
    
  } else {
    print(paste0(output, " already exists. Skipping..."))
  }
}

enrichment_plot <- function(pwy, gene_list, gene_set, out_dir, rerun = FALSE) {
  
  out_dir <- paste0(gsub("/$", "", out_dir), "/")
  output <- paste0(out_dir, pwy, ".png")
  
  if (!file.exists(output) | rerun == TRUE) {
    
    make_out_dir(out_dir)
    
    filtered_gs <- gene_set %>%
      filter(gs_name == pwy)
    
    x <- clusterProfiler::GSEA(gene_list, TERM2GENE = filtered_gs)
    y <- arrange(x, desc(NES))
    
    p <- enrichplot::gseaplot2(y, geneSetID = 1, title = pwy)
    
    ggsave(p, filename = output, height = 6, width = 6, dpi = 300)
    
  } else {
    print(paste0(output, " already exists. Skipping..."))
  }
}

import_gsea_inputs <- function(DE, GS_CAT, GS_SUB = "", METHOD = "DESeq2") {
  
  if (!(METHOD %in% c("DESeq2", "Novogene"))) {
    stop("METHOD must be one of \"DESeq2\" or \"Novogene\"")
  }
  
  if (METHOD == "DESeq2") {
    cols <- c(1, 3)
  } else {
    cols <- c(1, 8)
  }
  
  dat <- read_delim(DE, col_select = all_of(cols)) %>%
    setNames(c("gene_id", "logFC"))
  
  gene_list <- dat$logFC
  names(gene_list) <- dat$gene_id
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gene_set <- msigdbr::msigdbr(species = "Mus musculus",
                               category = GS_CAT,
                               subcategory = GS_SUB) %>%
    dplyr::select(3:6)
  
  output <- list("gene_list" = gene_list, "gene_set" = gene_set)
  
  return(output)
  
}

run_gsea <- function(DE,
                     GS_CAT,
                     GS_SUB = "",
                     OUT_DIR,
                     METHOD = c("DESeq2", "Novogene"),
                     rerun = FALSE
                     ) {
  
  GS <- paste0("MSigDB_", GS_CAT, "_", GS_SUB) %>%
    str_replace("_$", "") %>%
    str_replace("\\:", "_")
  
  out_dir <- gsub("/$", "", OUT_DIR)
  dir_t <- paste0(out_dir, "/tables/")
  output <- paste0(dir_t, GS, ".res.tsv")
  
  if (!file.exists(output) | rerun == TRUE) {
  
  make_out_dir(dir_t)
  
  gsea_inputs <- import_gsea_inputs(DE, GS_CAT, GS_SUB, METHOD)
  gene_list   <- gsea_inputs$gene_list
  gene_set    <- gsea_inputs$gene_set
  
  res <- clusterProfiler::GSEA(gene_list,
                               TERM2GENE = gene_set[, c(1, 4)],
                               pvalueCutoff = 1)
  
  result <- res@result %>%
    arrange(p.adjust)
  
  write.table(res,
              output,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  } else {
    print(paste0(output, " already exists. Skipping..."))
  }
}

make_gsea_plot <- function(result, de, out_dir, method = "DESeq2") {
  
  tmp <- result %>%
    str_replace(".*/", "") %>%
    str_replace("MSigDB_", "") %>%
    str_replace("\\..*", "")
  
  gs_cat <- tmp %>%
    str_replace("_.*", "")
  
  gs_sub <- ifelse(grepl("_", tmp), gsub(".*_", "", tmp), "")
  
  GS <- paste0("MSigDB_", gs_cat, "_", gs_sub) %>%
    str_replace("_$", "") %>%
    str_replace("\\:", "_")
  
  gsea_inputs <- import_gsea_inputs(de, gs_cat, gs_sub, method)
  gene_list   <- gsea_inputs$gene_list
  gene_set    <- gsea_inputs$gene_set
  
  sig_res <- read_delim(result) %>%
    as_tibble() %>%
    filter(p.adjust <= 0.05)
  
  pathways <- sig_res %>%
    top_n(-50, p.adjust) %>%
    pull(ID)
  
  n_sig <- length(pathways)
  
  if (n_sig > 0) {
    
    if (is.null(out_dir)) {
      out_dir <- paste0(dirname(dirname(result)), "/")
    } else {
      out_dir <- paste0(gsub("/$", "", out_dir), "/")
    }
    
    make_out_dir(out_dir)
    
    ######################
    ## enrichment plots ##
    ######################
    
    dir_e <- paste0(out_dir, "enrichment_plots/", GS, "/")
    
    lapply(pathways, enrichment_plot, gene_list, gene_set[, c(1, 4)], dir_e)
    
    ##############
    ## heatmaps ##
    ##############
    
    df <- read_delim(de, col_select = c(1, contains("_fp"))) %>%
      setNames(gsub("_fp.*", "", names(.))) %>%
      dplyr::rename(gene_id = 1)
    
    dir_h <- paste0(out_dir, "heatmaps/", GS, "/")
    
    lapply(pathways, gs_hmap, gene_set, df, dir_h)
    
    ##############
    ## bar plot ##
    ##############
    
    dir_b <- paste0(out_dir, "bar_plots/")
    
    p <- gsea_bar(sig_res, dir_b, GS)
    
  }
}

wrapper_gsea <- function(
    de,
    out_dir,
    gs_cat = "H",
    gs_sub = "",
    method = "DESeq2",
    rerun = FALSE
) {
  
  out_dir <- paste0(gsub("/$", "", out_dir), "/")
  make_out_dir(out_dir)
  sub_dir <- paste0(out_dir, gsub(".*\\/", "", dirname(de)), "/")
  
  run_gsea(
    DE = de,
    GS_CAT = gs_cat,
    GS_SUB = gs_sub,
    OUT_DIR = sub_dir,
    METHOD = method
  )
  
  res <- paste0("MSigDB_", gs_cat, "_", gs_sub) %>%
    str_replace("_$", "") %>%
    str_replace("\\:", "_") %>%
    str_c(sub_dir, "tables/", ., ".res.tsv")
  
  make_gsea_plot(
    result = res,
    de = de,
    out_dir = sub_dir
  )
  
}
