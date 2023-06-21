library(clusterProfiler)
library(DESeq2)
library(ggrepel)
library(tidyverse)

make_out_dir <- function(out_dir) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
}

get_group_names <- function() {
  
  group_names <- c("WD_Ctl",
                   "WD_LPS",
                   "CD_Ctl",
                   "CD_LPS")
  
  groups <- c(paste0("\U03B2", "-glucan diet (unstim.)"),
              paste0("\U03B2", "-glucan diet (stim.)"),
              "Control diet (unstim.)",
              "Control diet (stim.)")
  
  names(groups) <- group_names
  
  return(groups)
  
}

######################################
## differential expression analysis ##
######################################

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
      design = ~ group)
    
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
    
    write_tsv(tab, file = paste0(out_dir, "de.all_res.tsv"))
    
    sig <- tab %>%
      filter(padj <= 0.05 & abs(log2FoldChange) >= 1)
    
    write_tsv(sig, file = paste0(out_dir, "de.sig_res.tsv"))
    
  } else {
    
    print("Output detected. Skipping...")
    print(target)
    
  }
}

wrapper_deseq2 <- function(file, out_dir = "output/deseq2") {
  
  tmp <- file %>%
    gsub("\\.txt", "", .) %>%
    gsub(".*\\/", "", .)
  
  groups <- tmp %>%
    str_split("vs") %>%
    unlist()
  
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

import_deseq2_res <- function(path) {
  
  groups <- path %>%
    str_split("/") %>%
    unlist() %>%
    nth(length(.) - 1) %>%
    str_split("vs") %>%
    unlist()
  
  if (all(any(grepl("CD", groups)) & any(grepl("WD", groups)))) {
    
    value <- -1
    
    ref <- groups %>%
      purrr::map(\(x) if (grepl("WD", x)) x) %>%
      unlist()
    
  } else {
    
    value <- 1
    
    ref <- groups[1]
    
  }
  
  ref_group <- get_group_names()[[ref]]
  
  read_delim(path) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    mutate(log2FoldChange = value * log2FoldChange) %>%
    mutate(colour = case_when(
      padj <= 0.05 & log2FoldChange >= +1 ~ paste("\U02191 in", ref_group),
      padj <= 0.05 & log2FoldChange <= -1 ~ paste("\U02193 in", ref_group),
      padj > 0.05 | abs(log2FoldChange) < 1 ~ "NS"
    )) %>%
    arrange(colour)
  
}

make_volcano_plot <- function(de) {
  
  groups <- colnames(de) %>%
    purrr::map(\(x) if (grepl("^CD|^WD", x)) x) %>%
    unlist() %>%
    str_replace("_._fpm", "") %>%
    unique()
  
  if (all(any(grepl("CD", groups)) & any(grepl("WD", groups)))) {
    
    ref <- groups %>%
      purrr::map(\(x) if (grepl("WD", x)) x) %>%
      unlist()
    
  } else {
    
    ref <- groups[1]
    
  }
  
  ref_group <- get_group_names()[[ref]]
  
  tmp_ann <- de %>%
    filter(padj <= 0.05 & abs(log2FoldChange) >= 2) 
  
  if (nrow(tmp_ann) == 0) {
    tmp_ann <- de %>%
      filter(padj <= 0.05 & abs(log2FoldChange) >= 1) 
  }
  
  ann <- tmp_ann %>%
    select(ENSEMBL, log2FoldChange, pvalue, colour) %>%
    top_n(-20, pvalue) %>%
    inner_join(annotables::grcm38[,c(1, 3)], by = c("ENSEMBL" = "ensgene"))
  
  pal <- c("salmon", "steelblue", "darkgrey")
  
  group_names <- c(paste("\U02191 in", ref_group),
                   paste("\U02193 in", ref_group),
                   "NS")
  
  names(pal) <- group_names
  
  y_intercept <- de %>%
    filter(colour != "NS") %>%
    filter(pvalue == max(pvalue)) %>%
    pull(pvalue) %>%
    log10() * - 1
  
  ggplot(de, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_hline(yintercept = y_intercept, colour = "lightgrey", linetype = "dashed") +
    geom_vline(xintercept = c(-2, -1, 1, 2), colour = "lightgrey", linetype = "dashed") +
    geom_point(aes(colour = colour, fill = colour), pch = 21, alpha = 0.75) +
    ggrepel::geom_text_repel(data = ann,
                             aes(x = log2FoldChange, y = -log10(pvalue), label = symbol),
                             min.segment.length = 1,
                             size = 3,
                             fontface = "italic") +
    theme_classic() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = ggtext::element_markdown()) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    labs(x = "**log<sub>2</sub> fold change**",
         y = "**-log<sub>10</sub>(*p*)**")
  
}

############################
## GO enrichment analysis ##
############################

run_go_enrichment_analysis <- function(de) {
  
  gene <- de %>%
    filter(padj <= 0.05) %>%
    pull(ENSEMBL)
  
  universe <- de %>%
    pull(ENSEMBL)
  
  enriched_go <- enrichGO(gene = gene,
                          OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                          keyType = "ENSEMBL",
                          ont = "BP",
                          pvalueCutoff = 0.01,
                          universe = universe,
                          qvalueCutoff = 0.05,
                          readable = TRUE)
  
  enriched_go@result %>%
    as_tibble() %>%
    separate(GeneRatio, c("x", "y"), sep = "/") %>%
    mutate(GeneRatio = as.numeric(x) / as.numeric(y)) %>%
    dplyr::select(-c(x, y))
  
}

make_enriched_go_plot <- function(enriched_go, top = 10) {
  
  p_inp <- enriched_go %>%
    top_n(-top, qvalue) %>%
    mutate(GO = str_trunc(paste0(ID, ": ", Description), 50))
  
  write_tsv(p_inp, "tmp/")
  
  ggplot(p_inp, aes(x = GeneRatio, y = reorder(GO, GeneRatio))) +
    geom_point(aes(fill = -log10(pvalue), size = Count), pch = 21, colour = "black") +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme_classic() +
    theme(legend.title = ggtext::element_markdown(),
          axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown()) +
    labs(x = "**Gene ratio**",
         y = "**GO term**",
         fill = "**-log<sub>10</sub>(*p*)**",
         size = "**Count**")
  
}

cluster_axis <- function(data) {
  
  if (nrow(data) > 1) {
    
    data %>%
      dist() %>%
      hclust() %>%
      as.dendrogram() %>%
      labels()
    
  } else {
    
    data %>%
      rownames()
    
  }
}

make_hmap <- function(go_term, de, enriched_go, top = 25) {
  
  tmp <- enriched_go %>%
    as_tibble() %>%
    filter(ID == go_term)
  
  title <- paste0(tmp[[1, 1]], ": ", tmp[[1, 2]])
  
  genes <- tmp %>%
    pull(geneID) %>%
    str_split("/") %>%
    unlist()
  
  p_inp <- de %>%
    inner_join(annotables::grcm38, by = c("ENSEMBL" = "ensgene")) %>%
    filter(symbol %in% genes)
  
  if (nrow(p_inp) > top) {
    
    p_inp <- p_inp %>%
      top_n(top, log2FoldChange)
    
  }
  
  p_inp <- p_inp %>%
    select(symbol, contains("_fpm")) %>%
    pivot_longer(!symbol, names_to = "sample_id", values_to = "fpm") %>%
    mutate(sample_id = gsub("_fpm", "", sample_id)) %>%
    group_by(symbol) %>%
    mutate(z_score = scale(log(fpm + 1))[,1]) %>%
    ungroup() %>%
    distinct() %>%
    mutate(diet = ifelse(grepl("WD", sample_id), "\U03B2-glucan diet", "Control diet")) %>%
    mutate(stim = ifelse(grepl("Ctl", sample_id), "(unstim.)", "(stim.)")) %>%
    mutate(group = paste0("**", diet, " ", stim, "**"))
  
  group_order <- p_inp %>%
    select(group) %>%
    distinct() %>%
    mutate(n1 = ifelse(!grepl("Control", group), 0, 1)) %>%
    mutate(n2 = ifelse(!grepl("unstim", group), 0, 1)) %>%
    mutate(n = n1 + n2) %>%
    arrange(n) %>%
    pull(group)
  
  p_inp$group <- factor(p_inp$group, levels = group_order)
  
  tmp <- p_inp %>%
    select(symbol, sample_id, z_score) %>%
    pivot_wider(names_from = sample_id, values_from = z_score) %>%
    column_to_rownames("symbol") %>%
    as.matrix()
  
  x_ord <- tmp %>%
    t() %>%
    cluster_axis()
  
  y_ord <- cluster_axis(tmp)
  
  p_inp$sample_id <- factor(p_inp$sample_id, levels = x_ord)
  
  p_inp$symbol <- factor(p_inp$symbol, levels = y_ord)
  
  ggplot(p_inp, aes(x = sample_id, y = symbol)) +
    facet_grid(~group, scales = "free_x") +
    geom_tile(aes(fill = z_score), colour = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    scale_x_discrete(expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(mult = c(0, 0)),
                     position = "right") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"),
          legend.title = ggtext::element_markdown(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(face = "italic"),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          panel.spacing = unit(1.5, "lines"),
          strip.text.x = ggtext::element_markdown()) +
    ggtitle(title) +
    labs(fill = "***z*-score**")
  
}

##########
## GSEA ##
##########

import_gsea_inputs <- function(de, gs_cat, gs_sub = "") {
  
  dat <- de %>%
    select(1, 3) %>%
    setNames(c("gene_id", "logFC"))
  
  gene_list <- dat$logFC
  names(gene_list) <- dat$gene_id
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gene_set <- msigdbr::msigdbr(species = "Mus musculus",
                               category = gs_cat,
                               subcategory = gs_sub) %>%
    dplyr::select(3:6)
  
  output <- list("gene_list" = gene_list, "gene_set" = gene_set)
  
  return(output)
  
}

run_gsea <- function(de,
                     gs_cat,
                     gs_sub = "",
                     out_dir,
                     rerun = FALSE) {
  
  gs <- paste0("MSigDB_", gs_cat, "_", gs_sub) %>%
    str_replace("_$", "") %>%
    str_replace("\\:", "_")
  
  out_dir <- gsub("/$", "", out_dir)
  dir_t <- paste0(out_dir, "/tables/")
  output <- paste0(dir_t, gs, ".res.tsv")
  
  if (!file.exists(output) | rerun == TRUE) {
    
    make_out_dir(dir_t)
    
    gsea_inputs <- import_gsea_inputs(de, gs_cat, gs_sub)
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

make_gsea_fig <- function(path, gs_cat, gs_sub, pwy, out_dir) {
  
  groups <- path %>%
    str_split("/") %>%
    unlist() %>%
    nth(length(.) - 1) %>%
    str_split("vs") %>%
    unlist()
  
  if (all(any(grepl("CD", groups)) & any(grepl("WD", groups)))) {
    
    groups <- rev(groups)
    
    value <- -1
    
  }
  
  labels <- get_group_names()[groups] %>%
    as.character()
  
  de <- import_deseq2_res(path) %>%
    mutate(log2FoldChange = value * log2FoldChange)
  
  run_gsea(de, gs_cat, gs_sub, out_dir)
  
  gsea_inputs <- import_gsea_inputs(de, gs_cat, gs_sub)
  
  gene_list   <- gsea_inputs$gene_list
  gene_set    <- gsea_inputs$gene_set
  
  gs <- paste0("MSigDB_", gs_cat, "_", gs_sub) %>%
    str_replace("_$", "") %>%
    str_replace("\\:", "_")
  
  result <- paste0(out_dir, "/tables/", gs, ".res.tsv") %>%
    read_delim() %>%
    filter(ID == pwy)
  
  filtered_gs <- gene_set %>%
    filter(gs_name == pwy)
  
  x <- clusterProfiler::GSEA(gene_list,
                             TERM2GENE = filtered_gs[, c(1, 4)],
                             pvalueCutoff = 1)
  
  y <- arrange(x, desc(NES))
  
  tmp_p <- enrichplot::gseaplot2(y,
                                 geneSetID = 1,
                                 title = pwy,
                                 subplots = 1:2,
                                 color = "lightgreen")
  
  title <- str_split(pwy, "_") %>%
    unlist() %>%
    tail(-1) %>%
    str_c(collapse = " ") %>%
    str_to_sentence() %>%
    str_wrap()
  
  qval <- result %>% pull(p.adjust)
  
  if (qval < 0.001) {
    qval <- scales::scientific(qval)
  } else {
    qval <- round(qval, 3)
  }
  
  nes  <- result %>% pull(NES) %>% round(3)
  
  subtitle <- paste0("*q*-value = ", qval, " and NES = ", nes)
  
  p_a <- tmp_p[[1]] +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          plot.subtitle = ggtext::element_markdown()) +
    ggtitle(title) +
    labs(subtitle = subtitle)
  
  x_min <- min(tmp_p[[2]]$data$x)
  x_max <- max(tmp_p[[2]]$data$x)
  
  ann <- tibble(x = c(x_min, x_max),
                y = 0,
                label = labels,
                hjust = c(0, 1))
  
  p_b <- tmp_p[[2]] +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(b = 1, unit = "cm"),
          plot.background = element_rect(fill = "transparent")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    geom_text(data = ann,
              aes(x = x, y = y, label = label, hjust = hjust),
              vjust = 1.25,
              colour = c("salmon", "steelblue")) +
    coord_cartesian(clip = "off")
  
  p_a  / p_b +
    plot_layout(heights = c(1, 0.125), tag_level = "new")
  
}
