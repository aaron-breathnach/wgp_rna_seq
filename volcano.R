library(tidyverse)
library(annotables)

volcano_plot <- function(inp,
                         gene_id = c("ensgene", "entrez", "symbol"),
                         q = 0.05,
                         title = NULL) {
  
  if (!(gene_id %in% c("ensgene", "entrez", "symbol"))) {
    error <- paste("gene_id must be one of \"ensgene\", \"entrez\", or \"symbol\"")
    stop(error)
  }
  
  cols <- c(1, 3, 5, 6)
  
  col_names <- c("gene_id", "log2FoldChange", "pvalue", "padj")
  
  tab <- read_delim(inp) %>%
    dplyr::select(all_of(cols)) %>%
    setNames(col_names) %>%
    mutate(colour = case_when(
      padj >  0.05 & abs(log2FoldChange) < 1 ~ "NS",
      padj >  0.05 & abs(log2FoldChange) >= 1 ~ "log<sub>2</sub> FC",
      padj <= 0.05 & abs(log2FoldChange) < 1 ~ "*p*-value",
      padj <= 0.05 & abs(log2FoldChange) >= 1 ~ "*p*-value and log<sub>2</sub> FC"
    )) %>%
    drop_na()
  
  x_min <- min(tab$log2FoldChange)
  x_max <- max(tab$log2FoldChange)
  x_val <- max(abs(c(x_min, x_max)))
  x_lim <- c(-x_val, x_val)
  
  ann <- tab %>%
    drop_na() %>%
    filter(padj <= q) %>%
    inner_join(annotables::grcm38, by = c("gene_id" = gene_id)) %>%
    select(1:7) %>%
    top_n(-10, padj)
  
  pal <- c("grey", "#ff5100", "#f00080", "#7d5df9")
                 
  names(pal) <- c("NS",
                  "log<sub>2</sub> FC",
                  "*p*-value",
                  "*p*-value and log<sub>2</sub> FC")
  
  tmp_max_sig <- tab %>%
    filter(padj <= q)
  
  tmp_volcano <- ggplot(tab, aes(x = log2FoldChange, y = -log10(pvalue)))
  
  if (nrow(tmp_max_sig) > q) {
    
    max_sig <- filter(tmp_max_sig, pvalue == max(pvalue)) %>%
      .[[1, "pvalue"]]
    
    tmp_volcano <- tmp_volcano +
      geom_hline(yintercept = -log10(max_sig), linetype = "dashed", colour = "darkgrey")
    
  }
  
  volcano <- tmp_volcano +
    geom_vline(xintercept = c(-2, -1, 1, 2), linetype = "dashed", colour = "darkgrey") +
    geom_point(aes(colour = colour, fill = colour), pch = 21, alpha = 0.75) +
    theme_bw() +
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
    ggtitle(title) +
    labs(x = "**log<sub>2</sub> fold change**",
         y = "**-log<sub>10</sub>(*p*)**") +
    xlim(x_lim) 
  
  if (nrow(ann) > 0) {
    
    volcano <- volcano +
      ggrepel::geom_text_repel(
        data = ann,
        aes(x = log2FoldChange, y = -log10(pvalue), label = symbol)
      )
    
  }
  
  return(volcano)
  
}

wrapper_volcano <- function(inp, q = 0.05, save = TRUE, out_dir = NULL) {
  
  tmp <- inp %>%
    str_replace("\\/de.all_res.tsv", "") %>%
    str_replace(".*\\/", "")
  
  title <- tmp %>%
    str_replace_all("vs", " vs ") %>%
    str_replace_all("_", "-")
  
  p <- volcano_plot(inp, gene_id = "ensgene", q, title)
  
  out_dir <- str_replace(inp, "de.all_res.tsv", "")
  
  if (save) {
    
    if(is.null(out_dir)) {
      stop("You must specify the output directory!")
    }
    
    filename <- tmp %>%
      str_c(".volcano.png")
    
    outdir <- out_dir %>%
      str_replace("/$", "") %>%
      str_c("/")
    
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }
    
    output <- str_c(outdir, filename)
    
    ggsave(p, file = output, width = 5, height = 5, dpi = 300)
    
  }
  
  return(p)
  
}
