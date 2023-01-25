library(tidyverse)

setwd("~/Desktop/help/fred_rna_seq/")

read_counts <- read_delim("input/read_counts/readcount.xls")

x <- unique(gsub("_.$", "", colnames(read_counts)[-1]))

combos <- apply(combn(x, 2), 2, paste, collapse = "vs")

sub_cnt <- function(combo, read_counts, out_dir) {
  
  groups <- unlist(str_split(combo, "vs"))
  
  df <- read_counts %>%
    dplyr::select(1, contains(groups)) %>%
    dplyr::rename(ENSEMBL_ID = 1)
  
  outdir <- paste0(gsub("/$", "", out_dir), "/")
  
  filename <- paste0(outdir, combo, ".txt")
  
  write.table(df, filename, quote = FALSE, sep = "\t", row.names = FALSE)
  
}

lapply(combos, sub_cnt, read_counts, "input/read_counts")
