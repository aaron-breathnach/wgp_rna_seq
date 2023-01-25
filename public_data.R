setwd("~/Desktop/help/fred_rna_seq/")

source("src/run_deseq2.R")
source("src/volcano.R")
source("src/run_gsea.R")

ftp <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195750/suppl/GSE195750_GeneMatrix_RawCounts.txt.gz"

cnt <- paste0("input/read_counts/", gsub(".*\\/", "", gsub(".gz", "", ftp)))

if (!file.exists(cnt)) {
  cmd <- paste0("wget ", ftp, " -O ", cnt, ".gz; gunzip ", cnt)
  system(cmd)
}

count_data <- read_delim(cnt)

meta <- tibble(sample_id = colnames(count_data)[-1]) %>%
  mutate(group = ifelse(grepl("PBS", sample_id), "PBS", "WGP")) %>%
  mutate(group = as.factor(group))

run_deseq2(meta = meta,
           id = "sample_id",
           group = "group",
           count_data = count_data,
           out_dir = "output/public_data/deseq2",
           rerun = TRUE)

volcano <- volcano_plot(inp = "output/public_data/deseq2/de.all_res.tsv",
            gene_id = "ensgene")

ggsave(volcano,
       file = "output/public_data/deseq2/volcano.png",
       width = 5,
       height = 5,
       dpi = 300)

##########
## GSEA ##
##########

run_gsea(
  DE = "output/public_data/deseq2/de.all_res.tsv",
  GS_CAT = "H",
  GS_SUB = "",
  OUT_DIR = "output/public_data/gsea",
  METHOD = "DESeq2"
)

run_gsea(
  DE = "output/public_data/deseq2/de.all_res.tsv",
  GS_CAT = "C2",
  GS_SUB = "CP:KEGG",
  OUT_DIR = "output/public_data/gsea",
  METHOD = "DESeq2"
)
