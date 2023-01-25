setwd("~/Desktop/help/fred_rna_seq/")

source("src/run_gsea.R")
source("src/volcano.R")

files <- list.files("DEGlist",
                    pattern = "DEG.xls",
                    recursive = TRUE,
                    full.names = TRUE)

out_dir <- "output/"

volcanoes <- lapply(files, wrapper_volcano, "output/volcano_plots/pngs/")

pdf("output/volcano_plots/volcano_plots.pdf", width = 5, height = 5)
print(volcanoes)
dev.off()

## H
lapply(files[3], run_gsea, "H", "", out_dir)
## C2 KEGG
lapply(files[3], run_gsea, "C2", "CP:KEGG", out_dir)
## C2 REACTOME
lapply(files[3], run_gsea, "C2", "CP:REACTOME", out_dir)
## C7 IMMUNESIGDB
lapply(files[3], run_gsea, "C7", "IMMUNESIGDB", out_dir)


DE = files[3]
GS_CAT = "H"
GS_SUB = ""
OUT_DIR = "output/"
