setwd("~/Desktop/help/fred_rna_seq/")

source("src/run_deseq2.R")
source("src/volcano.R")
source("src/run_gsea.R")

##############################
## download public datasets ##
##############################

download_geo <- function(geo, base_dir = "input/read_counts/") {
  
  GEOquery::getGEOSuppFiles(geo,
                            makeDirectory = FALSE,
                            baseDir = base_dir)
  
}

accession_numbers <- c("GSE187464", "GSE195750")

lapply(accession_numbers, download_geo)

############
## DESeq2 ##
############

files <- lapply(tmp_files, function(x) {if (!grepl("xls", x)) { return(x) }})

lapply(files, wrapper_deseq2)

de_files <- list.files("output/deseq2",
                       pattern = "de.all_res.tsv",
                       recursive = TRUE,
                       full.names = TRUE)

##################
## volcano plot ##
##################

lapply(de_files, wrapper_volcano, out_dir)

##############
## run GSEA ##
##############

lapply(de_files[c(2, 6, 7)],
       wrapper_gsea,
       "output/gsea",
       gs_cat = "H",
       gs_sub = "",
       method = "DESeq2",
       rerun = FALSE)
