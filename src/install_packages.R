.install_packages <- function(i, packages, repos) {
  
  package <- packages[i]
  
  if (!package %in% installed.packages()[, "Package"]) {
    
    if (repos[i] == "Bioconductor") {
      
      if(!require(BiocManager)) install.packages("BiocManager")
      
      BiocManager::install(package, update = FALSE, ask = FALSE)
      
    } else if (repos[i] == "CRAN") {
      
      install.packages(package, repos = "http://cran.us.r-project.org")
      
    } else {
      
      if(!require(devtools)) install.packages("devtools")
      
      devtools::install_github("stephenturner/annotables")
      
    }
  } else {
    print(paste(package, "already installed"))
  }
}

install_packages <- function() {
  
  packages <- c("annotables",
                "clusterProfiler",
                "cowplot",
                "DESeq2",
                "enrichplot",
                "ggrepel",
                "ggtext",
                "org.Mm.eg.db",
                "patchwork",
                "msigdbr",
                "scales",
                "tidyverse")
  
  if (!all(packages %in% installed.packages()[, "Package"])) {
  
    repos <- c("GitHub",
               "Bioconductor",
               "CRAN",
               "Bioconductor",
               "Bioconductor",
               "CRAN",
               "CRAN",
               "Bioconductor",
               "CRAN",
               "CRAN",
               "CRAN",
               "CRAN")
    
    n <- length(packages)
    
    lapply(1:n, .install_packages, packages, repos)
    
  } else {
    
    print("All packages already installed")
    
  }
}
