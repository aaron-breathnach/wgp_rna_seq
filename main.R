source("src/install_packages.R")
source("src/source.R")

## installs packages if not installed already
install_packages()

## runs DESeq2
files <- list.files("input", full.names = TRUE)
lapply(files, wrapper_deseq2)

## import the DESeq2 results
de <- import_deseq2_res("output/deseq2/CD_CtlvsWD_Ctl/de.all_res.tsv")

## Fig S9D
fig_s9d <- make_volcano_plot(de) %>%
  cowplot::plot_grid(labels = "D.")

## Fig S9E
fig_s9e <- import_deseq2_res("output/deseq2/CD_LPSvsWD_LPS/de.all_res.tsv") %>%
  make_volcano_plot() %>%
  cowplot::plot_grid(labels = "E.")

## Fig 9E
enriched_go <- run_go_enrichment_analysis(de)

fig_9e <- make_enriched_go_plot(enriched_go) %>%
  cowplot::plot_grid(labels = "E.")

## Fig 9F
fig_9f <- make_gsea_fig("output/deseq2/CD_CtlvsWD_Ctl/de.all_res.tsv",
                        gs_cat = "H",
                        gs_sub = "",
                        pwy = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                        out_dir = "output/gsea/CD_CtlvsWD_Ctl") %>%
  cowplot::plot_grid(labels = "F.")

## Fig 9G
tmp_1 <- make_hmap("GO:0006119", de, enriched_go)

tmp_2 <- cowplot::plot_grid(plot.new(),
                          tmp_1,
                          rel_widths = c(0.05, 1),
                          nrow = 1)

fig_9g <- cowplot::plot_grid(plot.new(),
                             tmp_2,
                             labels = c("", "G."),
                             rel_heights = c(0.05, 1),
                             hjust = 0.1,
                             vjust = -0.1,
                             nrow = 2)

## add panel labels + save the plots
plot_list <- list(fig_s9d = fig_s9d,
                  fig_s9e = fig_s9e,
                  fig_9e = fig_9e,
                  fig_9f = fig_9f,
                  fig_9g = fig_9g)

plot_widths <- c(6.25, 6.25, 7.5, 5, 5)

n <- length(plot_list)

out_dir <- "output/figures/"

if (!dir.exists(out_dir)) dir.create(out_dir)

lapply(1:n, function(x) ggsave(paste0(out_dir, names(plot_list[x]), ".png"),
                               plot_list[[x]],
                               width = plot_widths[[x]],
                               height = 5,
                               dpi = 300,
                               bg = "white"))
