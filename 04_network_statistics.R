# Network statistics of the selected networks (Figures 3A, 4A, 6A; Supplementary Table 15) and sensitivity of
# the key results to the SLR rank (Supplementary Table 18, network-level columns).
# Run 03_*.Rmd first. Usage (from the folder with the 03 output):
#   Rscript 04_network_statistics.R <Edge_weights_...txt> <Meta_info_...txt> <label>
#   Rscript 04_network_statistics.R --ranks <SLR_model_..._rank> <Meta_info_...txt> <label>   (all saved ranks)
suppressMessages({library(igraph); library(SpiecEasi)})
source("00_network_functions.R")
a <- commandArgs(TRUE)

summarise_edges <- function(e, meta, hubs = NULL) {
  org <- setNames(meta$Organism, meta$ID)
  e$type <- ifelse(org[as.character(e$Source)] != org[as.character(e$Target)], "cross-kingdom",
                   paste("within", org[as.character(e$Source)]))
  out <- data.frame(edges = nrow(e), negative = sum(e$partial_correlation < 0),
                    negative_pct = round(100 * mean(e$partial_correlation < 0), 1),
                    cross_kingdom_pct = round(100 * mean(e$type == "cross-kingdom"), 1))
  print(out); print(table(e$sign, e$type))
  deg <- table(c(e$Source_name, e$Target_name)); negdeg <- table(c(e$Source_name, e$Target_name)[rep(e$sign == "negative", 2)])
  top <- head(sort(deg, decreasing = TRUE), 5)
  neg_top <- as.integer(negdeg[names(top)]); neg_top[is.na(neg_top)] <- 0
  print(data.frame(taxon = names(top), degree = as.integer(top), negative = neg_top))
  out
}

modularity_test <- function(e, label, nperm = 1000) {
  # fast greedy modularity of the unweighted positive, negative and absolute networks vs. degree-preserving
  # rewired networks; p = (k + 1) / (n + 1), k = random networks with modularity >= observed
  res <- data.frame()
  for (k in c("positive", "negative", "absolute")) {
    ee <- switch(k, positive = e[e$partial_correlation > 0, ], negative = e[e$partial_correlation < 0, ], absolute = e)
    g <- simplify(graph_from_data_frame(ee[, c("Source_name", "Target_name")], directed = FALSE)); set.seed(1)
    q <- modularity(cluster_fast_greedy(g))
    null <- replicate(nperm, modularity(cluster_fast_greedy(rewire(g, keeping_degseq(niter = 10 * ecount(g))))))
    res <- rbind(res, data.frame(network = label, type = k, nodes = vcount(g), edges = ecount(g), Q = round(q, 3),
                                 null_mean = round(mean(null), 3), p = round((sum(null >= q) + 1) / (nperm + 1), 3)))
  }
  res
}

if (a[1] == "--ranks") {
  meta <- read.table(a[3], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = ""); label <- a[4]
  files <- Sys.glob(paste0(a[2], "*.rds")); files <- files[grepl("_rank[0-9]+\\.rds$", files)]
  out <- data.frame()
  for (f in files) {
    m <- readRDS(f); r <- as.integer(sub(".*_rank([0-9]+)\\.rds$", "\\1", f))
    if (is_degenerate(m)) { out <- rbind(out, data.frame(rank = r, ebic = NA, degenerate = TRUE, edges = NA, negative = NA, negative_pct = NA, cross_kingdom_pct = NA)); next }
    e <- selected_edges(m, meta$Name)
    write.table(e, paste0(label, "_rank", r, "_edges.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    s <- summarise_edges(e, meta)
    out <- rbind(out, data.frame(rank = r, ebic = round(ebic_logdet(m), 1), degenerate = FALSE, s))
    rm(m); gc()
  }
  out <- out[order(out$rank), ]; print(out)
  write.csv(out, paste0(label, "_rank_sensitivity_network.csv"), row.names = FALSE)
} else {
  e <- read.table(a[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
  meta <- read.table(a[2], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = ""); label <- a[3]
  s <- summarise_edges(e, meta)
  mod <- modularity_test(e, label); print(mod)
  write.csv(mod, paste0(label, "_modularity.csv"), row.names = FALSE)
}
