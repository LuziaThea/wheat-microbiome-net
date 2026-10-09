# Pathogen suppressor assessment (Figure 5A, B; Supplementary Tables 8 and 16).
# Classification and plots follow LS041_Biocontrol_assessment_PS_ITS_Domi.Rmd (D. Stalder); edge weights are the
# partial correlations of the StARS-selected edges of the Pseudomonas-fungal network.
# Usage: Rscript 05_pathogen_suppressor_assessment.R <SLR_model.rds> <Meta_info_...txt from 03> <05_Meta_info_..._trophy.txt (pathogen annotation, in the repository)> <out_prefix>
#   main network: SLR_model_..._rank45.rds; Supplementary Table 16: SLR_model_..._noS44S45_rank40.rds (no plots needed)
#   The trophy table holds the manual phytopathogen annotation of the 30 fungal genera (column phytopathogen_manual).
suppressMessages({library(Matrix); library(SpiecEasi); library(igraph)})
args <- commandArgs(TRUE); OUT <- args[4]

# load info on pathogenicity
tro <- read.table(file = args[3], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
meta <- data.frame(name = tro$Name[1:30], is_pathogen = tro$phytopathogen_manual[1:30] == "pathogen")

# adjacency matrix: partial correlations of the selected edges, labelled with the network meta table
labs <- read.table(args[2], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")$Name
model <- readRDS(args[1]); icov <- as.matrix(model$est$icov[[model$select$stars$opt.index]])
adj <- -cov2cor(icov) * (as.matrix(getRefit(model)) != 0); adj <- (adj + t(adj)) / 2; diag(adj) <- 0
dimnames(adj) <- list(labs, labs)
adj <- adj[!(labs %in% c(NA, "NA", "unassigned_fungal_genus")), !(labs %in% c(NA, "NA", "unassigned_fungal_genus"))]
sel <- !(row.names(adj) %in% meta$name[is.na(meta$is_pathogen)])
adj <- adj[sel, sel]
# type of interactions
type <- data.frame(name=row.names(adj), type=character(nrow(adj)), n=rep(NA, nrow(adj)))

# pathogens
type$type[row.names(adj) %in% meta$name[meta$is_pathogen]] <- "pathogen"

# order 1
type$type[(type$type=="") & (rowSums(adj[, type$type=="pathogen"] < 0)>0) & (rowSums(adj[, type$type=="pathogen"] > 0)>0)] <- "mixed_o1"
type$type[(type$type=="") & (rowSums(adj[, type$type=="pathogen"] < 0)>0)] <- "supressor_o1"
type$type[(type$type=="") & (rowSums(adj[, type$type=="pathogen"] > 0)>0)] <- "facilitator_o1"

# order 2
type$type[(type$type=="") & ((colSums(adj[type$type=="mixed_o1",] > 0) > 0) | ((colSums(adj[type$type=="supressor_o1",] > 0) > 0) & (colSums(adj[type$type=="facilitator_o1",] > 0) > 0)))] <- "mixed_o2"
type$type[(type$type=="") & (colSums(adj[type$type=="facilitator_o1",] > 0) > 0)] <- "facilitator_o2"
type$type[(type$type=="") & (colSums(adj[type$type=="supressor_o1",] > 0) > 0)] <- "supressor_o2"

# order 3
type$type[(type$type=="") & ((colSums(adj[type$type=="mixed_o2",] > 0) > 0) | ((colSums(adj[type$type=="supressor_o2",] > 0) > 0) & (colSums(adj[type$type=="facilitator_o2",] > 0) > 0)))] <- "mixed_o3"
type$type[(type$type=="") & (colSums(adj[type$type=="facilitator_o2",] > 0) > 0)] <- "facilitator_o3"
type$type[(type$type=="") & (colSums(adj[type$type=="supressor_o2",] > 0) > 0)] <- "supressor_o3"

# those with not even 3 order interaction w. pathogens
type$type[type$type==""] <- "aprox_neutral"
type$type <- as.factor(type$type) 

# conveart to factor
type$type <- as.factor(type$type)

is_pathogen <- row.names(adj) %in% meta$name[meta$is_pathogen]
is_supressing <- !(is_pathogen) & (colSums(adj[is_pathogen,] < 0) > 0)

# FIX: fixed factor levels, so colours/positions stay correct even if a category is empty
type$type <- factor(as.character(type$type), levels=c("aprox_neutral","facilitator_o1","facilitator_o2","facilitator_o3","mixed_o1","mixed_o2","mixed_o3","pathogen","supressor_o1","supressor_o2","supressor_o3"))
write.csv(type, paste0(OUT,"_node_types.csv"), row.names=FALSE)
print(table(type$type))
# suppressors with their target pathogens (partial correlations), and suppressor stabilizers
pat <- type$name[type$type == "pathogen"]
for (s in type$name[type$type == "supressor_o1"]) { v <- adj[s, pat]; v <- v[v != 0]; cat(s, ":", paste(sprintf("%s %.4f", names(v), v), collapse = "; "), "\n") }
cat("stabilizers:", paste(type$name[type$type == "supressor_o2"], collapse = ", "), "\n")
library(igraph)

net <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=T)
V(net)$type <- type$type

# setting for the edges
E(net)$sign <- sign(E(net)$weight)
E(net)$weight <- abs(E(net)$weight)

E(net)$width <- 0.2 + 3 * abs(E(net)$weight) / max(abs(E(net)$weight))  # CHANGE: width proportional to |partial correlation|, thinner scale
E(net)$alpha <- 0.15 + 0.85*(((E(net)$weight - min(E(net)$weight)) / (max(E(net)$weight) - min(E(net)$weight)))^0.5)  # change tranparency by interactions strength
E(net)$color <- ifelse(E(net)$sign > 0, "#335d8c", "#cf3829")  # CHANGE: colours as in the published Figure 5 (blue positive, red negative)
E(net)$color <- mapply(function(color, alpha) {
  grDevices::adjustcolor(color, alpha.f = alpha)
}, color=E(net)$color, alpha=E(net)$alpha)  # adjust color transparency


# settings for the vertex
V(net)$name <- sub("^g__|^P_", "", type$name)
V(net)$name[V(net)$name == "Exobasidiomycetes_gen_Incertae_sedis"] <- "Exobasidiomycetes"  # FIX: by name instead of position

V(net)$organism <- factor(ifelse(grepl("^P_", type$name), "Pseudomonas", "Fungus"), levels = c("Fungus", "Pseudomonas"))  # FIX: from names instead of fixed counts
V(net)$type <- type$type
V(net)$importance <- c(5,2,3,4,2,3,4,1,2,3,4)[V(net)$type]

V(net)$shape <- c("square", "circle")[V(net)$organism]
V(net)$size <- 5
# CHANGE: node colours as in the published Figure 5; order = factor levels of type$type
V(net)$color <- c(aprox_neutral="#bfbdbd", facilitator_o1="#d7e046", facilitator_o2="#c6d27c", facilitator_o3="#c0c79d", mixed_o1="#5587b8", mixed_o2="#768dc4", mixed_o3="#9aa5c0", pathogen="#595a5b", supressor_o1="#6f59a6", supressor_o2="#9782bc", supressor_o3="#a99ec7")[as.character(V(net)$type)]


# create layout
set.seed(4)
lo <- layout_nicely(net)
lo[,1] <- ((V(net)$importance-1)/4) + runif(nrow(lo), -1, 1) * 0.085
lo[,2] <- runif(nrow(lo), -1, 1) 
lo[,2] <- lo[,2] + 2*c(0,-1,-1,-1,0,0,0,0,1,1,1)[type$type]
lo[V(net)$importance==1,2] <- lo[V(net)$importance==1,2] + ifelse(rowSums(adj[is_pathogen, is_supressing]<0)>0, 2, -2)
for(i in unique(V(net)$importance)){
  x <- rank(-lo[V(net)$importance==i,2])
  x <- max(x) - x + 1
  lo[V(net)$importance==i,2] <- 2*(((x-1)/(max(x)-1))-0.5)
}

# plot
pdf(file = paste0(OUT,"_Fig5A.pdf"), width = 4, height = 4)
par(mar=c(0,0,0,0))
plot(net, vertex.frame.color=rgb(0,0,0,0.15), vertex.label.color="black", vertex.label.cex=0.35, layout=lo, xlim=c(-1.1, 1.1), vertex.label.family="sans")
write.csv(data.frame(name=V(net)$name, type=as.character(V(net)$type), organism=as.character(V(net)$organism), x=lo[,1], y=lo[,2]), paste0(OUT,"_A_nodes.csv"), row.names=FALSE)
e <- igraph::as_data_frame(net, what="edges"); write.csv(data.frame(from=e$from, to=e$to, w=e$weight, sign=e$sign), paste0(OUT,"_A_edges.csv"), row.names=FALSE)
text(x=seq(-1,1,length.out=5)*0.85, y=1.12, label=c("Phytopathogens", "Direct Interactions", "Indirect Interactions\n(Order2)", "Indirect Interactions\n(Order 3)", "Marginal or No Effect"), font=2, cex=0.45, adj=0.5)
dev.off()
# create subgraph
net_biocontrol <- induced_subgraph(net, V(net)[V(net)$type %in% c("pathogen", "supressor_o1", "supressor_o2", "supressor_o3")])
sel <- components(net_biocontrol)$membership==1
# CHANGE: instead of removing Helgardiomyces and Bipolaris by name, drop pathogens without a negative edge to a
# first-order suppressor (with partial correlations Helgardiomyces is a target of Vishniacozyma and must stay)
sup1 <- V(net_biocontrol)$type == "supressor_o1"
Aneg <- as.matrix(as_adjacency_matrix(net_biocontrol, attr = "sign")) < 0
no_target <- V(net_biocontrol)$type == "pathogen" & rowSums(Aneg[, sup1, drop = FALSE]) == 0
sel[no_target] <- F
net_biocontrol <- induced_subgraph(net_biocontrol, V(net_biocontrol)[sel])

# get components of facilitating reactions within biocontrol strains
components_biocontrol <- induced_subgraph(net, V(net)[V(net)$type %in% c("supressor_o1", "supressor_o2", "supressor_o3")])
components_biocontrol <- components(net_biocontrol)$membership

# assign components
V(net_biocontrol)$component <- NA
V(net_biocontrol)$component[match(names(components_biocontrol), V(net_biocontrol)$name)] <- components_biocontrol
V(net_biocontrol)$component[V(net_biocontrol)$type=="pathogen"] <- 6

# color according to components
V(net_biocontrol)$color <- ifelse(V(net_biocontrol)$type=="pathogen", "#595a5b", c(supressor_o1="#6f59a6", supressor_o2="#9782bc", supressor_o3="#a99ec7")[as.character(V(net_biocontrol)$type)])  # CHANGE: colours as in published Figure 5B

# layout
set.seed(47)
lo <- layout_with_mds(net_biocontrol)
# lo[,1] <- jitter(lo[,1], factor=20)
# lo[,2] <- jitter(lo[,2], factor=20)
# CHANGE: the manual nudges of the original script (lo[4,2], lo[15,1], lo[25,1], lo[26,]) refer to node positions
# of the covariance-based subgraph and are not applied here; adjust label overlaps by hand if needed
if (isTRUE(as.logical(Sys.getenv("FIG5_ORIGINAL_NUDGES", "FALSE")))) {
lo[4,2] <- lo[4,2] - 0.5
lo[15,1] <- lo[15,1] + 0.5
lo[25,1] <- lo[25,1] - 0.35
lo[26,] <- lo[26,] - 0.5 }

# plot
pdf(file = paste0(OUT,"_Fig5B.pdf"), width = 3, height = 3)
par(mar=c(0,0,0,0))
plot(net_biocontrol, vertex.frame.color=rgb(0,0,0,0.15), vertex.label.color="black", vertex.label.cex=0.35, vertex.label.family="sans", layout=lo)
write.csv(data.frame(name=V(net_biocontrol)$name, type=as.character(V(net_biocontrol)$type), organism=as.character(V(net_biocontrol)$organism), x=lo[,1], y=lo[,2]), paste0(OUT,"_B_nodes.csv"), row.names=FALSE)
e <- igraph::as_data_frame(net_biocontrol, what="edges"); write.csv(data.frame(from=e$from, to=e$to, w=e$weight, sign=e$sign), paste0(OUT,"_B_edges.csv"), row.names=FALSE)
dev.off()
