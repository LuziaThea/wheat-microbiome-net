# Extension of the lambda path beyond the StARS-selected boundary (Supplementary Table 19).
# For the Pseudomonas-Z. tritici network, StARS edge instability stayed below 0.05 along the whole path
# (nlambda = 50, lambda.min.ratio = 1e-2), so the smallest lambda was selected. This script refits the
# selected rank on the full data along a longer path with the same log spacing, down to lambda.min.ratio = 1e-3,
# and summarises the network at selected path positions (no StARS; one fit, ~15 min for 400 taxa).
#
# Usage: Rscript 09_lambda_path_extension.R <selected model .rds from 03> <Meta_info_... .txt from 03> <rank> <out.csv>
#   The Meta_info table lists the taxa in model order (column Name; P_ = Pseudomonas, Z_ = Z. tritici).
suppressMessages({library(SpiecEasi); library(Matrix)})
a <- commandArgs(TRUE)
model <- readRDS(a[1]); r <- as.integer(a[3])
meta <- read.table(a[2], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
nm <- meta$Name[order(meta$ID)]
X <- model$est$data                                   # clr-transformed data exactly as used by spiec.easi()
stopifnot(ncol(X) == length(nm))
mx <- pulsar::getMaxCov(cov(X))
lam50 <- pulsar::getLamPath(mx, mx * 1e-2, 50, log = TRUE)          # path used for all networks
step <- log(lam50[1]) - log(lam50[2])
lam <- exp(log(lam50[1]) - step * (0:(50 + ceiling(log(10) / step) - 1)))   # same spacing, down to 1e-3
fit <- SpiecEasi:::sparseLowRankiCov(X, r = r, lambda = lam)
res <- data.frame()
for (k in unique(c(40, 45, 50, 55, 60, 67, 75, length(lam)))) {
  if (k > length(lam)) next
  A <- as.matrix(fit$path[[k]]) != 0; pc <- -cov2cor(as.matrix(fit$icov[[k]]))
  ut <- which(upper.tri(A) & A, arr.ind = TRUE); w <- pc[ut]
  x <- nm[ut[, 1]]; y <- nm[ut[, 2]]; cross <- substr(x, 1, 1) != substr(y, 1, 1)
  s70 <- w[(x == "P_S70" | y == "P_S70") & cross]
  res <- rbind(res, data.frame(index = k, lambda = signif(lam[k], 4), edges = nrow(ut),
                               negative_pct = round(100 * mean(w < 0), 1), cross_pct = round(100 * mean(cross), 1),
                               S70_Zt_pos = sum(s70 > 0), S70_Zt_neg = sum(s70 < 0)))
}
print(res)
write.csv(res, a[4], row.names = FALSE)
