# Helper functions used by the 03 scripts (source this file first).
# Tested with R 4.3.3, SpiecEasi 1.1.3, pulsar 0.3.11, igraph 2.0.3.

library(Matrix)

## Count matrices exactly as they enter spiec.easi(): one row per leaf, samples shared by both objects.
## merge = TRUE sums PCR replicates and sequencing runs per leaf (merge_samples(group = "Sample")).
counts_from_physeq <- function(physeq_first, physeq_second, merge_first = TRUE, merge_second = TRUE) {
  if (merge_first)  physeq_first  <- phyloseq::merge_samples(physeq_first,  group = "Sample", fun = sum)
  if (merge_second) physeq_second <- phyloseq::merge_samples(physeq_second, group = "Sample", fun = sum)
  m <- function(p) { M <- as(phyloseq::otu_table(p), "matrix"); if (phyloseq::taxa_are_rows(p)) M <- t(M); M }
  A <- m(physeq_first); B <- m(physeq_second)
  common <- intersect(rownames(A), rownames(B))
  list(A = A[common, , drop = FALSE], B = B[common, , drop = FALSE])
}

## Sparse precision matrix of the selected model and degenerate-fit check
## (a zero diagonal entry means a taxon received zero estimated precision; such fits are excluded).
sparse_icov <- function(model) as.matrix(model$est$icov[[model$select$stars$opt.index]])
is_degenerate <- function(model) any(diag(sparse_icov(model)) <= 0)

## EBIC as in SpiecEasi::ebic(), with the log-likelihood computed from the log-determinant.
## log(det()) overflows for large networks (all values were -Inf for the 400-taxon network).
## Note: like SpiecEasi, only the sparse component enters the likelihood and the degrees of freedom;
## the low-rank component is not penalised.
ebic_logdet <- function(model, gamma = 0.5) {
  o <- model$select$stars$opt.index; X <- model$est$data
  icov <- as.matrix(model$est$icov[[o]]); S <- cov(X)
  z <- Matrix::rowSums(model$est$path[[o]]) != 0; q <- sum(!z)
  Sz <- icov[z, z, drop = FALSE]
  ll <- as.numeric(determinant(Sz, logarithm = TRUE)$modulus) - sum(diag(Sz %*% S[z, z])) - (ncol(X) - q)
  n <- nrow(X); df <- sum(model$refit$stars) / 2
  -n * ll + log(n) * df + 4 * gamma * log(ncol(X)) * df
}

## StARS-selected edges with partial correlations (edge weight) and the covariance used previously.
## rho_ij = -theta_ij / sqrt(theta_ii * theta_jj), Theta = sparse component at the selected lambda.
selected_edges <- function(model, names) {
  icov <- sparse_icov(model); sel <- as.matrix(SpiecEasi::getRefit(model)) != 0
  stopifnot(length(names) == nrow(icov))
  pc <- -cov2cor(icov)
  cv <- tryCatch(solve(icov), error = function(e) solve(icov + diag(1e-6, nrow(icov))))
  ut <- which(upper.tri(sel) & sel, arr.ind = TRUE)
  data.frame(Source = ut[, 1], Target = ut[, 2], Source_name = names[ut[, 1]], Target_name = names[ut[, 2]],
             partial_correlation = pc[ut], covariance = cv[ut],
             sign = ifelse(pc[ut] > 0, "positive", "negative"), stringsAsFactors = FALSE)
}

## Robustness: networks estimated at the lambda selected for the full data, for the full data (id "full") and
## for random subsamples of 80% of the samples (taxa without reads in a subsample are dropped).
## The full lambda path is computed because the SLR solver uses warm starts along the path.
fixed_lambda_networks <- function(counts, r, opt, ids = c("full", 1:5), nlambda = 50, lambda.min.ratio = 1e-2) {
  A <- counts[[1]]; B <- counts[[2]]
  Xf <- SpiecEasi:::.spiec.easi.norm(list(A, B))
  mx <- pulsar::getMaxCov(cov(Xf))
  lam <- pulsar::getLamPath(mx, mx * lambda.min.ratio, nlambda, log = TRUE)
  out <- list()
  for (id in ids) {
    if (id == "full") { X <- Xf; cols <- seq_len(ncol(Xf)) } else {
      set.seed(as.integer(id)); keep <- sample(rownames(A), size = round(nrow(A) * 0.8))
      ka <- which(colSums(A[keep, ]) > 0); kb <- which(colSums(B[keep, ]) > 0); cols <- c(ka, ncol(A) + kb)
      X <- SpiecEasi:::.spiec.easi.norm(list(A[keep, ka], B[keep, kb]))
    }
    e <- SpiecEasi:::sparseLowRankiCov(X, r = r, lambda = lam)
    out[[as.character(id)]] <- list(adj = as.matrix(e$path[[opt]]) != 0, cols = cols)
  }
  out
}

## Relative Hamming distance (FP+FN)/P (Kurtz et al. 2019). all_taxa = TRUE (default, as in the manuscript) counts
## the reference edges of taxa absent from a subsample as not recovered; FALSE uses only taxa present in both.
hamming_distances <- function(nets, all_taxa = TRUE) {
  ref <- nets[["full"]]$adj
  sapply(setdiff(names(nets), "full"), function(id) {
    s <- nets[[id]]
    if (all_taxa) { S <- matrix(FALSE, nrow(ref), ncol(ref)); S[s$cols, s$cols] <- s$adj; R <- ref } else { R <- ref[s$cols, s$cols]; S <- s$adj }
    ut <- upper.tri(S)
    (sum(S[ut] & !R[ut]) + sum(R[ut] & !S[ut])) / sum(R[ut])
  })
}
