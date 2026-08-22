# ==============================================================================
# directionality_test.R
#
# Decomposes SW-SA into the part that measures CLUSTERING (same-phase cells
# placed near each other) and the part that measures ORDERING (cells actually
# progressing G1 -> S -> G2M rather than the reverse).
#
#   SW-SA = P(phase increases) + P(tie)
#   An ordering that clusters but has no directional preference scores (1+T)/2.
#   DIRECTIONAL = SW-SA - (1+T)/2 = (P(increase) - P(decrease)) / 2
#
# THREE EXPERIMENTS
#   A. Window sweep      -- is near-zero directionality an artefact of w = 10?
#   B. Exact sign test   -- is the directional component distinguishable from 0?
#   C. Order-blindness   -- do Rand Index / NMI / accuracy change when an
#                           ordering is REVERSED? (They should not. That is the
#                           point.)
#
#   Rscript analysis/directionality_test.R
# ==============================================================================

PROJ <- getwd()   # run from the repository root
CACHE_DIR  <- file.path(PROJ, "cache")
RESULT_DIR <- file.path(PROJ, "results")
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)
suppressPackageStartupMessages({library(Seurat); library(igraph)})

WINDOWS <- c(1, 2, 5, 10, 20, 50, 100, 200)

# ------------------------------------------------------------------ engine
find_adaptive_origin <- function(x, y)
  optim(c(median(x), median(y)),
        function(p) var(sqrt((x-p[1])^2 + (y-p[2])^2)))$par
unwrap_theta <- function(th) {
  d <- diff(th); d <- d - 2*pi*round(d/(2*pi)); c(th[1], th[1] + cumsum(d))
}
acpl_arc <- function(coords, span = 0.25) {
  uv <- find_adaptive_origin(coords[,1], coords[,2])
  th <- atan2(coords[,2]-uv[2], coords[,1]-uv[1]); ord <- order(th)
  tu <- unwrap_theta(th[ord])
  arc <- c(0, cumsum(sqrt(diff(cos(tu))^2 + diff(sin(tu))^2)))
  predict(loess(arc ~ tu, span = span))[order(seq_along(coords[,1])[ord])]
}
mst_pseudotime <- function(coords, root_phase, phases) {
  g  <- igraph::graph_from_adjacency_matrix(as.matrix(dist(coords)),
          weighted = TRUE, mode = "undirected")
  tr <- igraph::mst(g); ri <- which(phases == root_phase)[1]
  as.numeric(igraph::distances(tr, v = ri, weights = igraph::E(tr)$weight))
}
slingshot_pseudotime <- function(coords, phases, start_phase = "G1") {
  keep <- as.character(phases) %in% c("G1","S","G2M")
  sds <- slingshot::slingshot(data = coords[keep,,drop=FALSE],
           clusterLabels = as.character(phases)[keep],
           start.clus = start_phase, reducedDim = "UMAP")
  pt <- rep(NA_real_, nrow(coords)); pt[keep] <- rowMeans(
    slingshot::slingPseudotime(sds), na.rm = TRUE); pt
}

# phase-rank vector sorted by pseudotime, NAs dropped
ranked <- function(tau, phases) {
  pm <- c(G1=1L, S=2L, G2M=3L)
  d <- data.frame(t = tau, p = as.character(phases))
  d <- d[!is.na(d$t) & d$p %in% names(pm), , drop = FALSE]
  d <- d[order(d$t), , drop = FALSE]
  list(pn = pm[d$p], t = d$t)
}

# ============================== A + B ========================================
decompose <- function(pn, w) {
  n <- length(pn); if (n <= w) return(NULL)
  a <- pn[1:(n-w)]; b <- pn[(w+1):n]
  inc <- mean(a < b); dec <- mean(a > b); tie <- mean(a == b)
  swsa <- inc + tie
  # exact sign test on NON-OVERLAPPING windows only (independent pairs)
  idx <- seq(1, n-w, by = w)
  a2 <- pn[idx]; b2 <- pn[idx + w]
  nd <- sum(a2 != b2); k <- sum(a2 < b2)
  bt <- if (nd > 0) binom.test(k, nd, 0.5) else NULL
  list(w = w, swsa = 100*swsa, tie = 100*tie,
       noDir = 100*(1+tie)/2, direc = 100*(swsa - (1+tie)/2),
       n_indep = nd, k = k,
       p = if (is.null(bt)) NA_real_ else bt$p.value,
       ci_lo = if (is.null(bt)) NA_real_ else 100*bt$conf.int[1],
       ci_hi = if (is.null(bt)) NA_real_ else 100*bt$conf.int[2])
}

# ================================ C ==========================================
# order-blind clustering metrics, computed from a contingency table
cont_metrics <- function(true_lab, pred_lab) {
  tab <- table(true_lab, pred_lab); n <- sum(tab)
  ch2 <- function(x) x*(x-1)/2
  a  <- sum(ch2(tab)); ai <- sum(ch2(rowSums(tab))); bj <- sum(ch2(colSums(tab)))
  tot <- ch2(n)
  RI  <- (tot + 2*a - ai - bj) / tot
  exp_a <- ai*bj/tot
  ARI <- (a - exp_a) / ((ai+bj)/2 - exp_a)
  p  <- tab/n; pr <- rowSums(p); pc <- colSums(p)
  mi <- sum(ifelse(p > 0, p*log(p/outer(pr, pc)), 0))
  hr <- -sum(pr[pr>0]*log(pr[pr>0])); hc <- -sum(pc[pc>0]*log(pc[pc>0]))
  NMI <- if (hr > 0 && hc > 0) mi/sqrt(hr*hc) else NA_real_
  # accuracy under the best of the 6 label matchings
  perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  acc <- max(vapply(perms, function(pm)
    mean(true_lab == pm[pred_lab]), numeric(1)))
  c(RI = RI, ARI = ARI, NMI = NMI, Acc = acc)
}

sets <- c(nestorowa="Nestorowa HSC", buettner="Buettner mESC",
          leng="Leng hESC", richard="Richard T cells")

sweep_rows <- list(); blind_rows <- list()

for (tag in names(sets)) {
  f <- file.path(CACHE_DIR, paste0("cache_", tag, ".rds"))
  if (!file.exists(f)) { cat("\n## missing", f, "\n"); next }
  obj <- readRDS(f); coords <- obj$coords; phases <- obj$phases
  cat("\n\n################", sets[[tag]], " (N =", length(phases), ")################\n")

  pts <- list(ACPL = acpl_arc(coords))
  pts$MST       <- tryCatch(mst_pseudotime(coords,"G1",phases), error=function(e) NULL)
  pts$Slingshot <- tryCatch(slingshot_pseudotime(coords,phases,"G1"), error=function(e) NULL)
  pts <- pts[!vapply(pts, is.null, logical(1))]

  # ---------------- A + B ----------------
  cat("\n[A+B] WINDOW SWEEP and EXACT SIGN TEST on independent windows\n")
  cat("      DIRECTIONAL = SW-SA minus the score a non-directional ordering\n")
  cat("      with the same tie rate would get. Sign test H0: P(inc)=P(dec).\n\n")
  for (m in names(pts)) {
    r <- ranked(pts[[m]], phases)
    cat(sprintf("  -- %s --\n", m))
    cat("     w   SW-SA   ties   no-dir   DIRECTIONAL   n_indep   P(inc|indep)      95% CI        p\n")
    for (w in WINDOWS) {
      d <- decompose(r$pn, w); if (is.null(d)) next
      cat(sprintf("  %4d %7.2f %6.2f %8.2f %13.2f %9d %10.3f  [%.3f,%.3f] %8.4f\n",
          d$w, d$swsa, d$tie, d$noDir, d$direc, d$n_indep,
          if (d$n_indep>0) d$k/d$n_indep else NA, d$ci_lo/100, d$ci_hi/100, d$p))
      sweep_rows[[length(sweep_rows)+1]] <- data.frame(
        Dataset=sets[[tag]], Method=m, w=d$w, SWSA=round(d$swsa,3),
        Ties=round(d$tie,3), NoDir=round(d$noDir,3), Directional=round(d$direc,3),
        n_indep=d$n_indep, P_inc=round(if(d$n_indep>0) d$k/d$n_indep else NA,4),
        CI_lo=round(d$ci_lo/100,4), CI_hi=round(d$ci_hi/100,4), p=signif(d$p,4))
    }
    cat("\n")
  }

  # ---------------- C ----------------
  cat("[C] ORDER-BLINDNESS: same ordering, then REVERSED\n")
  cat("    Pseudotime cut into 3 equal-size bins -> predicted phase labels.\n")
  cat("    If a metric is order-blind its value will not move.\n\n")
  cat("  method       direction     RI      ARI      NMI      Acc    SW-SA(w=10)\n")
  for (m in names(pts)) {
    r <- ranked(pts[[m]], phases)
    for (dir in c("forward","REVERSED")) {
      tt <- if (dir == "forward") r$t else -r$t
      o  <- order(tt); pn_o <- r$pn[o]
      pred <- cut(seq_along(pn_o), breaks = 3, labels = FALSE)
      cm <- cont_metrics(as.integer(pn_o), as.integer(pred))
      dd <- decompose(pn_o, 10)
      cat(sprintf("  %-12s %-10s %7.4f %8.4f %8.4f %8.4f %10.2f\n",
                  m, dir, cm["RI"], cm["ARI"], cm["NMI"], cm["Acc"], dd$swsa))
      blind_rows[[length(blind_rows)+1]] <- data.frame(
        Dataset=sets[[tag]], Method=m, Direction=dir,
        RI=round(cm["RI"],5), ARI=round(cm["ARI"],5),
        NMI=round(cm["NMI"],5), Acc=round(cm["Acc"],5),
        SWSA_w10=round(dd$swsa,3), row.names=NULL)
    }
  }
}

sw <- do.call(rbind, sweep_rows); bl <- do.call(rbind, blind_rows)
if (is.null(sw) || !nrow(sw)) stop("No rows produced: the dataset cache is missing or every method failed. Refusing to write an empty table.")
write.csv(sw, file.path(RESULT_DIR, "directionality_sweep.csv"), row.names=FALSE)
if (is.null(bl) || !nrow(bl)) stop("No rows produced: the dataset cache is missing or every method failed. Refusing to write an empty table.")
write.csv(bl, file.path(RESULT_DIR, "order_blindness.csv"),      row.names=FALSE)

cat("\n\n================= SUMMARY: directional component at each w =================\n")
if (!is.null(sw)) {
  agg <- tapply(sw$Directional, list(sw$Dataset, sw$w), function(x) round(mean(x),2))
  print(agg)
  cat("\n(mean across the three methods, percentage points)\n")
}
cat("\nWritten to ", RESULT_DIR, "\n", sep="")
