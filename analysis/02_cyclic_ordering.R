# ==============================================================================
# cyclic_ordering_test.R   -- SUPERSEDES global_ordering_test.R
#
# Two faults in the previous script are fixed here.
#
# FAULT 1 (mine). Jammalamadaka circular correlation centres on the circular
#   mean of the phase angles. Three phases at 0, 2pi/3, 4pi/3 have mean
#   resultant length 0.000 when balanced (Buettner) and 0.054 for Leng, so the
#   circular mean is undefined and the statistic is noise. Every rcirc value in
#   the previous run must be discarded.
#
# FAULT 2. Kendall tau-b treats G1 < S < G2M as a LINE. The same perfect cyclic
#   ordering scores +0.82 or -0.27 depending only on where it starts. MST and
#   Slingshot were given start = "G1"; ACPL is root-free and was given nothing.
#   The comparison was unfair to ACPL.
#
# WHAT THIS COMPUTES INSTEAD
#
#   1. Per-phase circular mean and resultant of the INFERRED angle. Rotation
#      invariant. Shows directly whether a method separates phases on the circle.
#
#   2. Cyclic order verdict. Forward gaps mu_S - mu_G1, mu_G2M - mu_S,
#      mu_G1 - mu_G2M each wrapped into [0, 2pi). They sum to 2pi if the phases
#      sit in the correct cyclic order and to 4pi if reversed. Origin-free.
#
#   3. Concentration C = mean over phases of the within-phase resultant length
#      of inferred angle. Zero if a phase is smeared around the circle, one if
#      it is a point. Permutation null.
#
#   4. tau_cyc = max Kendall tau over a grid of ORIGIN ROTATIONS and both
#      directions, with the identical max applied under permutation so the
#      selection is priced into the null. This is the origin-free version of
#      what the previous script computed unfairly.
#
#   Rscript analysis/cyclic_ordering_test.R
# ==============================================================================

PROJ <- getwd()   # run from the repository root
CACHE_DIR  <- file.path(PROJ, "cache")
RESULT_DIR <- file.path(PROJ, "results")
NROT  <- 36      # origin rotations
NPERM <- 1000
set.seed(2026)
suppressPackageStartupMessages({library(Seurat); library(igraph)})

find_adaptive_origin <- function(x,y)
  optim(c(median(x),median(y)), function(p) var(sqrt((x-p[1])^2+(y-p[2])^2)))$par
unwrap_theta <- function(th){d<-diff(th);d<-d-2*pi*round(d/(2*pi));c(th[1],th[1]+cumsum(d))}
acpl_arc <- function(cd, span=0.25){
  uv<-find_adaptive_origin(cd[,1],cd[,2]); th<-atan2(cd[,2]-uv[2],cd[,1]-uv[1])
  o<-order(th); tu<-unwrap_theta(th[o])
  arc<-c(0,cumsum(sqrt(diff(cos(tu))^2+diff(sin(tu))^2)))
  predict(loess(arc~tu,span=span))[order(seq_along(cd[,1])[o])]
}
mst_pt <- function(cd, phases, root="G1"){
  g<-igraph::graph_from_adjacency_matrix(as.matrix(dist(cd)),weighted=TRUE,mode="undirected")
  tr<-igraph::mst(g); ri<-which(phases==root)[1]
  as.numeric(igraph::distances(tr,v=ri,weights=igraph::E(tr)$weight))
}
sling_pt <- function(cd, phases, root="G1"){
  keep<-phases %in% c("G1","S","G2M")
  sds<-slingshot::slingshot(data=cd[keep,,drop=FALSE],clusterLabels=phases[keep],
        start.clus=root, reducedDim="UMAP")
  pt<-rep(NA_real_,nrow(cd)); pt[keep]<-rowMeans(slingshot::slingPseudotime(sds),na.rm=TRUE); pt
}

pm <- c(G1=1L,S=2L,G2M=3L)
circ_mean <- function(a) atan2(mean(sin(a)), mean(cos(a)))
circ_R    <- function(a) sqrt(mean(sin(a))^2 + mean(cos(a))^2)
wrap2pi   <- function(x) (x %% (2*pi))

# concordant minus discordant, O(n), positions all distinct so no x-ties
cd_stat <- function(pn){
  b1 <- c(0, head(cumsum(pn==1L),-1)); b2 <- c(0, head(cumsum(pn==2L),-1))
  b3 <- c(0, head(cumsum(pn==3L),-1))
  less <- ifelse(pn==1L, 0, ifelse(pn==2L, b1, b1+b2))
  more <- ifelse(pn==1L, b2+b3, ifelse(pn==2L, b3, 0))
  sum(less) - sum(more)
}
tau_denom <- function(pn){
  n<-length(pn); n0<-n*(n-1)/2
  n2<-sum(vapply(1:3,function(k){m<-sum(pn==k); m*(m-1)/2},numeric(1)))
  sqrt(n0*(n0-n2))
}
# max tau over origin rotations and both directions
tau_cyc <- function(pn, nrot = NROT){
  n <- length(pn); den <- tau_denom(pn)
  cuts <- unique(round(seq(0, n-1, length.out = nrot+1)[1:nrot]))
  best <- -Inf; arg <- NA_integer_; dir <- NA_character_
  for (r in cuts) {
    rolled <- if (r == 0) pn else c(pn[(r+1):n], pn[1:r])
    for (d in c(1,-1)) {
      s <- if (d == 1) rolled else rev(rolled)
      v <- d * 0 + cd_stat(s)/den
      if (v > best) { best <- v; arg <- r; dir <- if (d==1) "forward" else "reversed" }
    }
  }
  list(tau = best, rot = arg, dir = dir)
}

sets <- c(nestorowa="Nestorowa HSC", buettner="Buettner mESC",
          leng="Leng hESC", richard="Richard T cells")
rows <- list()

for (tag in names(sets)) {
  f <- file.path(CACHE_DIR, paste0("cache_",tag,".rds"))
  if (!file.exists(f)) { cat("\n## missing",f,"\n"); next }
  o <- readRDS(f); cd <- o$coords; ph <- as.character(o$phases)

  pts <- list(ACPL = acpl_arc(cd))
  pts$MST       <- tryCatch(mst_pt(cd, ph),   error=function(e) NULL)
  pts$Slingshot <- tryCatch(sling_pt(cd, ph), error=function(e) NULL)
  pts <- pts[!vapply(pts, is.null, logical(1))]

  cat("\n\n################", sets[[tag]], "################\n")

  for (m in names(pts)) {
    ok <- ph %in% names(pm) & !is.na(pts[[m]])
    t  <- pts[[m]][ok]; pn <- pm[ph[ok]]; n <- length(pn)
    ang <- 2*pi*(rank(t, ties.method="first") - 0.5)/n     # inferred angle

    cat(sprintf("\n  -- %s  (n = %d) --\n", m, n))

    # 1. per-phase circular mean and resultant
    mu <- vapply(1:3, function(k) circ_mean(ang[pn==k]), numeric(1))
    Rk <- vapply(1:3, function(k) circ_R(ang[pn==k]),    numeric(1))
    cat("     phase   mean angle   resultant R   (R near 0 = smeared round the circle)\n")
    for (k in 1:3)
      cat(sprintf("     %-6s %10.3f %13.3f\n", names(pm)[k], mu[k], Rk[k]))

    # 2. cyclic order verdict, origin-free
    g <- c(wrap2pi(mu[2]-mu[1]), wrap2pi(mu[3]-mu[2]), wrap2pi(mu[1]-mu[3]))
    # sum of forward gaps is 2pi if the phases sit in cyclic order G1,S,G2M and
    # 4pi if the same cycle is traversed the other way round. A "reversed"
    # verdict means correct ORDER with opposite HANDEDNESS, not a wrong order.
    verdict <- if (abs(sum(g) - 2*pi) < 1e-6) "correct" else
               if (abs(sum(g) - 4*pi) < 1e-6) "reversed-handedness" else "indeterminate"
    cat(sprintf("     cyclic order: %-20s gaps/(2pi/3) = %.2f %.2f %.2f\n",
                verdict, g[1]/(2*pi/3), g[2]/(2*pi/3), g[3]/(2*pi/3)))

    # 3. concentration with permutation null
    C  <- mean(Rk)
    nullC <- vapply(seq_len(NPERM), function(i){
      s <- pn[sample(n)]
      mean(vapply(1:3, function(k) circ_R(ang[s==k]), numeric(1)))
    }, numeric(1))
    pC <- (1 + sum(nullC >= C))/(NPERM+1)

    # 4. origin-free tau with matched null
    # BUG FIX: tau must see the phase sequence ORDERED BY PSEUDOTIME. The
    # previous version passed pn in original file order, so every method
    # returned an identical value that did not depend on the ordering at all.
    ord    <- order(t)
    pn_ord <- pn[ord]
    tc <- tau_cyc(pn_ord)
    nullT <- vapply(seq_len(NPERM), function(i) tau_cyc(pn_ord[sample(n)])$tau, numeric(1))
    pT <- (1 + sum(nullT >= tc$tau))/(NPERM+1)

    cat(sprintf("     concentration C = %.4f   null %.4f [%.4f,%.4f]   p = %.4f\n",
                C, mean(nullC), quantile(nullC,.025), quantile(nullC,.975), pC))
    cat(sprintf("     tau_cyc         = %.4f   null %.4f [%.4f,%.4f]   p = %.4f   (%s, rot %d)\n",
                tc$tau, mean(nullT), quantile(nullT,.025), quantile(nullT,.975),
                pT, tc$dir, tc$rot))

    rows[[length(rows)+1]] <- data.frame(
      Dataset=sets[[tag]], Method=m, n=n,
      R_G1=round(Rk[1],4), R_S=round(Rk[2],4), R_G2M=round(Rk[3],4),
      CyclicOrder=verdict,
      Concentration=round(C,4), C_null=round(mean(nullC),4), p_C=round(pC,4),
      tau_cyc=round(tc$tau,4), tau_null=round(mean(nullT),4), p_tau=round(pT,4),
      tau_dir=tc$dir, row.names=NULL)
  }
}

out <- do.call(rbind, rows)
write.csv(out, file.path(RESULT_DIR,"cyclic_ordering.csv"), row.names=FALSE)
cat("\n\n===================== SUMMARY =====================\n")
print(out[,c("Dataset","Method","CyclicOrder","Concentration","p_C","tau_cyc","p_tau")],
      row.names=FALSE)
cat("\nWritten to ", RESULT_DIR, "\n", sep="")
