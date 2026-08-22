# ==============================================================================
# 09_multiplicity.R -- family-wise adjustment for the sixteen ACPL comparisons
#
# analysis/05_native_tau_ci.R reports sixteen paired comparisons, four datasets
# by four baseline configurations, and seven of the sixteen intervals exclude
# zero. Reporting seven of sixteen unadjusted intervals as findings is the
# multiplicity objection the manuscript already concedes. This script prices it.
#
# WHAT IT DOES
#   1. Recomputes the same paired bootstrap on the same cache, with more
#      replicates (the tails matter here and 400 is too few).
#   2. Derives a two-sided bootstrap p-value per comparison.
#   3. Applies Holm across all sixteen. Holm is valid under arbitrary
#      dependence, which matters because the four comparisons within a dataset
#      share resamples and are strongly correlated. Westfall-Young step-down
#      would be more powerful but needs a joint resampling scheme across
#      datasets that does not exist here, so Holm is the honest choice.
#   4. Reports Bonferroni simultaneous intervals from the bootstrap standard
#      error rather than from extreme empirical quantiles. At m = 16 the
#      Bonferroni percentile endpoints sit at 0.156% and 99.84%, which even
#      2000 replicates estimate from about three points per tail. The normal
#      approximation is stable; the empirical quantile is not.
#
#   Rscript analysis/09_multiplicity.R
#
# RESOLUTION FLOOR. The +1-corrected bootstrap p-value cannot go below
# 2/(NBOOT+1), so the smallest attainable Holm-adjusted p is 16 * 2/(NBOOT+1).
# At NBOOT = 2000 that floor is 0.016: enough to claim significance at 0.05,
# not enough to claim anything stronger. Raise NBOOT to 3200 or more before
# writing "p < 0.01 adjusted" anywhere.
#
# Runtime scales with NBOOT and is dominated by Nestorowa. Budget 1-3 hours at
# the default. Set NBOOT lower for a smoke test, but do not report those numbers.
# ==============================================================================

PROJ <- getwd()
CACHE_DIR <- file.path(PROJ,"cache"); RESULT_DIR <- file.path(PROJ,"results")
NPC <- 10; NROT <- 24; NBOOT <- 2000; ALPHA <- 0.05
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
sling_pt <- function(cd, phases, rd, root="G1"){
  keep<-phases %in% c("G1","S","G2M")
  sds<-slingshot::slingshot(data=cd[keep,,drop=FALSE],clusterLabels=phases[keep],
        start.clus=root, reducedDim=rd)
  pt<-rep(NA_real_,nrow(cd)); pt[keep]<-rowMeans(slingshot::slingPseudotime(sds),na.rm=TRUE); pt
}
pm <- c(G1=1L,S=2L,G2M=3L)
cd_stat <- function(pn){
  b1<-c(0,head(cumsum(pn==1L),-1)); b2<-c(0,head(cumsum(pn==2L),-1)); b3<-c(0,head(cumsum(pn==3L),-1))
  sum(ifelse(pn==1L,0,ifelse(pn==2L,b1,b1+b2))) - sum(ifelse(pn==1L,b2+b3,ifelse(pn==2L,b3,0)))
}
tau_cyc <- function(pn, nrot=NROT){
  n<-length(pn); n0<-n*(n-1)/2
  n2<-sum(vapply(1:3,function(k){m<-sum(pn==k); m*(m-1)/2},numeric(1)))
  den<-sqrt(n0*(n0-n2)); if(!is.finite(den) || den<=0) return(NA_real_)
  cuts<-unique(round(seq(0,n-1,length.out=nrot+1)[1:nrot])); best<--Inf
  for(r in cuts){
    rolled <- if(r==0) pn else c(pn[(r+1):n], pn[1:r])
    for(s in list(rolled, rev(rolled))){ v<-cd_stat(s)/den; if(v>best) best<-v }
  }
  best
}
tau_from <- function(t,pn,idx){ tt<-t[idx]; pp<-pn[idx]; o<-order(tt); tau_cyc(pp[o]) }

sets <- c(nestorowa="Nestorowa HSC", buettner="Buettner mESC",
          leng="Leng hESC", richard="Richard T cells")
rows <- list()

for(tag in names(sets)){
  f<-file.path(CACHE_DIR,paste0("cache_",tag,".rds"))
  if(!file.exists(f)) stop("Missing cache: ", f, "\nRun run_all.R first.")
  o<-readRDS(f); um<-o$coords; ph<-as.character(o$phases)
  pc<-o$pca[, seq_len(min(NPC, ncol(o$pca))), drop=FALSE]

  P <- list()
  P[["ACPL (UMAP)"]]        <- acpl_arc(um)
  P[["MST (UMAP)"]]         <- tryCatch(mst_pt(um,ph), error=function(e) NULL)
  P[["MST (PCA native)"]]   <- tryCatch(mst_pt(pc,ph), error=function(e) NULL)
  P[["Sling (UMAP)"]]       <- tryCatch(sling_pt(um,ph,"UMAP"), error=function(e) NULL)
  P[["Sling (PCA native)"]] <- tryCatch(sling_pt(pc,ph,"PCA"),  error=function(e) NULL)
  P <- P[!vapply(P,is.null,logical(1))]

  ok <- ph %in% names(pm)
  for(m in names(P)) ok <- ok & !is.na(P[[m]])
  pn <- pm[ph[ok]]; n <- sum(ok); P <- lapply(P, function(z) z[ok])

  obs <- vapply(names(P), function(m) tau_from(P[[m]],pn,seq_len(n)), numeric(1))
  cat(sprintf("\n######## %s (n = %d) ########\n", sets[[tag]], n))
  cat("  bootstrapping", NBOOT, "paired resamples ...\n")

  B <- matrix(NA_real_, NBOOT, length(P), dimnames=list(NULL,names(P)))
  for(b in seq_len(NBOOT)){
    idx <- sample.int(n,n,replace=TRUE)
    for(m in names(P)) B[b,m] <- tau_from(P[[m]],pn,idx)
  }

  for(m in setdiff(names(P),"ACPL (UMAP)")){
    d    <- B[,"ACPL (UMAP)"] - B[,m]
    d    <- d[is.finite(d)]
    obsd <- obs[["ACPL (UMAP)"]] - obs[[m]]
    Bn   <- length(d)
    # two-sided bootstrap p-value, +1 correction so p is never exactly zero
    p_lo <- (1 + sum(d <= 0)) / (Bn + 1)
    p_hi <- (1 + sum(d >= 0)) / (Bn + 1)
    p    <- min(1, 2 * min(p_lo, p_hi))
    se   <- sd(d)
    rows[[length(rows)+1]] <- data.frame(
      Dataset=sets[[tag]], Comparison=paste("ACPL (UMAP) -",m), n=n, B=Bn,
      Diff=round(obsd,4),
      CI_lo=round(quantile(d,ALPHA/2,names=FALSE),4),
      CI_hi=round(quantile(d,1-ALPHA/2,names=FALSE),4),
      SE=round(se,5), p_boot=signif(p,4), row.names=NULL)
  }
}

out <- do.call(rbind, rows)
if (is.null(out) || !nrow(out)) stop("No comparisons produced. Refusing to write an empty table.")

m <- nrow(out)
out$p_holm     <- signif(p.adjust(out$p_boot, method = "holm"), 4)
out$Sig_raw    <- out$p_boot  < ALPHA
out$Sig_holm   <- out$p_holm  < ALPHA
zb             <- qnorm(1 - ALPHA / (2 * m))
out$Bonf_lo    <- round(out$Diff - zb * out$SE, 4)
out$Bonf_hi    <- round(out$Diff + zb * out$SE, 4)
out$Bonf_excl0 <- (out$Bonf_lo > 0) | (out$Bonf_hi < 0)

# cross-check the point estimates against the published table
ref <- file.path(RESULT_DIR, "native_tau_ci.csv")
if (file.exists(ref)) {
  r <- read.csv(ref, stringsAsFactors = FALSE)
  k <- merge(out[,c("Dataset","Comparison","Diff")],
             r[,c("Dataset","Comparison","Diff")],
             by=c("Dataset","Comparison"), suffixes=c("_new","_pub"))
  if (nrow(k)) {
    drift <- max(abs(k$Diff_new - k$Diff_pub))
    cat(sprintf("\nPoint-estimate check against native_tau_ci.csv: max |drift| = %.5f over %d rows\n",
                drift, nrow(k)))
    if (drift > 1e-4)
      cat("  WARNING: point estimates moved. The cache is not the one the table was built from.\n")
  }
}

write.csv(out, file.path(RESULT_DIR,"multiplicity_adjusted.csv"), row.names=FALSE)

cat("\n================ FAMILY-WISE ADJUSTMENT (m =", m, ") ================\n")
print(out[,c("Dataset","Comparison","Diff","CI_lo","CI_hi","p_boot","p_holm",
             "Sig_raw","Sig_holm","Bonf_lo","Bonf_hi","Bonf_excl0")], row.names=FALSE)
cat(sprintf("\nUnadjusted, %d of %d exclude zero.\n", sum(out$Sig_raw), m))
cat(sprintf("After Holm,  %d of %d survive at alpha = %.2f.\n", sum(out$Sig_holm), m, ALPHA))
cat(sprintf("Bonferroni simultaneous intervals exclude zero in %d of %d.\n",
            sum(out$Bonf_excl0), m))
if (sum(out$Sig_holm) < sum(out$Sig_raw))
  cat("\nThe difference between those counts is the claim the unadjusted table\n",
      "was overstating. Report the adjusted count.\n", sep="")
cat("\nWritten to ", file.path(RESULT_DIR,"multiplicity_adjusted.csv"), "\n", sep="")
