# ==============================================================================
# tau_pairwise_ci.R
#
# Confidence intervals on DIFFERENCES in tau_cyc between methods.
# Without these, "ACPL leads on three of four datasets" is an unsupported claim.
#
# Cells are resampled with replacement; both methods are re-scored on the SAME
# resample, so the comparison stays paired. Percentile CI over NBOOT resamples.
#
# NROT reduced to 12 (from 36) and NBOOT to 500 to keep runtime tractable. The
# rotation grid shifts the point estimate slightly but affects both arms of
# every comparison identically.
#
#   Rscript analysis/tau_pairwise_ci.R
# ==============================================================================

PROJ <- getwd()   # run from the repository root
CACHE_DIR <- file.path(PROJ,"cache"); RESULT_DIR <- file.path(PROJ,"results")
# ROTATION GRID. 12 here, against the 24 used by analysis/05_native_tau_ci.R,
# which is the script the manuscripts quote: 500 bootstrap resamples on a finer
# grid is not tractable. The coarser grid shifts the point estimate slightly and
# shifts both arms of every comparison identically, so the DIFFERENCES and
# INTERVALS below are unaffected; only the per-method tau values move. Do not
# quote the point estimates from this script.
# See "Why tau_cyc differs between scripts" in the README.
NROT <- 12; NBOOT <- 500
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
cd_stat <- function(pn){
  b1<-c(0,head(cumsum(pn==1L),-1)); b2<-c(0,head(cumsum(pn==2L),-1)); b3<-c(0,head(cumsum(pn==3L),-1))
  sum(ifelse(pn==1L,0,ifelse(pn==2L,b1,b1+b2))) - sum(ifelse(pn==1L,b2+b3,ifelse(pn==2L,b3,0)))
}
tau_cyc <- function(pn, nrot=NROT){
  n<-length(pn); n0<-n*(n-1)/2
  n2<-sum(vapply(1:3,function(k){m<-sum(pn==k); m*(m-1)/2},numeric(1)))
  den<-sqrt(n0*(n0-n2)); if (!is.finite(den) || den<=0) return(NA_real_)
  cuts<-unique(round(seq(0,n-1,length.out=nrot+1)[1:nrot])); best<--Inf
  for (r in cuts){
    rolled <- if (r==0) pn else c(pn[(r+1):n], pn[1:r])
    for (s in list(rolled, rev(rolled))){
      v <- cd_stat(s)/den; if (v>best) best<-v
    }
  }
  best
}
tau_from <- function(t, pn, idx){
  tt<-t[idx]; pp<-pn[idx]; o<-order(tt); tau_cyc(pp[o])
}

sets <- c(nestorowa="Nestorowa HSC", buettner="Buettner mESC",
          leng="Leng hESC", richard="Richard T cells")
rows <- list()

for (tag in names(sets)){
  f<-file.path(CACHE_DIR,paste0("cache_",tag,".rds"))
  if(!file.exists(f)){cat("\n## missing",f,"\n"); next}
  o<-readRDS(f); cd<-o$coords; ph<-as.character(o$phases)
  pts<-list(ACPL=acpl_arc(cd))
  pts$MST       <- tryCatch(mst_pt(cd,ph),   error=function(e) NULL)
  pts$Slingshot <- tryCatch(sling_pt(cd,ph), error=function(e) NULL)
  pts<-pts[!vapply(pts,is.null,logical(1))]

  ok <- ph %in% names(pm)
  for (m in names(pts)) ok <- ok & !is.na(pts[[m]])
  pn <- pm[ph[ok]]; n <- sum(ok)
  P  <- lapply(pts, function(z) z[ok])

  cat("\n\n########", sets[[tag]], " (n =", n, ") ########\n")
  obs <- vapply(names(P), function(m) tau_from(P[[m]], pn, seq_len(n)), numeric(1))
  cat("  point estimates (NROT =", NROT, "):",
      paste(sprintf("%s=%.4f", names(obs), obs), collapse="  "), "\n")
  cat("  bootstrapping", NBOOT, "paired resamples ...\n")

  boots <- matrix(NA_real_, NBOOT, length(P), dimnames=list(NULL, names(P)))
  for (b in seq_len(NBOOT)){
    idx <- sample.int(n, n, replace=TRUE)
    for (m in names(P)) boots[b,m] <- tau_from(P[[m]], pn, idx)
  }

  cat("\n  comparison                  diff     95% CI              crosses 0\n")
  ms <- names(P)
  for (i in seq_along(ms)) for (j in seq_along(ms)) if (i<j){
    d    <- boots[,ms[i]] - boots[,ms[j]]
    obsd <- obs[ms[i]] - obs[ms[j]]
    ci   <- quantile(d, c(.025,.975), na.rm=TRUE)
    cross <- (ci[1] <= 0 && ci[2] >= 0)
    cat(sprintf("  %-25s %+7.4f  [%+.4f, %+.4f]    %s\n",
        paste(ms[i],"-",ms[j]), obsd, ci[1], ci[2], if (cross) "YES" else "no"))
    rows[[length(rows)+1]] <- data.frame(
      Dataset=sets[[tag]], Comparison=paste(ms[i],"-",ms[j]), n=n,
      Diff=round(obsd,4), CI_lo=round(ci[1],4), CI_hi=round(ci[2],4),
      Crosses_zero=cross, row.names=NULL)
  }
}

out <- do.call(rbind, rows)
if (is.null(out) || !nrow(out)) stop("No rows produced: the dataset cache is missing or every method failed. Refusing to write an empty table.")
write.csv(out, file.path(RESULT_DIR,"tau_pairwise_ci.csv"), row.names=FALSE)
cat("\n\n===================== SUMMARY =====================\n")
print(out, row.names=FALSE)
cat("\nAny comparison marked YES is not evidence of a difference.\n")
cat("Written to ", RESULT_DIR, "\n", sep="")
