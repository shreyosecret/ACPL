# ==============================================================================
# native_tau_ci.R   -- THE DECIDING RUN
#
# native_baselines.R answered the JCB objection using chance-corrected SW-SA.
# That is the wrong instrument: this project has established that the statistic
# is miscalibrated, direction-blind and insensitive. This script asks the same
# question with tau_cyc, the origin-free rank statistic, and puts paired
# bootstrap intervals on every comparison against ACPL.
#
# Five configurations per dataset:
#   ACPL       on 2D UMAP        (the method is defined on a planar embedding)
#   MST        on 2D UMAP  and on 10 PCs (native)
#   Slingshot  on 2D UMAP  and on 10 PCs (native)
#
#   Rscript analysis/native_tau_ci.R
# ==============================================================================

PROJ <- getwd()   # run from the repository root
CACHE_DIR <- file.path(PROJ,"cache"); RESULT_DIR <- file.path(PROJ,"results")
# ROTATION GRID. tau_cyc maximises over a grid of origin rotations, so its value
# depends on how fine that grid is. NROT = 24 here, and THIS SCRIPT IS THE SOURCE
# OF EVERY tau_cyc VALUE QUOTED IN THE MANUSCRIPTS. Script 02 uses 36 and script
# 03 uses 12, so their per-method tau values differ slightly from these. See
# "Why tau_cyc differs between scripts" in the README.
NPC <- 10; NROT <- 24; NBOOT <- 400
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
  den<-sqrt(n0*(n0-n2)); if(!is.finite(den)||den<=0) return(NA_real_)
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
rows<-list(); pw<-list()

for(tag in names(sets)){
  f<-file.path(CACHE_DIR,paste0("cache_",tag,".rds"))
  if(!file.exists(f)){cat("\n## missing",f,"\n"); next}
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

  cat("\n\n########", sets[[tag]], " (n =", n, ") ########\n")
  obs <- vapply(names(P), function(m) tau_from(P[[m]],pn,seq_len(n)), numeric(1))
  cat("  configuration              tau_cyc\n")
  for(m in names(obs)) cat(sprintf("  %-25s %7.4f\n", m, obs[[m]]))
  for(m in names(obs)) rows[[length(rows)+1]] <- data.frame(
    Dataset=sets[[tag]], Config=m, n=n, tau_cyc=round(obs[[m]],4), row.names=NULL)

  cat("\n  bootstrapping", NBOOT, "paired resamples ...\n")
  B <- matrix(NA_real_, NBOOT, length(P), dimnames=list(NULL,names(P)))
  for(b in seq_len(NBOOT)){
    idx <- sample.int(n,n,replace=TRUE)
    for(m in names(P)) B[b,m] <- tau_from(P[[m]],pn,idx)
  }

  cat("\n  ACPL vs                     diff     95% CI              crosses 0\n")
  for(m in setdiff(names(P),"ACPL (UMAP)")){
    d <- B[,"ACPL (UMAP)"] - B[,m]
    obsd <- obs[["ACPL (UMAP)"]] - obs[[m]]
    ci <- quantile(d,c(.025,.975),na.rm=TRUE)
    cross <- (ci[1]<=0 && ci[2]>=0)
    cat(sprintf("  %-25s %+7.4f  [%+.4f, %+.4f]    %s\n",
        m, obsd, ci[1], ci[2], if(cross) "YES" else "no"))
    pw[[length(pw)+1]] <- data.frame(
      Dataset=sets[[tag]], Comparison=paste("ACPL (UMAP) -",m), n=n,
      Diff=round(obsd,4), CI_lo=round(ci[1],4), CI_hi=round(ci[2],4),
      Crosses_zero=cross, row.names=NULL)
  }
}

t1<-do.call(rbind,rows); t2<-do.call(rbind,pw)
if (is.null(t1) || !nrow(t1)) stop("No rows produced: the dataset cache is missing or every method failed. Refusing to write an empty table.")
write.csv(t1, file.path(RESULT_DIR,"native_tau.csv"), row.names=FALSE)
if (is.null(t2) || !nrow(t2)) stop("No rows produced: the dataset cache is missing or every method failed. Refusing to write an empty table.")
write.csv(t2, file.path(RESULT_DIR,"native_tau_ci.csv"), row.names=FALSE)
cat("\n\n============== tau_cyc, NATIVE vs EMBEDDED ==============\n")
print(t1,row.names=FALSE)
cat("\n============== ACPL vs each, paired bootstrap ==============\n")
print(t2,row.names=FALSE)
cat("\nThe row that decides the paper is Nestorowa, ACPL vs Sling (PCA native).\n")
cat("Written to ", RESULT_DIR, "\n", sep="")
