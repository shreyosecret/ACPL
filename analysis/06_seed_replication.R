# ==============================================================================
# seed_replication.R   -- the last threat to the Nestorowa result
#
# Every ACPL number so far rests on ONE UMAP embedding. UMAP is stochastic. A
# single draw producing a significant effect is exactly the shape of a result
# that evaporates on replication.
#
# This re-runs UMAP from the cached PCA under NSEED different seeds, recomputes
# ACPL each time, and compares against the PCA-based baselines, which do not
# depend on the UMAP seed at all. The comparison therefore isolates exactly the
# quantity in doubt.
#
#   Rscript analysis/seed_replication.R
# ==============================================================================

PROJ <- getwd()   # run from the repository root
CACHE_DIR <- file.path(PROJ,"cache"); RESULT_DIR <- file.path(PROJ,"results")
NSEED <- 10; NPC <- 10; NROT <- 24
suppressPackageStartupMessages({library(Seurat); library(igraph); library(uwot)})

# UMAP settings must match those used to build the cache
UMAP_PARAMS <- list(
  nestorowa = list(n_neighbors=30, min_dist=0.3, dims=10),
  buettner  = list(n_neighbors=15, min_dist=0.4, dims=10),
  leng      = list(n_neighbors=20, min_dist=0.4, dims=15),
  richard   = list(n_neighbors=30, min_dist=0.3, dims=15))

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
tau_of <- function(t, ph){
  ok <- ph %in% names(pm) & !is.na(t)
  tau_cyc(pm[ph[ok]][order(t[ok])])
}

sets <- c(nestorowa="Nestorowa HSC", buettner="Buettner mESC",
          leng="Leng hESC", richard="Richard T cells")
rows <- list()

for(tag in names(sets)){
  f<-file.path(CACHE_DIR,paste0("cache_",tag,".rds"))
  if(!file.exists(f)){cat("\n## missing",f,"\n"); next}
  o<-readRDS(f); ph<-as.character(o$phases); P<-o$pca
  prm <- UMAP_PARAMS[[tag]]
  pc_umap <- P[, seq_len(min(prm$dims, ncol(P))), drop=FALSE]   # input to UMAP
  pc_base <- P[, seq_len(min(NPC,      ncol(P))), drop=FALSE]   # input to baselines

  cat("\n\n########", sets[[tag]], "########\n")

  # baselines in native form: independent of the UMAP seed, computed once
  tS <- tryCatch(tau_of(sling_pt(pc_base, ph, "PCA"), ph), error=function(e) NA_real_)
  tM <- tryCatch(tau_of(mst_pt(pc_base, ph),          ph), error=function(e) NA_real_)
  cat(sprintf("  fixed native baselines:  Slingshot(PCA) %.4f   MST(PCA) %.4f\n", tS, tM))

  cat("  re-running UMAP under", NSEED, "seeds ...\n")
  cat("   seed   ACPL tau   ACPL - Sling(PCA)   ACPL - MST(PCA)\n")
  tA <- numeric(NSEED)
  for(s in seq_len(NSEED)){
    set.seed(s)
    um <- uwot::umap(pc_umap, n_neighbors=prm$n_neighbors, min_dist=prm$min_dist,
                     n_components=2, verbose=FALSE)
    tA[s] <- tau_of(acpl_arc(um), ph)
    cat(sprintf("   %4d   %8.4f   %17.4f   %15.4f\n", s, tA[s], tA[s]-tS, tA[s]-tM))
  }

  dS <- tA - tS; dM <- tA - tM
  cat(sprintf("\n  ACPL tau over seeds: mean %.4f  SD %.4f  range [%.4f, %.4f]\n",
              mean(tA), sd(tA), min(tA), max(tA)))
  cat(sprintf("  ACPL - Slingshot(PCA): mean %+.4f  SD %.4f  seeds where ACPL wins: %d/%d\n",
              mean(dS), sd(dS), sum(dS>0), NSEED))
  cat(sprintf("  ACPL - MST(PCA)      : mean %+.4f  SD %.4f  seeds where ACPL wins: %d/%d\n",
              mean(dM), sd(dM), sum(dM>0), NSEED))

  rows[[length(rows)+1]] <- data.frame(
    Dataset=sets[[tag]], Sling_PCA=round(tS,4), MST_PCA=round(tM,4),
    ACPL_mean=round(mean(tA),4), ACPL_sd=round(sd(tA),4),
    ACPL_min=round(min(tA),4), ACPL_max=round(max(tA),4),
    diff_Sling_mean=round(mean(dS),4), wins_vs_Sling=sum(dS>0),
    diff_MST_mean=round(mean(dM),4),  wins_vs_MST=sum(dM>0),
    n_seeds=NSEED, row.names=NULL)
}

out<-do.call(rbind,rows)
if (is.null(out) || !nrow(out)) stop("No rows produced: the dataset cache is missing or every method failed. Refusing to write an empty table.")
write.csv(out, file.path(RESULT_DIR,"seed_replication.csv"), row.names=FALSE)
cat("\n\n===================== SUMMARY =====================\n")
print(out,row.names=FALSE)
cat("\nIf ACPL wins on 10/10 seeds against native Slingshot on Nestorowa, the\n")
cat("result is not an artefact of one embedding. If it wins on 5-7, it is.\n")
cat("Written to ", RESULT_DIR, "\n", sep="")
