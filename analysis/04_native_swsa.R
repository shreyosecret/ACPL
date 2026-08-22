# ==============================================================================
# native_baselines.R
#
# Answers the JCB desk-rejection objection: baselines must run on their NATIVE
# input, not on this method's UMAP embedding.
#
# For each dataset, Slingshot and Monocle-style MST are run twice:
#   (a) on the 2D UMAP  -- the handicapped configuration used in the preprint
#   (b) on the first 10 PCs -- native form
# ACPL runs on the 2D UMAP only; it is defined on a planar embedding.
#
# Scored with chance-corrected windowed concordance AND circular correlation.
#
#   Rscript analysis/native_baselines.R
# ==============================================================================

PROJ <- getwd()   # run from the repository root
CACHE_DIR  <- file.path(PROJ, "cache")
RESULT_DIR <- file.path(PROJ, "results")
NPC <- 10; WINDOW <- 10; NPERM <- 2000
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
sling_pt <- function(cd, phases, rd="PCA", root="G1"){
  keep <- phases %in% c("G1","S","G2M")
  sds <- slingshot::slingshot(data=cd[keep,,drop=FALSE],
          clusterLabels=phases[keep], start.clus=root, reducedDim=rd)
  pt <- rep(NA_real_, nrow(cd))
  pt[keep] <- rowMeans(slingshot::slingPseudotime(sds), na.rm=TRUE); pt
}

pm <- c(G1=1L,S=2L,G2M=3L)
chance <- function(ph){p<-prop.table(table(factor(ph,levels=names(pm))));100*(1+sum(p^2))/2}
swsa <- function(t,ph,w=WINDOW){
  d<-data.frame(t=t,p=ph); d<-d[!is.na(d$t)&d$p %in% names(pm),,drop=FALSE]
  d$pn<-pm[d$p]; d<-d[order(d$t),,drop=FALSE]; n<-nrow(d)
  100*mean(d$pn[1:(n-w)] <= d$pn[(w+1):n])
}
circ_mean <- function(a) atan2(mean(sin(a)),mean(cos(a)))
circ_cor <- function(a,b){da<-a-circ_mean(a);db<-b-circ_mean(b)
  sum(sin(da)*sin(db))/sqrt(sum(sin(da)^2)*sum(sin(db)^2))}
t_ang <- function(t) 2*pi*(rank(t,ties.method="average")-0.5)/length(t)
p_ang <- function(pn) c(0,2*pi/3,4*pi/3)[pn]

score <- function(t, phases, label, ds, cfg){
  ok <- phases %in% names(pm) & !is.na(t)
  ch <- chance(phases[ok]); s <- swsa(t[ok], phases[ok])
  ta <- t_ang(t[ok]); pa <- p_ang(pm[phases[ok]])
  rc <- circ_cor(ta, pa)
  nul <- replicate(NPERM, circ_cor(ta, pa[sample(length(pa))]))
  p <- (1+sum(abs(nul)>=abs(rc)))/(NPERM+1)
  cat(sprintf("  %-12s %-10s %7.2f %8.2f %+9.4f %8.4f\n",
              label, cfg, s, 100*(s-ch)/(100-ch), rc, p))
  data.frame(Dataset=ds, Method=label, Input=cfg, Chance=round(ch,2),
             SWSA=round(s,2), Corrected=round(100*(s-ch)/(100-ch),2),
             Circ_cor=round(rc,4), p_circ=round(p,4))
}

sets <- c(nestorowa="Nestorowa HSC", buettner="Buettner mESC",
          leng="Leng hESC", richard="Richard T cells")
rows <- list()

for (tag in names(sets)) {
  f <- file.path(CACHE_DIR, paste0("cache_",tag,".rds"))
  if (!file.exists(f)) { cat("\n## missing",f,"\n"); next }
  o <- readRDS(f); um <- o$coords; ph <- as.character(o$phases)
  pc <- o$pca[, seq_len(min(NPC, ncol(o$pca))), drop=FALSE]

  cat("\n\n########", sets[[tag]], "(N =", length(ph), ") ########\n")
  cat("  method       input       SW-SA corrected     rcirc   p(circ)\n")

  rows[[length(rows)+1]] <- score(acpl_arc(um), ph, "ACPL", sets[[tag]], "UMAP-2D")

  rows[[length(rows)+1]] <- score(mst_pt(um, ph), ph, "MST", sets[[tag]], "UMAP-2D")
  rows[[length(rows)+1]] <- score(mst_pt(pc, ph), ph, "MST", sets[[tag]],
                                  paste0("PCA-", ncol(pc), " (native)"))

  s1 <- tryCatch(sling_pt(um, ph, "UMAP"), error=function(e){message("   sling UMAP: ",conditionMessage(e)); NULL})
  if (!is.null(s1)) rows[[length(rows)+1]] <- score(s1, ph, "Slingshot", sets[[tag]], "UMAP-2D")
  s2 <- tryCatch(sling_pt(pc, ph, "PCA"),  error=function(e){message("   sling PCA: ",conditionMessage(e)); NULL})
  if (!is.null(s2)) rows[[length(rows)+1]] <- score(s2, ph, "Slingshot", sets[[tag]],
                                  paste0("PCA-", ncol(pc), " (native)"))
}

out <- do.call(rbind, rows)
if (is.null(out) || !nrow(out)) stop("No rows produced: the dataset cache is missing or every method failed. Refusing to write an empty table.")
write.csv(out, file.path(RESULT_DIR,"native_baselines.csv"), row.names=FALSE)
cat("\n\n===================== NATIVE vs EMBEDDED =====================\n")
print(out, row.names=FALSE)
cat("\nIf a baseline improves in native form, the preprint's comparison was\n")
cat("handicapping it, which is exactly the JCB objection.\n")
cat("\nWritten to ", RESULT_DIR, "\n", sep="")
