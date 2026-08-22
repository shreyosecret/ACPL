# ==============================================================================
# yeast_reference_free.R   -- does ACPL have a niche tricycle cannot fill?
#
# tricycle projects onto a reference learned from mouse neurosphere data. ACPL
# uses only the geometry of the data in front of it. If tricycle cannot run on
# a cyclic system outside its reference organism while ACPL recovers the cycle,
# that is a reason for ACPL to exist that does not depend on beating tricycle
# where tricycle works.
#
# Spellman alpha-factor synchronised yeast: 18 timepoints at 7 min intervals,
# cell cycle period roughly 66 min, so about 1.8 cycles. Ground truth is the
# EXPERIMENTAL TIMEPOINT, not an inferred phase label, which makes it a stronger
# reference than anything used so far in this project. Because ground truth is
# continuous, Jammalamadaka circular correlation is well defined here: the
# three-point degeneracy that broke it on phase labels does not arise.
#
# Install first if needed:  BiocManager::install("yeastCC")
#
#   Rscript analysis/yeast_reference_free.R
# ==============================================================================

PROJ <- getwd()   # run from the repository root
RESULT_DIR <- file.path(PROJ,"results")
PERIOD <- 66      # minutes; sensitivity swept below
NPERM  <- 5000
set.seed(2026)

for (p in c("yeastCC","uwot","igraph")) if (!requireNamespace(p, quietly=TRUE))
  stop("Install: BiocManager::install(\"", p, "\")")
suppressPackageStartupMessages({library(yeastCC); library(uwot); library(igraph)})

# ---------------------------------------------------------------- engine
find_adaptive_origin <- function(x,y)
  optim(c(median(x),median(y)), function(p) var(sqrt((x-p[1])^2+(y-p[2])^2)))$par
acpl_theta <- function(cd){
  uv <- find_adaptive_origin(cd[,1],cd[,2]); atan2(cd[,2]-uv[2], cd[,1]-uv[1])
}
mst_pt <- function(cd, root_idx){
  g<-igraph::graph_from_adjacency_matrix(as.matrix(dist(cd)),weighted=TRUE,mode="undirected")
  tr<-igraph::mst(g)
  as.numeric(igraph::distances(tr,v=root_idx,weights=igraph::E(tr)$weight))
}
circ_mean <- function(a) atan2(mean(sin(a)), mean(cos(a)))
circ_cor  <- function(a,b){ da<-a-circ_mean(a); db<-b-circ_mean(b)
  sum(sin(da)*sin(db))/sqrt(sum(sin(da)^2)*sum(sin(db)^2)) }
rank_angle <- function(v) 2*pi*(rank(v, ties.method="first")-0.5)/length(v)

# ---------------------------------------------------------------- data
# The object in current yeastCC is spYCCES, not "spellman". The original ACPL
# repository called data(spellman), which no longer resolves. Discover it rather
# than hardcode.
avail <- tryCatch(as.character(data(package="yeastCC")$results[,"Item"]),
                  error=function(e) character(0))
cat("objects in yeastCC: ", paste(avail, collapse=", "), "\n", sep="")
eset <- NULL
for (nm in c("spYCCES", avail)) {
  ok <- tryCatch({ data(list=nm, package="yeastCC", envir=environment()); TRUE },
                 error=function(e) FALSE, warning=function(w) FALSE)
  if (!ok) next
  obj <- get(nm, envir=environment())
  if (inherits(obj, "ExpressionSet")) { eset <- obj
    cat("using object: ", nm, "\n", sep=""); break }
}
if (is.null(eset)) stop("No ExpressionSet found in yeastCC.")

X <- Biobase::exprs(eset)
cat(sprintf("full matrix: %d genes x %d samples\n", nrow(X), ncol(X)))
cat("sample names (first 12): ", paste(head(colnames(X),12), collapse=", "), "\n", sep="")

alpha <- grep("^alpha", colnames(X), value=TRUE)
if (!length(alpha)) {
  cat("no 'alpha' prefix found; falling back to the cdc15 series\n")
  alpha <- grep("^cdc15", colnames(X), value=TRUE)
}
if (!length(alpha)) stop("Could not identify a synchronised series in the column names.")
X <- X[, alpha, drop=FALSE]
tmin <- as.numeric(gsub("[^0-9]", "", colnames(X)))
ord  <- order(tmin); X <- X[,ord,drop=FALSE]; tmin <- tmin[ord]

X <- X[rowMeans(is.na(X)) < 0.2, , drop=FALSE]
for (i in seq_len(nrow(X))) { na<-is.na(X[i,]); if(any(na)) X[i,na]<-mean(X[i,],na.rm=TRUE) }
cat(sprintf("selected series: %d genes x %d timepoints, t = %s min\n",
            nrow(X), ncol(X), paste(range(tmin), collapse=" to ")))

E  <- t(scale(t(X)))                      # genes x samples, gene-standardised
pc <- prcomp(t(E), scale.=FALSE)$x[, 1:min(5, ncol(X)-1), drop=FALSE]
um <- uwot::umap(pc, n_neighbors=min(5, ncol(X)-1), min_dist=0.5,
                 n_components=2, verbose=FALSE)

true_ang <- 2*pi * ((tmin %% PERIOD) / PERIOD)

# ---------------------------------------------------------------- tricycle
cat("\n================ tricycle on yeast ================\n")
tri_ok <- FALSE
if (requireNamespace("tricycle", quietly=TRUE)) {
  res <- tryCatch({
    sce <- SingleCellExperiment::SingleCellExperiment(
             list(logcounts = X))
    out <- tricycle::estimate_cycle_position(sce, gname.type="SYMBOL", species="mouse")
    tri_ok <- TRUE
    as.numeric(SummarizedExperiment::colData(out)$tricyclePosition)
  }, error=function(e){ cat("  tricycle FAILED, verbatim error:\n    ",
                           conditionMessage(e), "\n", sep=""); NULL })
  if (tri_ok) cat("  tricycle returned positions -- it did NOT fail. Report this honestly.\n")
} else cat("  tricycle not installed; install to complete the comparison.\n")
cat("  (rownames here are yeast systematic ORF names, e.g. ",
    paste(head(rownames(X),3), collapse=", "), ")\n", sep="")

# ---------------------------------------------------------------- scoring
# v may be an ANGLE (ACPL) or a scalar pseudotime (MST, Slingshot).
# Ranking an angle destroys it, which was the fault in the first version.
score <- function(v, label, is_angle = FALSE){
  a  <- if (is_angle) v else rank_angle(v)
  rc <- circ_cor(a, true_ang)
  nul <- replicate(NPERM, circ_cor(a, true_ang[sample(length(true_ang))]))
  p  <- (1 + sum(abs(nul) >= abs(rc)))/(NPERM+1)
  sp <- suppressWarnings(cor(rank(v), tmin, method="spearman"))
  cat(sprintf("  %-22s circ_cor %+7.4f  p %.4f   Spearman vs time %+7.4f\n",
              label, rc, p, sp))
  data.frame(Method=label, Circ_cor=round(rc,4), p=round(p,4),
             Spearman_time=round(sp,4), row.names=NULL)
}

cat("\n================ reference-free methods ================\n")
cat("  ground truth is experimental timepoint; period assumed", PERIOD, "min\n\n")
rows <- list()
rows[[1]] <- score(acpl_theta(um), "ACPL (UMAP)", is_angle = TRUE)
root <- which.min(tmin)
rows[[2]] <- score(mst_pt(um, root),             "MST (UMAP)")
rows[[3]] <- score(mst_pt(pc, root),             "MST (PCA)")
if (requireNamespace("slingshot", quietly=TRUE)) {
  cl <- kmeans(um, centers=min(4, ncol(X)-1), nstart=25)$cluster
  s <- tryCatch({
    sds <- slingshot::slingshot(data=um, clusterLabels=as.character(cl),
             start.clus=as.character(cl[root]), reducedDim="UMAP")
    rowMeans(slingshot::slingPseudotime(sds), na.rm=TRUE)
  }, error=function(e){ message("   slingshot failed: ", conditionMessage(e)); NULL })
  if(!is.null(s)) rows[[length(rows)+1]] <- score(s, "Slingshot (UMAP)")
}
if (tri_ok) rows[[length(rows)+1]] <- score(res, "tricycle")

# --------------------------------------------------------------------------
# The full series spans about 1.8 cycles. ACPL's planar polar construction
# cannot separate lap 1 from lap 2: every cell is placed on one circle. So the
# multi-lap series is outside the method's stated scope. Restrict to the first
# full cycle, where the comparison is meaningful.
cat("\n================ restricted to ONE cycle ================\n")
one <- tmin < PERIOD
cat(sprintf("  timepoints kept: %d of %d  (t = %s min)\n",
            sum(one), length(tmin), paste(tmin[one], collapse=", ")))
if (sum(one) >= 8) {
  E1  <- t(scale(t(X[, one, drop=FALSE])))
  pc1 <- prcomp(t(E1), scale.=FALSE)$x[, 1:min(4, sum(one)-1), drop=FALSE]
  um1 <- uwot::umap(pc1, n_neighbors=min(4, sum(one)-1), min_dist=0.5,
                    n_components=2, verbose=FALSE)
  t1  <- tmin[one]
  ta1 <- 2*pi*(t1/PERIOD)
  sc1 <- function(v, label, is_angle=FALSE){
    a  <- if (is_angle) v else 2*pi*(rank(v,ties.method="first")-0.5)/length(v)
    rc <- circ_cor(a, ta1)
    nul <- replicate(NPERM, circ_cor(a, ta1[sample(length(ta1))]))
    pv <- (1+sum(abs(nul)>=abs(rc)))/(NPERM+1)
    cat(sprintf("  %-22s circ_cor %+7.4f  p %.4f\n", label, rc, pv))
    data.frame(Method=paste0(label," [1 cycle]"), Circ_cor=round(rc,4),
               p=round(pv,4), Spearman_time=NA_real_, row.names=NULL)
  }
  r1 <- which.min(t1)
  rows[[length(rows)+1]] <- sc1(acpl_theta(um1), "ACPL (UMAP)", is_angle=TRUE)
  rows[[length(rows)+1]] <- sc1(mst_pt(um1, r1), "MST (UMAP)")
  rows[[length(rows)+1]] <- sc1(mst_pt(pc1, r1), "MST (PCA)")
} else cat("  too few timepoints in one cycle to test\n")

cat("\n================ period sensitivity ================\n")
cat("  ACPL circular correlation across assumed period:\n")
a_acpl <- rank_angle(acpl_theta(um))
for (P in seq(56, 76, by=2)) {
  ta <- 2*pi*((tmin %% P)/P)
  cat(sprintf("   period %2d min   circ_cor %+7.4f\n", P, circ_cor(a_acpl, ta)))
}

out <- do.call(rbind, rows)
if (is.null(out) || !nrow(out)) stop("No rows produced: the dataset cache is missing or every method failed. Refusing to write an empty table.")
write.csv(out, file.path(RESULT_DIR,"yeast_reference_free.csv"), row.names=FALSE)
cat("\n===================== SUMMARY =====================\n")
print(out, row.names=FALSE)
cat("\nThe question this answers: can tricycle run at all on a cyclic system\n")
cat("outside its reference organism, and does ACPL recover the cycle there?\n")
cat("Written to ", RESULT_DIR, "\n", sep="")
