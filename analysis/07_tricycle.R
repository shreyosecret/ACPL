# ==============================================================================
# tricycle_benchmark.R   -- the missing baseline
#
# tricycle (Zheng et al., Genome Biology 2022) infers cell-cycle position as a
# polar angle by projecting onto a pre-learned reference. It is the method that
# owns this problem, and like ACPL it requires no starting cluster. Omitting it
# is the first thing a reviewer will ask about.
#
# NOTE: this reloads the raw datasets. The cache_*.rds files hold only the
# embedding, PCA and phase labels, and tricycle needs the expression matrix.
# First run is slow; results are cached to cache_tricycle_<tag>.rds.
#
# Install first if needed:
#   BiocManager::install("tricycle")
#
#   Rscript analysis/tricycle_benchmark.R
# ==============================================================================

PROJ <- getwd()   # run from the repository root
CACHE_DIR <- file.path(PROJ,"cache"); RESULT_DIR <- file.path(PROJ,"results")
NROT <- 24; NBOOT <- 400
set.seed(2026)

if (!requireNamespace("tricycle", quietly=TRUE))
  stop("Install first:  BiocManager::install(\"tricycle\")")
suppressPackageStartupMessages({
  library(Seurat); library(igraph); library(scRNAseq)
  library(SummarizedExperiment); library(SingleCellExperiment); library(tricycle)})

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

find_adaptive_origin <- function(x,y)
  optim(c(median(x),median(y)), function(p) var(sqrt((x-p[1])^2+(y-p[2])^2)))$par
unwrap_theta <- function(th){d<-diff(th);d<-d-2*pi*round(d/(2*pi));c(th[1],th[1]+cumsum(d))}
acpl_arc <- function(cd, span=0.25){
  uv<-find_adaptive_origin(cd[,1],cd[,2]); th<-atan2(cd[,2]-uv[2],cd[,1]-uv[1])
  o<-order(th); tu<-unwrap_theta(th[o])
  arc<-c(0,cumsum(sqrt(diff(cos(tu))^2+diff(sin(tu))^2)))
  predict(loess(arc~tu,span=span))[order(seq_along(cd[,1])[o])]
}

# Ensembl -> symbol, same mapping used to build the cache
ens_to_symbol <- function(sce, db_name="org.Mm.eg.db"){
  if(!requireNamespace(db_name, quietly=TRUE)) stop("need ", db_name)
  db  <- getExportedValue(db_name, db_name)
  ids <- sub("\\..*$","",rownames(sce))
  sym <- suppressMessages(AnnotationDbi::mapIds(db, keys=ids, column="SYMBOL",
           keytype="ENSEMBL", multiVals="first"))
  keep <- !is.na(sym) & !duplicated(sym)
  sce <- sce[keep,]; rownames(sce) <- as.character(sym[keep]); sce
}
pick_assay <- function(sce){
  an<-assayNames(sce); for(a in c("counts","normcounts","logcounts","tpm")) if(a %in% an) return(a); an[1]
}

# species and loader per dataset
# Each loader must return cells in EXACTLY the order the cache stored phases in,
# otherwise tricycle positions get paired with the wrong labels.
SPEC <- list(
  nestorowa = list(sp="mouse", fn=function() ens_to_symbol(NestorowaHSCData(),"org.Mm.eg.db")),
  buettner  = list(sp="mouse", fn=function(){
      sce <- BuettnerESCData()
      sce <- sce[, !is.na(colData(sce)$phase)]        # cache kept only labelled cells
      if (!any(grepl("^ENS", head(rownames(sce)))))   # already symbols?
        return(sce)
      ens_to_symbol(sce, "org.Mm.eg.db") }),
  leng      = list(sp="human", fn=function(){
      sce <- LengESCData()
      raw <- as.character(colData(sce)$Phase)
      keep <- raw %in% c("G1","S","G2","G2M")         # cache kept 247 of 460
      sce[, keep] }),
  richard   = list(sp="mouse", fn=function() ens_to_symbol(RichardTCellData(),"org.Mm.eg.db")))
sets <- c(nestorowa="Nestorowa HSC", buettner="Buettner mESC",
          leng="Leng hESC", richard="Richard T cells")

rows<-list(); pw<-list()

for(tag in names(sets)){
  cf <- file.path(CACHE_DIR, paste0("cache_",tag,".rds"))
  if(!file.exists(cf)){ cat("\n## missing",cf,"\n"); next }
  base <- readRDS(cf)

  tf <- file.path(CACHE_DIR, paste0("cache_tricycle_",tag,".rds"))
  if(file.exists(tf)){
    tri <- readRDS(tf); cat("\n[cache] tricycle position for", sets[[tag]], "\n")
  } else {
    cat("\n[compute] tricycle for", sets[[tag]], "-- reloading raw data\n")
    tri <- tryCatch({
      sce <- SPEC[[tag]]$fn()
      assay(sce,"counts") <- as.matrix(assay(sce, pick_assay(sce)))
      sce <- scuttle::logNormCounts(sce)
      sce <- tricycle::estimate_cycle_position(sce, gname.type="SYMBOL",
                                               species=SPEC[[tag]]$sp)
      as.numeric(colData(sce)$tricyclePosition)
    }, error=function(e){ message("   tricycle failed: ",conditionMessage(e)); NULL })
    if(!is.null(tri)) saveRDS(tri, tf)
  }
  if(is.null(tri)){ cat("  skipped\n"); next }

  ph <- as.character(base$phases)
  if(length(tri) != length(ph)){
    cat(sprintf("  LENGTH MISMATCH: tricycle %d vs cached phases %d -- skipped\n",
                length(tri), length(ph)))
    cat("  (the cache dropped cells during phase scoring; align before using)\n"); next
  }

  P <- list(ACPL = acpl_arc(base$coords), tricycle = tri)
  ok <- ph %in% names(pm)
  for(m in names(P)) ok <- ok & !is.na(P[[m]])
  pn <- pm[ph[ok]]; n <- sum(ok); P <- lapply(P, function(z) z[ok])

  cat(sprintf("\n######## %s  (n = %d) ########\n", sets[[tag]], n))
  obs <- vapply(names(P), function(m) tau_from(P[[m]],pn,seq_len(n)), numeric(1))

  # ---- the same orderings scored by the OLD statistic ----------------------
  # If windowed concordance cannot see a difference that tau_cyc measures at
  # 0.23, that is the sharpest demonstration of its insensitivity available.
  swsa <- function(t, phv, w=10){
    d <- data.frame(t=t, p=phv); d <- d[!is.na(d$t) & d$p %in% names(pm),,drop=FALSE]
    d$pn <- pm[d$p]; d <- d[order(d$t),,drop=FALSE]; nn <- nrow(d)
    100*mean(d$pn[1:(nn-w)] <= d$pn[(w+1):nn])
  }
  pp <- prop.table(table(factor(ph[ok], levels=names(pm))))
  chance <- 100*(1+sum(pp^2))/2
  cat(sprintf("  chance level for windowed concordance = %.2f%%\n", chance))
  cat("  method       tau_cyc    SW-SA   SW-SA corrected\n")
  for(m in names(obs)){
    sv <- swsa(P[[m]], ph[ok])
    cat(sprintf("  %-12s %7.4f %8.2f %14.2f\n",
                m, obs[[m]], sv, 100*(sv-chance)/(100-chance)))
  }
  for(m in names(obs)) rows[[length(rows)+1]] <- data.frame(
    Dataset=sets[[tag]], Method=m, n=n, tau_cyc=round(obs[[m]],4), row.names=NULL)

  B <- matrix(NA_real_, NBOOT, 2, dimnames=list(NULL,names(P)))
  for(b in seq_len(NBOOT)){
    idx <- sample.int(n,n,replace=TRUE)
    for(m in names(P)) B[b,m] <- tau_from(P[[m]],pn,idx)
  }
  d <- B[,"ACPL"] - B[,"tricycle"]; obsd <- obs[["ACPL"]] - obs[["tricycle"]]
  ci <- quantile(d,c(.025,.975),na.rm=TRUE)
  cat(sprintf("  ACPL - tricycle = %+.4f  95%% CI [%+.4f, %+.4f]  %s\n",
      obsd, ci[1], ci[2], if(ci[1]<=0 && ci[2]>=0) "crosses zero" else "excludes zero"))
  pw[[length(pw)+1]] <- data.frame(Dataset=sets[[tag]], n=n, Diff=round(obsd,4),
    CI_lo=round(ci[1],4), CI_hi=round(ci[2],4),
    Crosses_zero=(ci[1]<=0 && ci[2]>=0), row.names=NULL)
}

t1<-do.call(rbind,rows); t2<-do.call(rbind,pw)
if(!is.null(t1)) write.csv(t1, file.path(RESULT_DIR,"tricycle_tau.csv"), row.names=FALSE)
if(!is.null(t2)) write.csv(t2, file.path(RESULT_DIR,"tricycle_ci.csv"), row.names=FALSE)
cat("\n\n===================== SUMMARY =====================\n")
if(!is.null(t1)) print(t1,row.names=FALSE)
cat("\n")
if(!is.null(t2)) print(t2,row.names=FALSE)
cat("\nWritten to ", RESULT_DIR, "\n", sep="")
