# Called by support_compare.py with output, did source, and DRDID source paths.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3L)
p <- normalizePath(args[[1L]], mustWork=TRUE)
did_path <- normalizePath(args[[2L]], mustWork=TRUE)
drdid_path <- normalizePath(args[[3L]], mustWork=TRUE)
source_info <- function(path) {
  status <- system2("git", c("-C", shQuote(path), "status", "--porcelain", "--untracked-files=no"), stdout=TRUE)
  stopifnot(length(status) == 0L)
  list(commit=system2("git", c("-C", shQuote(path), "rev-parse", "HEAD"), stdout=TRUE), tracked_source_clean=TRUE)
}
sources <- list(did=source_info(did_path), DRDID=source_info(drdid_path))
suppressPackageStartupMessages(pkgload::load_all(drdid_path, quiet=TRUE))
suppressPackageStartupMessages(pkgload::load_all(did_path, quiet=TRUE))
stopifnot(normalizePath(getNamespaceInfo("did", "path")) == did_path,
          normalizePath(getNamespaceInfo("DRDID", "path")) == drdid_path)
sources$did$version <- as.character(packageVersion("did"))
sources$DRDID$version <- as.character(packageVersion("DRDID"))
jsonlite::write_json(sources, file.path(p,"r-sources.json"), auto_unbox=TRUE, pretty=TRUE)
data.table::setDTthreads(1L)
cases <- jsonlite::fromJSON(file.path(p,'cases.json'),simplifyVector=FALSE)
results <- list()
for (case in cases) {
 d <- read.csv(file.path(p,paste0(case$name,'.csv')))
 for (mode in c(FALSE,TRUE)) {
  fit <- suppressWarnings(suppressMessages(did::att_gt(yname='Y',tname='period',idname='id',gname='G',data=d,xformla=~X,
    weightsname='weight',faster_mode=mode,panel=case$route!='rc',allow_unbalanced_panel=case$route=='unbalanced',fix_weights=if(case$route=='varying')'varying' else NULL,
    est_method=case$method,bstrap=FALSE,cband=FALSE)))
  influence <- as.matrix(fit$inffunc)
  results[[paste0(case$name,'_',mode)]] <- list(att=fit$att,se=fit$se,
    if_ids=if(case$route=='rc')d$id[as.integer(rownames(influence))] else as.integer(rownames(influence)),
    if_shape=dim(influence),if_values=as.vector(t(influence)))
 }
}
jsonlite::write_json(results,file.path(p,'r-results.json'),auto_unbox=TRUE,digits=17,na='null')
