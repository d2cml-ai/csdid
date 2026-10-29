# Generate numerical references from pinned, clean R sources by default.
# Rscript generate.R --did /path/to/did --drdid /path/to/DRDID
# Deliberate upstream comparison: add --source-lock /path/to/source-lock.json.
# Existing inputs are reused unless --regenerate-inputs is supplied.
args <- commandArgs(trailingOnly = TRUE)
flag <- function(name, default = NULL) {
  at <- match(paste0("--", name), args)
  if (is.na(at)) return(default)
  if (at == length(args)) stop("Missing value for --", name)
  args[at + 1L]
}
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
output <- flag("output", dirname(script))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output)
did_path <- normalizePath(flag("did"), mustWork = TRUE)
drdid_path <- normalizePath(flag("drdid"), mustWork = TRUE)
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
git <- function(path, ...) system2("git", c("-C", shQuote(path), ...), stdout = TRUE)
expected_sources <- list(
  did = list(commit = "48863b646e657f2aaa5cc9e9d2e7b9d764a6a44e", version = "2.5.1.904"),
  DRDID = list(commit = "eaac3b49cd225d0d7706c633858555a5f204e21e", version = "1.4.0"))
valid_object <- function(value, fields) {
  is.list(value) && !is.null(names(value)) && !anyDuplicated(names(value)) &&
    length(value) == length(fields) && setequal(names(value), fields)
}
source_lock <- flag("source-lock")
if (!is.null(source_lock)) {
  lock <- jsonlite::fromJSON(normalizePath(source_lock, mustWork = TRUE), simplifyVector = FALSE)
  if (!valid_object(lock, c("schema_version", "sources")) ||
      !is.integer(lock$schema_version) || length(lock$schema_version) != 1L ||
      is.na(lock$schema_version) || lock$schema_version != 1 ||
      !valid_object(lock$sources, c("did", "DRDID"))) {
    stop("Invalid source lock: expected schema_version 1 and exactly did/DRDID sources.")
  }
  for (package in c("did", "DRDID")) {
    entry <- lock$sources[[package]]
    if (!valid_object(entry, c("commit", "version")) ||
        !is.character(entry$commit) || length(entry$commit) != 1L || is.na(entry$commit) ||
        !grepl("^[0-9a-fA-F]{40}$", entry$commit) ||
        !is.character(entry$version) || length(entry$version) != 1L || is.na(entry$version) ||
        !grepl("^[0-9]+([.-][0-9]+)+$", entry$version)) {
      stop("Invalid source lock entry for ", package, ": expected a 40-hex commit and R package version.")
    }
    entry$commit <- tolower(entry$commit)
    expected_sources[[package]] <- entry
  }
}
source_info <- function(path, package) {
  head <- git(path, "rev-parse", "HEAD")
  status <- git(path, "status", "--porcelain", "--untracked-files=no")
  description <- read.dcf(file.path(path, "DESCRIPTION"), fields = c("Package", "Version"))
  stopifnot(identical(head, expected_sources[[package]]$commit), length(status) == 0L,
            nrow(description) == 1L, description[1L, "Package"] == package,
            description[1L, "Version"] == expected_sources[[package]]$version)
  files <- git(path, "ls-files", "DESCRIPTION", "NAMESPACE", "R", "src")
  files <- files[file.exists(file.path(path, files)) & !dir.exists(file.path(path, files))]
  hashes <- setNames(lapply(file.path(path, files), sha256), files)
  list(commit = head, tracked_source_clean = TRUE, files_sha256 = hashes,
       source_tree_sha256 = digest::digest(paste(files, unlist(hashes), collapse = "\n"),
                                          algo = "sha256", serialize = FALSE))
}
sources <- list(did = source_info(did_path, "did"), DRDID = source_info(drdid_path, "DRDID"))
suppressPackageStartupMessages(pkgload::load_all(drdid_path, quiet = TRUE))
suppressPackageStartupMessages(pkgload::load_all(did_path, quiet = TRUE))
stopifnot(normalizePath(getNamespaceInfo("did", "path")) == did_path,
          normalizePath(getNamespaceInfo("DRDID", "path")) == drdid_path)
for (package in c("did", "DRDID")) {
  stopifnot(identical(as.character(getNamespaceVersion(package)),
                     as.character(package_version(expected_sources[[package]]$version))))
}
data.table::setDTthreads(1L)

# These choices are fixed before numerical comparison. No input is selected or
# changed in response to a numerical failure. CSVs are the shared authority.
input_files <- c("panel.csv", "rc.csv", "unbalanced.csv",
                 "multipliers_unit.csv", "multipliers_cluster.csv")
regenerate_inputs <- "--regenerate-inputs" %in% args
missing_inputs <- input_files[!file.exists(file.path(output, input_files))]
if (length(missing_inputs) && !regenerate_inputs) {
  stop("Missing frozen inputs: ", paste(missing_inputs, collapse = ", "),
       ". Copy the committed input files to the output directory, or explicitly ",
       "use --regenerate-inputs to recreate all inputs from the fixed specification.")
}
if (regenerate_inputs) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20261023L)
  n <- 240L
  ids <- seq_len(n)
  groups <- rep(c(0L, 3L, 4L, 5L), each = n / 4L)
  x <- tanh(rnorm(n, sd = .7)) + .035 * groups
  alpha <- rnorm(n, sd = .4)
  panel <- do.call(rbind, lapply(seq_len(5L), function(t) {
    xt <- x + .08 * sin(ids / 9 + t)
    effect <- ifelse(groups > 0L & t >= groups,
                     .6 + .16 * (t - groups) + .09 * xt, 0)
    data.frame(id = ids, period = t, G = groups, X = xt,
      Y = alpha + .27 * xt + .13 * t + .06 * t * xt + rnorm(n, sd = .45) + effect,
      weight = exp(.2 * x + .12 * sin(t + ids / 7)),
      cluster = 1L + (ids %% 17L))
  }))
  panel <- panel[sample.int(nrow(panel)), ]
  unbalanced <- panel[(panel$id * 7L + panel$period * 3L) %% 17L != 0L, ]
  rc <- panel
  rc$id <- 100000L + seq_len(nrow(rc))
  rc$X <- tanh(rnorm(nrow(rc), sd = .7)) + .035 * rc$G + .04 * rc$period
  rc$Y <- .27 * rc$X + .13 * rc$period + .06 * rc$period * rc$X +
    rnorm(nrow(rc), sd = .6) + ifelse(rc$G > 0L & rc$period >= rc$G,
      .6 + .16 * (rc$period - rc$G) + .09 * rc$X, 0)
  rc$weight <- exp(.2 * rc$X + .12 * sin(rc$period + rc$id / 7))
  rc$cluster <- 1L + (rc$id %% 17L)
  for (kind in c("panel", "rc", "unbalanced")) {
    write.csv(get(kind), file.path(output, paste0(kind, ".csv")), row.names = FALSE)
  }
  # Explicit, frozen N(0,1) multipliers, reused for every fit/aggregation call.
  # This tests identical draws, not equivalence of the R and NumPy RNGs.
  for (kind in c("unit", "cluster")) {
    keys <- if (kind == "unit") ids else sort(unique(panel$cluster))
    u <- matrix(rnorm(99L * length(keys)), nrow = 99L)
    colnames(u) <- as.character(keys)
    write.csv(u, file.path(output, paste0("multipliers_", kind, ".csv")), row.names = FALSE)
  }
}

make_case <- function(design, method, setting, fix_weights = NULL, extra = NULL,
                      bootstrap = FALSE, clustered = FALSE) {
  universal <- setting == "notyet_universal"
  name <- paste(c(design, method, setting, extra), collapse = "_")
  ctor <- list(yname = "Y", tname = "period", idname = "id", gname = "G",
    xformla = "~ X", panel = design != "rc", allow_unbalanced_panel = design == "unbalanced",
    weights_name = "weight", control_group = if (universal) "notyettreated" else "nevertreated",
    anticipation = 0L, cband = bootstrap, biters = 99L, alp = .05,
    compute_inffunc = TRUE, fix_weights = fix_weights,
    clustervar = if (clustered) "cluster" else NULL, pl = FALSE, cores = 1L)
  fit <- list(est_method = method, base_period = if (universal) "universal" else "varying",
              bstrap = bootstrap)
  aggregation <- list(bstrap = bootstrap, cband = bootstrap, na_rm = FALSE)
  r_ctor <- ctor
  names(r_ctor)[names(r_ctor) == "weights_name"] <- "weightsname"
  names(r_ctor)[names(r_ctor) == "clustervar"] <- "clustervars"
  list(name = name, data = paste0(design, ".csv"), constructor = ctor,
    fit = fit, aggregation = aggregation,
    bootstrap_draws = if (bootstrap) list(file = paste0("multipliers_", if (clustered) "cluster" else "unit", ".csv"),
      kind = if (clustered) "cluster" else "unit", biters = 99L,
      column_identity = if (clustered) "cluster" else "id") else NULL,
    expected_r_kwargs = c(r_ctor, fit, list(faster_mode = TRUE)),
    expected_r_aggregation_kwargs = list(bstrap = bootstrap, cband = bootstrap, na.rm = FALSE))
}
cases <- list()
for (design in c("panel", "rc", "unbalanced")) {
  for (method in c("dr", "reg", "ipw")) {
    fixed <- if (design == "panel") switch(method, dr = NULL, reg = "first_period", ipw = "base_period") else
      if (design == "unbalanced") switch(method, dr = "first_period", reg = "base_period", ipw = NULL) else NULL
    for (setting in c("never_varying", "notyet_universal")) {
      cases[[length(cases) + 1L]] <- make_case(design, method, setting, fixed)
    }
  }
}
cases[[19L]] <- make_case("panel", "dr", "never_varying", "varying", "varying_weights")
cases[[20L]] <- make_case("panel", "reg", "notyet_universal", "varying", "varying_weights")
cases[[21L]] <- make_case("panel", "dr", "never_varying", NULL, "bootstrap_unit", TRUE, FALSE)
cases[[22L]] <- make_case("panel", "reg", "notyet_universal", "base_period", "bootstrap_cluster", TRUE, TRUE)
cases[[23L]] <- make_case("panel", "dr", "never_varying", NULL, "analytic_cluster", FALSE, TRUE)
cases[[24L]] <- make_case("rc", "reg", "never_varying", NULL, "analytic_cluster", FALSE, TRUE)
cases[[25L]] <- make_case("unbalanced", "ipw", "never_varying", NULL, "analytic_cluster", FALSE, TRUE)
jsonlite::write_json(cases, file.path(output, "cases.json"), pretty = TRUE, digits = NA,
                    auto_unbox = TRUE, null = "null")

# Explicit dimensions and zero-based nonfinite masks prevent null/NaN or vector
# coercion from silently dropping comparisons. Matrix values are row-major.
numeric_payload <- function(value, influence = FALSE) {
  if (is.null(value)) return(NULL)
  if (influence) value <- as.matrix(value)
  shape <- if (is.null(dim(value))) length(value) else dim(value)
  flat <- if (length(shape) == 2L) as.numeric(t(value)) else as.numeric(value)
  mask <- list(nan = as.list(which(is.na(flat)) - 1L),
               posinf = as.list(which(is.infinite(flat) & flat > 0) - 1L),
               neginf = as.list(which(is.infinite(flat) & flat < 0) - 1L))
  flat[!is.finite(flat)] <- 0
  c(list(shape = as.list(as.integer(shape)), values = as.list(flat)), mask)
}

ns <- asNamespace("did")
original_bootstrap <- get("run_multiplier_bootstrap", envir = ns)
replace_bootstrap <- function(value) {
  unlockBinding("run_multiplier_bootstrap", ns)
  assign("run_multiplier_bootstrap", value, envir = ns)
  lockBinding("run_multiplier_bootstrap", ns)
}
active <- new.env(parent = emptyenv())
shared_bootstrap <- function(inf.func, biters, pl = FALSE, cores = 1L) {
  inf.func <- as.matrix(inf.func)
  n <- nrow(inf.func)
  stopifnot(n == ncol(active$draws), biters <= nrow(active$draws))
  ids <- suppressWarnings(as.numeric(rownames(inf.func)))
  if (active$kind == "unit") {
    if (is.null(active$permutation)) {
      stopifnot(length(ids) == n, all(is.finite(ids)), !anyDuplicated(ids),
                identical(sort(ids), active$draw_ids))
      active$permutation <- match(active$draw_ids, ids)
    }
    ord <- active$permutation
  } else {
    stopifnot(length(ids) == n, all(is.finite(ids)), !anyDuplicated(ids),
              identical(sort(ids), active$draw_ids))
    ord <- match(active$draw_ids, ids)
  }
  active$calls <- active$calls + 1L
  (active$draws[seq_len(biters), , drop = FALSE] %*% inf.func[ord, , drop = FALSE]) / n
}

execute <- function(z, fast) {
  d <- read.csv(file.path(output, z$data))
  if (!is.null(z$bootstrap_draws)) {
    active$draws <- as.matrix(read.csv(file.path(output, z$bootstrap_draws$file), check.names = FALSE))
    active$draw_ids <- as.numeric(colnames(active$draws))
    active$kind <- z$bootstrap_draws$kind
    active$permutation <- NULL
    active$calls <- 0L
    replace_bootstrap(shared_bootstrap)
  } else replace_bootstrap(original_bootstrap)
  on.exit(replace_bootstrap(original_bootstrap), add = TRUE)
  kw <- z$expected_r_kwargs
  kw$xformla <- as.formula(kw$xformla)
  kw$data <- d
  kw$faster_mode <- fast
  warnings <- character()
  result <- withCallingHandlers({
    fit <- do.call(did::att_gt, kw)
    ids <- suppressWarnings(as.numeric(rownames(fit$inffunc)))
    stopifnot(length(ids) == nrow(fit$inffunc), all(is.finite(ids)), !anyDuplicated(ids))
    if (!z$constructor$panel) {
      # R creates .rowid after filtering early-treated/missing observations.
      # It is not generally an index into the original shared CSV. Recover the
      # original unique observation ID from the retained authoritative frame.
      retained <- fit$DIDparams$data
      internal_ids <- retained[[fit$DIDparams$idname]]
      original_ids <- retained[[z$constructor$idname]]
      stopifnot(!anyDuplicated(internal_ids), !anyDuplicated(original_ids))
      position <- match(ids, internal_ids)
      stopifnot(!anyNA(position))
      ids <- as.numeric(original_ids[position])
    }
    ord <- order(ids)
    aggs <- setNames(lapply(c("simple", "group", "dynamic", "calendar"), function(kind) {
      a <- do.call(did::aggte, c(list(MP = fit, type = kind), z$expected_r_aggregation_kwargs))
      influence <- lapply(a$inf.function, function(v) {
        m <- as.matrix(v)
        stopifnot(nrow(m) == length(ids))
        numeric_payload(m[ord, , drop = FALSE], influence = TRUE)
      })
      names(influence) <- gsub(".", "_", names(influence), fixed = TRUE)
      list(overall_att = numeric_payload(a$overall.att), overall_se = numeric_payload(a$overall.se),
        egt = numeric_payload(a$egt), att_egt = numeric_payload(a$att.egt),
        se_egt = numeric_payload(a$se.egt), crit_val_egt = numeric_payload(a$crit.val.egt),
        inf_function = influence)
    }), c("simple", "group", "dynamic", "calendar"))
    list(unit_ids = as.list(ids[ord]), fit = list(group = numeric_payload(fit$group),
      t = numeric_payload(fit$t), att = numeric_payload(fit$att), se = numeric_payload(fit$se),
      c = numeric_payload(fit$c), l_se = numeric_payload(fit$att - fit$c * fit$se),
      u_se = numeric_payload(fit$att + fit$c * fit$se), W = numeric_payload(fit$W),
      Wpval = numeric_payload(fit$Wpval),
      inffunc = numeric_payload(as.matrix(fit$inffunc)[ord, , drop = FALSE], influence = TRUE),
      cband = isTRUE(fit$DIDparams$cband)), aggs = aggs)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  result$warnings <- as.list(warnings)
  result$bootstrap_calls <- if (is.null(z$bootstrap_draws)) 0L else active$calls
  result
}

# Record all R mode differences, including full IF. The authoritative reference
# is faster_mode=TRUE. This check does not tune or suppress any Python mismatch.
compare_modes <- function(a, b, path = "") {
  if (is.list(a) && !is.null(a$shape) && !is.null(a$values)) {
    same_support <- identical(a$shape, b$shape) && identical(a$nan, b$nan) &&
      identical(a$posinf, b$posinf) && identical(a$neginf, b$neginf)
    av <- unlist(a$values); bv <- unlist(b$values)
    max_abs <- if (same_support && length(av)) max(abs(av - bv)) else if (same_support) 0 else Inf
    return(list(list(field = path, same_support = same_support, max_abs = max_abs)))
  }
  if (is.list(a) && !is.null(names(a))) {
    stopifnot(identical(names(a), names(b)))
    out <- list()
    for (nm in setdiff(names(a), c("warnings", "bootstrap_calls"))) {
      out <- c(out, compare_modes(a[[nm]], b[[nm]], paste0(path, if (nzchar(path)) "." else "", nm)))
    }
    return(out)
  }
  stopifnot(identical(a, b))
  list()
}
references <- list()
mode_checks <- list()
for (z in cases) {
  cat("Generating", z$name, "\n")
  fast <- execute(z, TRUE)
  standard <- execute(z, FALSE)
  checks <- compare_modes(fast, standard)
  stopifnot(all(vapply(checks, function(x) x$same_support, logical(1))))
  references[[z$name]] <- fast
  mode_checks[[z$name]] <- checks
}
con <- gzfile(file.path(output, "references.json.gz"), "wt", compression = 9L)
jsonlite::write_json(references, con, digits = NA, auto_unbox = TRUE, null = "null")
close(con)
jsonlite::write_json(mode_checks, file.path(output, "r_mode_checks.json"), digits = NA,
                    auto_unbox = TRUE, pretty = TRUE, null = "null")
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
artifacts <- c(input_files, "cases.json", "references.json.gz", "r_mode_checks.json", "sessionInfo.txt")
hashes <- setNames(lapply(file.path(output, artifacts), sha256), artifacts)
hashes[["generate.R"]] <- sha256(script)
manifest <- list(schema_version = 1L, generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  sources = sources, packages = setNames(lapply(c("did", "DRDID", "fastglm", "Matrix", "data.table", "pkgload", "jsonlite", "digest"),
    function(p) as.character(packageVersion(p))), c("did", "DRDID", "fastglm", "Matrix", "data.table", "pkgload", "jsonlite", "digest")),
  R = R.version.string, platform = R.version$platform, external = as.list(extSoftVersion()),
  threads = as.list(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))),
  data_table_threads = data.table::getDTthreads(), case_count = length(cases),
  reference_faster_mode = TRUE, r_modes_evaluated = list(TRUE, FALSE),
  tolerance_policy = list(public = list(atol = 1e-7, rtol = 1e-6),
    wald = list(atol = 1e-5, rtol = 1e-5), influence_function = list(atol = 1e-7, rtol = 1e-6),
    missing_masks_and_support = "exact", declared_before_comparison = TRUE),
  bootstrap = list(draws = "Explicit frozen N(0,1) matrices; columns keyed by unit/cluster identity",
    hook = "did namespace run_multiplier_bootstrap(inf.func, biters, pl, cores)",
    reuse = "Same complete matrix on every bootstrap invocation; no RNG equivalence claim"),
  numeric_encoding = list(order = "row-major", masks = "zero-based flat positions", nonfinite_placeholder = 0,
    absent_field = "null", influence_rows = "unit_ids, sorted by original shared-input id"),
  artifacts_sha256 = hashes)
jsonlite::write_json(manifest, file.path(output, "manifest.json"), pretty = TRUE,
                    digits = NA, auto_unbox = TRUE, null = "null")
cat("Generated", length(cases), "references in", output, "\n")
