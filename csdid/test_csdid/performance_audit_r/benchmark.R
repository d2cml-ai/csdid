# Current-source R comparison for benchmark_v043.py. See README.md for replay.
args <- commandArgs(trailingOnly = TRUE)
flag <- function(name, default = NULL) {
  at <- match(paste0("--", name), args)
  if (is.na(at)) return(default)
  if (at == length(args)) stop("Missing value for --", name)
  args[at + 1L]
}
did_path <- normalizePath(flag("did"), mustWork = TRUE)
drdid_path <- normalizePath(flag("drdid"), mustWork = TRUE)
input <- normalizePath(flag("input"), mustWork = TRUE)
output <- flag("output")
if (is.null(output)) stop("--output is required")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output)
reps <- as.integer(flag("reps", "5"))
stopifnot(length(reps) == 1L, is.finite(reps), reps >= 1L)
profiles <- "--profiles" %in% args
case_names <- strsplit(flag("cases", "analytic,clustered,bootstrap"), ",", fixed = TRUE)[[1]]
stopifnot(all(case_names %in% c("analytic", "clustered", "bootstrap")))

suppressPackageStartupMessages(pkgload::load_all(drdid_path, quiet = TRUE))
suppressPackageStartupMessages(pkgload::load_all(did_path, quiet = TRUE))
data.table::setDTthreads(1L)

source_info <- function(path) {
  files <- c(file.path(path, "DESCRIPTION"), list.files(file.path(path, "R"),
             pattern = "[.]R$", full.names = TRUE),
             list.files(file.path(path, "src"), full.names = TRUE))
  files <- files[file.exists(files) & !dir.exists(files)]
  md5 <- tools::md5sum(files)
  names(md5) <- substring(names(md5), nchar(path) + 2L)
  list(path = path,
       head = system2("git", c("-C", shQuote(path), "rev-parse", "HEAD"), stdout = TRUE),
       status = system2("git", c("-C", shQuote(path), "status", "--short"), stdout = TRUE),
       md5 = as.list(md5))
}
env <- list(R = R.version.string, platform = R.version$platform,
            packages = lapply(c("did", "DRDID", "data.table", "fastglm", "Matrix",
                                "BMisc", "pkgload"), function(p) {
              list(name = p, version = as.character(packageVersion(p)),
                   path = find.package(p))
            }),
            external = as.list(extSoftVersion()),
            threads = as.list(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"))),
            data_table_threads = data.table::getDTthreads(),
            sources = list(did = source_info(did_path), DRDID = source_info(drdid_path)),
            scope = "att_gt preprocessing + fitting + simple/group/dynamic/calendar aggte",
            bootstrap = "399 independent R multipliers; R and Python RNG draws are NOT shared",
            input_json = jsonlite::fromJSON(file.path(input, "inputs.json"), simplifyVector = FALSE))
jsonlite::write_json(env, file.path(output, "environment.json"), pretty = TRUE,
                     auto_unbox = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))

inputs <- env$input_json
datasets <- lapply(case_names, function(name) {
  filename <- file.path(input, inputs[[name]]$filename)
  stopifnot(unname(tools::md5sum(filename)) == inputs[[name]]$md5)
  read.csv(filename)
})
names(datasets) <- case_names

execute <- function(data, fast, boot, clustered) {
  set.seed(1234)
  warnings <- character()
  withCallingHandlers({
    start <- proc.time()[["elapsed"]]
    fit <- did::att_gt(yname = "Y", tname = "period", idname = "id", gname = "G",
        data = data, control_group = "nevertreated", xformla = ~X,
        faster_mode = fast, weightsname = "wt", est_method = "dr",
        clustervars = if (clustered) "cluster" else NULL,
        bstrap = boot, cband = boot, biters = 399L,
        base_period = "varying", pl = FALSE, cores = 1L)
    fitted <- proc.time()[["elapsed"]]
    aggs <- list()
    for (kind in c("simple", "group", "dynamic", "calendar")) {
      aggs[[kind]] <- did::aggte(fit, type = kind, bstrap = boot, cband = boot)
    }
    end <- proc.time()[["elapsed"]]
    list(fit = fit, aggs = aggs, warnings = warnings,
         timing = c(att_gt_s = fitted - start, aggregation_s = end - fitted,
                    total_s = end - start))
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}

# Normalize full influence functions using authoritative MP row IDs. Fast mode
# orders cohorts differently from standard mode, so raw row order is not valid.
snapshot <- function(run) {
  r <- run$fit
  ids <- suppressWarnings(as.numeric(rownames(r$inffunc)))
  stopifnot(length(ids) == nrow(r$inffunc), all(is.finite(ids)), !anyDuplicated(ids))
  ord <- order(ids)
  aggs <- lapply(run$aggs, function(a) {
    inf <- lapply(a$inf.function, function(v) {
      m <- as.matrix(v)
      stopifnot(nrow(m) == length(ids))
      m[ord, , drop = FALSE]
    })
    names(inf) <- gsub(".", "_", names(inf), fixed = TRUE)
    list(overall_att = a$overall.att, overall_se = a$overall.se, egt = a$egt,
         att_egt = a$att.egt, se_egt = a$se.egt, crit_val_egt = a$crit.val.egt,
         inf_function = inf)
  })
  list(fit = list(group = r$group, t = r$t, att = r$att, se = r$se,
       c = r$c, l_se = r$att - r$c * r$se, u_se = r$att + r$c * r$se,
       W = r$W, Wpval = r$Wpval, inffunc = as.matrix(r$inffunc)[ord, , drop = FALSE]),
       aggs = aggs, unit_ids = ids[ord], warnings = run$warnings)
}

raw <- list()
references <- list()
for (name in case_names) {
  data <- datasets[[name]]
  boot <- isTRUE(inputs[[name]]$bootstrap)
  clustered <- isTRUE(inputs[[name]]$clustered)
  for (fast in c(FALSE, TRUE)) invisible(execute(data, fast, boot, clustered))
  for (rep in seq_len(reps)) {
    modes <- if (rep %% 2L) c(FALSE, TRUE) else c(TRUE, FALSE)
    for (fast in modes) {
      invisible(gc())
      run <- execute(data, fast, boot, clustered)
      key <- paste0(name, "/fast=", if (fast) "True" else "False")
      current <- snapshot(run)
      if (is.null(references[[key]])) {
        references[[key]] <- current
      } else if (!isTRUE(all.equal(references[[key]], current, tolerance = 0))) {
        stop("Repeat changed full result: ", key)
      }
      raw[[length(raw) + 1L]] <- data.frame(case = name, fast = fast, rep = rep,
            att_gt_s = run$timing[["att_gt_s"]],
            aggregation_s = run$timing[["aggregation_s"]],
            total_s = run$timing[["total_s"]])
      cat(sprintf("%s fast=%s rep=%d total=%.3f s\n", name, fast, rep,
                  run$timing[["total_s"]]))
    }
  }
}
raw <- do.call(rbind, raw)
write.csv(raw, file.path(output, "timings.csv"), row.names = FALSE)
summary <- aggregate(cbind(att_gt_s, aggregation_s, total_s) ~ case + fast,
                     data = raw, FUN = median)
write.csv(summary, file.path(output, "medians.csv"), row.names = FALSE)
jsonlite::write_json(list(raw = raw, medians = summary), file.path(output, "results.json"),
                    auto_unbox = TRUE, pretty = TRUE, digits = NA)
con <- gzfile(file.path(output, "snapshot.json.gz"), "wt")
jsonlite::write_json(references, con, digits = NA, na = "string", null = "null",
                    auto_unbox = FALSE, matrix = "rowmajor")
close(con)

# Profiling is a separate execution after every measured run; it cannot inflate
# the medians above. Memory columns describe allocations, not process peak RSS.
if (profiles) for (name in case_names) for (fast in c(FALSE, TRUE)) {
  stem <- file.path(output, paste0(name, "-fast", as.integer(fast)))
  invisible(gc())
  Rprof(paste0(stem, ".Rprof"), interval = .001, memory.profiling = TRUE)
  tryCatch(invisible(execute(datasets[[name]], fast,
           isTRUE(inputs[[name]]$bootstrap), isTRUE(inputs[[name]]$clustered))),
           finally = Rprof(NULL))
  prof <- summaryRprof(paste0(stem, ".Rprof"), memory = "both")
  write.csv(prof$by.total, paste0(stem, "-total.csv"))
  write.csv(prof$by.self, paste0(stem, "-self.csv"))
}
print(summary)
