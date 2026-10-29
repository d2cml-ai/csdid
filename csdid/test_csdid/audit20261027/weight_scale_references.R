# Export complete output for the independent weight-scale audit. This script
# never sources generate.R and never writes to the committed fixture directory.
args <- commandArgs(trailingOnly = TRUE)
flag <- function(name) {
  at <- match(paste0('--', name), args)
  if (is.na(at) || at == length(args)) stop('Missing --', name)
  args[at + 1L]
}
output <- normalizePath(flag('output'), mustWork = TRUE)
did_path <- normalizePath(flag('did'), mustWork = TRUE)
drdid_path <- normalizePath(flag('drdid'), mustWork = TRUE)
suppressPackageStartupMessages(pkgload::load_all(drdid_path, quiet = TRUE))
suppressPackageStartupMessages(pkgload::load_all(did_path, quiet = TRUE))
stopifnot(normalizePath(getNamespaceInfo('did', 'path')) == did_path,
          normalizePath(getNamespaceInfo('DRDID', 'path')) == drdid_path)
data.table::setDTthreads(1L)

# The R-side portable DGP is also used in did's regression test added with its
# weight-scale correction. Seed, sample and all outcomes are fixed in advance.
if ('--generate-inputs' %in% args) {
  RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
  set.seed(20261027L)
  n <- 160L
  units <- data.frame(id = seq_len(n), G = rep(c(0, 2, 3, 4), each = n / 4),
    X = rnorm(n), a = rnorm(n), weight = .5 + (seq_len(n) %% 7) / 10)
  d <- merge(units, data.frame(period = 1:4))
  d$Y <- d$a + .2 * d$period + .1 * d$X * d$period +
    (d$G > 0 & d$period >= d$G) * (.5 + .2 * d$G) * (d$period - d$G + 1) +
    .1 * rnorm(nrow(d))
  write.csv(d, file.path(output, 'multi_panel.csv'), row.names = FALSE)
  d$id <- seq_len(nrow(d))
  write.csv(d, file.path(output, 'multi_rc.csv'), row.names = FALSE)
  quit(status = 0L)
}

# Numeric schema and retained-ID mapping follow r_ref/current/generate.R.
# Shapes and explicit nonfinite masks prevent coercion from hiding a mismatch.
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
execute <- function(z, fast) {
  data <- read.csv(file.path(output, z$data))
  kw <- z$expected_r_kwargs
  kw$xformla <- as.formula(kw$xformla)
  kw$data <- data
  kw$faster_mode <- fast
  recorded_warnings <- character()
  result <- withCallingHandlers({
    fit <- do.call(did::att_gt, kw)
    ids <- suppressWarnings(as.numeric(rownames(fit$inffunc)))
    stopifnot(length(ids) == nrow(fit$inffunc), all(is.finite(ids)), !anyDuplicated(ids))
    if (!z$constructor$panel) {
      retained <- fit$DIDparams$data
      internal_ids <- retained[[fit$DIDparams$idname]]
      original_ids <- retained[[z$constructor$idname]]
      stopifnot(!anyDuplicated(internal_ids), !anyDuplicated(original_ids))
      position <- match(ids, internal_ids)
      stopifnot(!anyNA(position))
      ids <- as.numeric(original_ids[position])
    }
    ord <- order(ids)
    aggs <- setNames(lapply(c('simple', 'group', 'dynamic', 'calendar'), function(kind) {
      a <- do.call(did::aggte, c(list(MP = fit, type = kind), z$expected_r_aggregation_kwargs))
      influence <- lapply(a$inf.function, function(value) {
        matrix <- as.matrix(value)
        stopifnot(nrow(matrix) == length(ids))
        numeric_payload(matrix[ord, , drop = FALSE], influence = TRUE)
      })
      names(influence) <- gsub('.', '_', names(influence), fixed = TRUE)
      list(overall_att = numeric_payload(a$overall.att), overall_se = numeric_payload(a$overall.se),
        egt = numeric_payload(a$egt), att_egt = numeric_payload(a$att.egt),
        se_egt = numeric_payload(a$se.egt), crit_val_egt = numeric_payload(a$crit.val.egt),
        inf_function = influence)
    }), c('simple', 'group', 'dynamic', 'calendar'))
    list(unit_ids = as.list(ids[ord]), fit = list(group = numeric_payload(fit$group),
      t = numeric_payload(fit$t), att = numeric_payload(fit$att), se = numeric_payload(fit$se),
      c = numeric_payload(fit$c), l_se = numeric_payload(fit$att - fit$c * fit$se),
      u_se = numeric_payload(fit$att + fit$c * fit$se), W = numeric_payload(fit$W),
      Wpval = numeric_payload(fit$Wpval),
      inffunc = numeric_payload(as.matrix(fit$inffunc)[ord, , drop = FALSE], influence = TRUE),
      cband = isTRUE(fit$DIDparams$cband)), aggs = aggs)
  }, warning = function(warning) {
    recorded_warnings <<- c(recorded_warnings, conditionMessage(warning))
    invokeRestart('muffleWarning')
  })
  result$warnings <- as.list(recorded_warnings)
  result
}
cases <- jsonlite::fromJSON(file.path(output, 'cases.json'), simplifyVector = FALSE)
references <- list()
for (z in cases) for (fast in c(FALSE, TRUE)) {
  name <- paste0(z$name, if (fast) '__fast' else '__standard')
  cat('Exporting ', name, '\n', sep = '')
  references[[name]] <- execute(z, fast)
}
connection <- gzfile(file.path(output, 'r_references.json.gz'), 'wt', compression = 9L)
jsonlite::write_json(references, connection, digits = NA, auto_unbox = TRUE, null = 'null')
close(connection)
writeLines(capture.output(sessionInfo()), file.path(output, 'sessionInfo.txt'))
