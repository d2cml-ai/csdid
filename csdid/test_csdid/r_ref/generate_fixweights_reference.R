## Generate R `did` reference values for fix_weights parity tests.
## Uses a shared time-varying-weight dataset (mpdta_tvw.csv) so R and Python
## consume identical inputs.
suppressMessages(library(did))

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."

d <- read.csv(file.path(outdir, "mpdta_tvw.csv"))

rows <- list()
for (fw in list(NULL, "base_period", "first_period", "varying")) {
  tag <- if (is.null(fw)) "none" else fw
  out <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                gname = "first.treat", xformla = ~1, data = d,
                control_group = "nevertreated", est_method = "reg",
                weightsname = "wt", fix_weights = fw,
                bstrap = FALSE, cband = FALSE)
  rows[[tag]] <- data.frame(fix_weights = tag, group = out$group,
                            t = out$t, att = out$att, se = out$se)
  cat("ran fix_weights =", tag, "\n")
}
ref <- do.call(rbind, rows)
write.csv(ref, file.path(outdir, "ref_fixweights.csv"), row.names = FALSE)
cat("did version:", as.character(packageVersion("did")), "\n")
cat("wrote ref_fixweights.csv (", nrow(ref), "rows)\n")
