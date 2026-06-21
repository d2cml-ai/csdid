## Generate R `did` references for previously-unverified scenarios
## (repeated cross sections, universal base period, anticipation, weighted,
## clustered analytical SE). Shared data: mpdta_extra.csv (mpdta + wt + clust).
suppressMessages(library(did))
args <- commandArgs(trailingOnly = TRUE)
simdir <- if (length(args) >= 1) args[1] else "."
d <- read.csv(file.path(simdir, "mpdta_extra.csv"))
common <- list(yname = "lemp", tname = "year", idname = "countyreal",
               gname = "first.treat", data = d, bstrap = FALSE, cband = FALSE)
run <- function(scn, ...) {
  o <- do.call(att_gt, c(common, list(...)))
  data.frame(scenario = scn, group = o$group, t = o$t, att = o$att, se = o$se)
}
rows <- rbind(
  run("rc",            panel = FALSE, control_group = "nevertreated", est_method = "reg"),
  run("universal",     base_period = "universal", control_group = "nevertreated", est_method = "reg"),
  run("anticipation1", anticipation = 1, control_group = "nevertreated", est_method = "reg"),
  run("weighted",      weightsname = "wt", control_group = "nevertreated", est_method = "reg"),
  run("clustered",     clustervars = "clust", control_group = "nevertreated", est_method = "reg")
)
write.csv(rows, file.path(simdir, "ref_gaps.csv"), row.names = FALSE)
cat("wrote ref_gaps.csv:", nrow(rows), "rows; did", as.character(packageVersion("did")), "\n")
