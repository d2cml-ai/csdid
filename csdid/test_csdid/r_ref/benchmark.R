## Benchmark R `did` faster_mode vs standard. Writes bench_r.csv.
## Run: Rscript csdid/test_csdid/r_ref/benchmark.R
suppressMessages(library(did))

here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) here <- args[1]
sim <- file.path(here, "sim")

datasets <- list(
  c("tp8_dyn",    file.path(sim, "tp8_dyn.csv")),
  c("tp10_const", file.path(sim, "tp10_const.csv")),
  c("big_tp10",   file.path(sim, "big_tp10.csv"))
)

tm <- function(d, fm, reps = 5) {
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- Sys.time()
    att_gt(yname = "Y", tname = "period", idname = "id", gname = "G", data = d,
           control_group = "nevertreated", xformla = ~X, est_method = "dr",
           bstrap = FALSE, faster_mode = fm)
    ts[i] <- as.numeric(Sys.time() - t0, units = "secs")
  }
  min(ts) * 1000
}

rows <- list()
cat(sprintf("%-14s %7s %13s %10s %8s\n", "dataset", "rows", "R(Standard)", "R(Fast)", "speedup"))
for (ds in datasets) {
  d <- read.csv(ds[2])
  std <- tm(d, FALSE); fast <- tm(d, TRUE)
  rows[[ds[1]]] <- data.frame(dataset = ds[1], rows = nrow(d),
                              r_standard_ms = round(std, 1), r_fast_ms = round(fast, 1),
                              r_speedup = round(std / fast, 2))
  cat(sprintf("%-14s %7d %13.1f %10.1f %7.2fx\n", ds[1], nrow(d), std, fast, std / fast))
}
write.csv(do.call(rbind, rows), file.path(here, "bench_r.csv"), row.names = FALSE)
cat("\nwrote bench_r.csv\n")
