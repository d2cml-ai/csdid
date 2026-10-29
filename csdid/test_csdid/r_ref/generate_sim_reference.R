suppressMessages(library(did))
simdir <- "C:/Users/h_min/repos/csdid/test_csdid/r_ref/sim"
names <- c("tp2_const","tp4_const","tp4_dyn","tp5_dyn","tp8_dyn","tp10_const")
rows <- list()
for (nm in names) {
  d <- read.csv(file.path(simdir, paste0(nm, ".csv")))
  for (cg in c("nevertreated","notyettreated")) {
    for (em in c("dr","reg")) {
      out <- tryCatch(att_gt(yname="Y", tname="period", idname="id", gname="G",
                    xformla=~X, data=d, control_group=cg, est_method=em,
                    bstrap=FALSE, cband=FALSE), error=function(e) NULL)
      if (is.null(out)) next
      rows[[paste(nm,cg,em)]] <- data.frame(dataset=nm, control=cg, est=em,
            group=out$group, t=out$t, att=out$att, se=out$se)
    }
  }
  cat("ran", nm, "\n")
}
ref <- do.call(rbind, rows)
write.csv(ref, file.path(simdir, "ref_sim.csv"), row.names=FALSE)
cat("wrote ref_sim.csv (", nrow(ref), "rows )\n")