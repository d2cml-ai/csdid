## Generate R `did` reference for a FACTOR covariate (parity test).
## Shared dataset factor_cov.csv so R and Python consume identical inputs.
suppressMessages(library(did))
args <- commandArgs(trailingOnly = TRUE)
simdir <- if (length(args) >= 1) args[1] else "."
d <- read.csv(file.path(simdir, "factor_cov.csv"))
d$cat <- factor(d$cat)
out <- att_gt(yname = "Y", tname = "period", idname = "id", gname = "G",
              xformla = ~cat, data = d, control_group = "nevertreated",
              est_method = "reg", bstrap = FALSE)
write.csv(data.frame(group = out$group, t = out$t, att = out$att, se = out$se),
          file.path(simdir, "ref_factor.csv"), row.names = FALSE)
cat("wrote ref_factor.csv (", length(out$att), "cells); did", as.character(packageVersion("did")), "\n")
