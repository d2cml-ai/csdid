## Generate R `did` reference values for csdid Python parity tests.
## Outputs tidy CSVs into the same directory.
suppressMessages(library(did))

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."
datadir <- if (length(args) >= 2) args[2] else "."

mpdta <- read.csv(file.path(datadir, "mpdta.csv"))
simdt <- read.csv(file.path(datadir, "sim_data.csv"))

attgt_rows <- list()
aggte_rows <- list()

run_scenario <- function(scn, data, yname, tname, idname, gname,
                         xformla = ~1, control_group = "nevertreated",
                         est_method = "dr", base_period = "varying",
                         panel = TRUE) {
  out <- att_gt(yname = yname, tname = tname, idname = idname, gname = gname,
                xformla = xformla, data = data, control_group = control_group,
                est_method = est_method, base_period = base_period,
                bstrap = FALSE, cband = FALSE, panel = panel)
  attgt_rows[[scn]] <<- data.frame(
    scenario = scn, group = out$group, t = out$t,
    att = out$att, se = out$se
  )
  for (tp in c("simple", "group", "dynamic", "calendar")) {
    ag <- aggte(out, type = tp, bstrap = FALSE)
    egt <- if (is.null(ag$egt)) NA else ag$egt
    att_egt <- if (is.null(ag$att.egt)) NA else ag$att.egt
    se_egt <- if (is.null(ag$se.egt)) NA else ag$se.egt
    aggte_rows[[paste(scn, tp)]] <<- data.frame(
      scenario = scn, type = tp,
      egt = egt, att_egt = att_egt, se_egt = se_egt,
      overall_att = ag$overall.att, overall_se = ag$overall.se
    )
  }
  cat("ran:", scn, "\n")
}

run_scenario("mpdta_nev_dr", mpdta, "lemp", "year", "countyreal", "first.treat",
             xformla = ~1, control_group = "nevertreated", est_method = "dr")
run_scenario("mpdta_nyt_dr", mpdta, "lemp", "year", "countyreal", "first.treat",
             xformla = ~1, control_group = "notyettreated", est_method = "dr")
run_scenario("mpdta_nev_reg_cov", mpdta, "lemp", "year", "countyreal", "first.treat",
             xformla = ~lpop, control_group = "nevertreated", est_method = "reg")
run_scenario("mpdta_nev_ipw", mpdta, "lemp", "year", "countyreal", "first.treat",
             xformla = ~1, control_group = "nevertreated", est_method = "ipw")
run_scenario("sim_nev_dr", simdt, "Y", "period", "id", "G",
             xformla = ~X, control_group = "nevertreated", est_method = "dr")

attgt <- do.call(rbind, attgt_rows)
aggte_df <- do.call(rbind, aggte_rows)
write.csv(attgt, file.path(outdir, "ref_attgt.csv"), row.names = FALSE)
write.csv(aggte_df, file.path(outdir, "ref_aggte.csv"), row.names = FALSE)
cat("did version:", as.character(packageVersion("did")), "\n")
cat("wrote ref_attgt.csv (", nrow(attgt), "rows) and ref_aggte.csv (", nrow(aggte_df), "rows)\n")
