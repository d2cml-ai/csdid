# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% [markdown]
# # Medicaid Expansion Results (R)
# Compute ATT and SE for unweighted and weighted configurations, then write CSV output.

# %% [markdown]
# ## Libraries

# %%
suppressPackageStartupMessages({
  library(did)
  library(dplyr)
  library(readr)
  library(stringr)
  library(broom)
})

# %% [markdown]
# ## Configuration and Paths

# %%
BITERS <- 25000
SEED <- 20240924
COVS <- c("perc_female","perc_white","perc_hispanic","unemp_rate","poverty_rate","median_income")

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}

if (basename(script_dir) == "scripts") {
  root_dir <- normalizePath(file.path(script_dir, ".."))
  out_dir <- script_dir
} else {
  root_dir <- normalizePath(script_dir)
  out_dir <- file.path(root_dir, "scripts")
}

data_path <- file.path(root_dir, "tests", "fixtures", "county_mortality_data.csv")
out_path <- file.path(out_dir, "medicaid_r_results.csv")

set.seed(SEED)

# %% [markdown]
# ## Data Preparation

# %%
data <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(state = str_sub(county, nchar(county) - 1, nchar(county))) %>%
  filter(!(state %in% c("DC", "DE", "MA", "NY", "VT"))) %>%
  filter(yaca == 2014 | is.na(yaca) | yaca > 2019) %>%
  mutate(
    perc_white = population_20_64_white / population_20_64 * 100,
    perc_hispanic = population_20_64_hispanic / population_20_64 * 100,
    perc_female = population_20_64_female / population_20_64 * 100,
    unemp_rate = unemp_rate * 100,
    median_income = median_income / 1000
  )

keep_cols <- c(
  "state","county","county_code","year","population_20_64","yaca",
  names(data)[startsWith(names(data), "perc_")], "crude_rate_20_64", COVS
)
keep_cols <- keep_cols[!duplicated(keep_cols)]
data <- data[, keep_cols]

cols_except_yaca <- setdiff(names(data), "yaca")
data <- data %>% filter(if_all(all_of(cols_except_yaca), ~ !is.na(.)))

data <- data %>%
  group_by(county_code) %>%
  filter(sum(year %in% c(2013, 2014)) == 2) %>%
  ungroup()

data <- data %>%
  group_by(county_code) %>%
  filter(sum(!is.na(crude_rate_20_64)) == 11) %>%
  ungroup()

# %% [markdown]
# ## Estimation Helper

# %%
short_data <- data %>%
  mutate(
    Treat = if_else(yaca == 2014 & !is.na(yaca), 1, 0),
    Post = if_else(year == 2014, 1, 0)
  ) %>%
  filter(year %in% c(2013, 2014)) %>%
  group_by(county_code) %>%
  mutate(set_wt = population_20_64[which(year == 2013)]) %>%
  ungroup()

data_cs <- short_data %>%
  mutate(treat_year = if_else(yaca == 2014 & !is.na(yaca), 2014, 0),
         county_code = as.numeric(county_code))

run_cs <- function(method, wt) {
  atts <- att_gt(
    yname = "crude_rate_20_64",
    tname = "year",
    idname = "county_code",
    gname = "treat_year",
    xformla = as.formula(paste("~", paste(COVS, collapse = "+"))),
    data = data_cs,
    panel = TRUE,
    control_group = "nevertreated",
    bstrap = TRUE,
    cband = TRUE,
    est_method = method,
    weightsname = wt,
    base_period = "universal",
    biters = BITERS
  )

  agg <- aggte(atts, type = "group", na.rm = TRUE, biters = BITERS)
  tidy_df <- broom::tidy(agg)
  row <- tidy_df[tidy_df$group == 2014, ]
  if (nrow(row) == 0) {
    stop("Missing group 2014 in aggte output")
  }

  data.frame(
    weighting = ifelse(is.null(wt), "unweighted", "weighted"),
    method = method,
    att = as.numeric(row$estimate[1]),
    se = as.numeric(row$std.error[1]),
    source = "r"
  )
}

# %% [markdown]
# ## Run and Save

# %%
out <- bind_rows(
  lapply(c("reg","ipw","dr"), function(m) run_cs(m, NULL)),
  lapply(c("reg","ipw","dr"), function(m) run_cs(m, "set_wt"))
)

write.csv(out, out_path, row.names = FALSE)
cat("R results:", out_path, "\n")
