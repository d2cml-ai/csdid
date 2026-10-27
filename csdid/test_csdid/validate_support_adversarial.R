# Independent estimator and influence-function equations. No did/DRDID package
# is loaded or called. jsonlite only serializes arrays and the case manifest.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
output <- normalizePath(args[[1L]], mustWork = TRUE)
cases <- jsonlite::fromJSON(file.path(output, "cases.json"), simplifyVector = FALSE)$targeted

panel_oracle <- function(data, method, intercept_only) {
  ids <- sort(unique(data$id))
  groups <- sort(unique(data$G[data$G > 0]))
  times <- sort(unique(data$period))
  result <- list(unit_ids = ids, group = numeric(), t = numeric(), att = numeric())
  influences <- list()
  oracle_warnings <- character()
  for (group in groups) {
    for (time_index in seq_along(times)[-1L]) {
      time <- times[[time_index]]
      pretime <- if (time < group) times[[time_index - 1L]] else max(times[times < group])
      pre <- data[data$period == pretime & data$G %in% c(0, group), ]
      post <- data[data$period == time & data$G %in% c(0, group), ]
      pre <- pre[order(pre$id), ]
      post <- post[order(post$id), ]
      stopifnot(identical(pre$id, post$id))
      early <- if (time > pretime) pre else post
      delta <- post$Y - pre$Y
      D <- as.numeric(pre$G == group)
      w <- early$weight / mean(early$weight)
      X <- if (intercept_only) matrix(1, nrow(pre), 1L) else cbind(1, early$X)
      m <- nrow(X)
      if (intercept_only) {
        # Exact weighted-DID ratio and its quotient-rule derivative: this does
        # not estimate a propensity score or outcome regression.
        treated_mean <- sum(w * D * delta) / sum(w * D)
        control_mean <- sum(w * (1 - D) * delta) / sum(w * (1 - D))
        att <- treated_mean - control_mean
        influence <- w * D * (delta - treated_mean) / mean(w * D) -
          w * (1 - D) * (delta - control_mean) / mean(w * (1 - D))
      } else {
        if (method %in% c("dr", "reg")) {
          wc <- w * (1 - D)
          gram <- crossprod(X, wc * X) / m
          beta <- solve(gram, colMeans(wc * delta * X))
          residual <- delta - drop(X %*% beta)
          beta_if <- (wc * residual * X) %*% solve(gram)
        } else residual <- delta
        if (method == "reg") {
          mass <- mean(w * D)
          att <- mean(w * D * residual) / mass
          influence <- (w * D * (residual - att) -
                         drop(beta_if %*% colMeans(w * D * X))) / mass
        } else {
          # A base-R weighted logit, independently refit at tighter convergence;
          # its observed score and information define the nuisance derivative.
          propensity <- withCallingHandlers(
            glm.fit(X, D, weights = w, family = binomial(),
                    control = glm.control(epsilon = 1e-12, maxit = 200L)),
            warning = function(condition) {
              oracle_warnings <<- c(oracle_warnings, conditionMessage(condition))
              invokeRestart("muffleWarning")
            })
          stopifnot(propensity$converged)
          ps <- pmin(propensity$fitted.values, 1 - 1e-6)
          trim <- rep(TRUE, length(D))
          trim[D == 0] <- ps[D == 0] < .995
          wt <- trim * w * D
          wc <- trim * w * ps * (1 - D) / (1 - ps)
          eta_t <- mean(wt * residual) / mean(wt)
          eta_c <- mean(wc * residual) / mean(wc)
          att <- eta_t - eta_c
          information <- crossprod(X, w * ps * (1 - ps) * X) / m
          gamma_if <- (w * (D - ps) * X) %*% solve(information)
          t_if <- wt * (residual - eta_t)
          c_if <- wc * (residual - eta_c) +
            drop(gamma_if %*% colMeans(wc * (residual - eta_c) * X))
          if (method == "dr") {
            t_if <- t_if - drop(beta_if %*% colMeans(wt * X))
            c_if <- c_if - drop(beta_if %*% colMeans(wc * X))
          }
          influence <- t_if / mean(wt) - c_if / mean(wc)
        }
      }
      full <- rep(0, length(ids))
      full[match(pre$id, ids)] <- influence * length(ids) / nrow(pre)
      result$group <- c(result$group, group)
      result$t <- c(result$t, time)
      result$att <- c(result$att, att)
      influences[[length(influences) + 1L]] <- full
    }
  }
  result$influence <- unname(do.call(cbind, influences))
  result$warnings <- unique(oracle_warnings)
  result
}

rc_oracle <- function(data) {
  # With no covariates every built-in method reduces to four weighted means;
  # the influence of each disjoint group-period sample follows directly from
  # differentiating its weighted ratio.
  ids <- sort(unique(data$id))
  groups <- sort(unique(data$G[data$G > 0]))
  times <- sort(unique(data$period))
  result <- list(unit_ids = ids, group = numeric(), t = numeric(), att = numeric())
  influences <- list()
  for (group in groups) {
    for (time_index in seq_along(times)[-1L]) {
      time <- times[[time_index]]
      pretime <- if (time < group) times[[time_index - 1L]] else max(times[times < group])
      att <- 0
      influence <- rep(0, length(ids))
      for (g in c(0, group)) {
        for (period in c(pretime, time)) {
          cell <- data[data$G == g & data$period == period, ]
          mean_y <- sum(cell$weight * cell$Y) / sum(cell$weight)
          sign <- if (xor(g == group, period == time)) -1 else 1
          att <- att + sign * mean_y
          influence[match(cell$id, ids)] <- length(ids) * sign * cell$weight *
            (cell$Y - mean_y) / sum(cell$weight)
        }
      }
      result$group <- c(result$group, group)
      result$t <- c(result$t, time)
      result$att <- c(result$att, att)
      influences[[length(influences) + 1L]] <- influence
    }
  }
  result$influence <- unname(do.call(cbind, influences))
  result
}

results <- list()
for (case in cases) {
  if (!isTRUE(case$oracle)) next
  data <- read.csv(file.path(output, case$data))
  result <- if (case$route == "rc") rc_oracle(data) else
    panel_oracle(data, case$method, case$formula == "~1")
  results[[case$name]] <- result
}
jsonlite::write_json(results, file.path(output, "oracles.json"), auto_unbox = TRUE,
                     digits = 17, na = "null")
writeLines(capture.output(sessionInfo()), file.path(output, "oracle-sessionInfo.txt"))
