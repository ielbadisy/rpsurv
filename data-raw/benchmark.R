# Benchmarks rpsurv against rstpm2::stpm2 fit time across sample size.
# Not run as part of R CMD check; re-run manually to refresh the vignette
# figures and cached results in data-raw/benchmark_results.rds.

suppressMessages({
  library(survival)
  library(rstpm2)
  devtools::load_all(".", quiet = TRUE)
})

simulate_data <- function(n, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  lp <- 0.5 * x1 - 0.3 * x2
  u <- runif(n)
  event_time <- (-log(u) / (0.01 * exp(lp)))^(1 / 1.2)
  cens_time <- rexp(n, 0.008)
  time <- pmin(event_time, cens_time)
  status <- as.numeric(event_time <= cens_time)
  data.frame(time = time, status = status, x1 = x1, x2 = x2)
}

ns <- c(1000, 5000, 20000, 50000, 100000, 200000, 500000)
reps <- 3

results <- do.call(rbind, lapply(ns, function(n) {
  d <- simulate_data(n)

  t_rpsurv <- replicate(reps, {
    st <- system.time(rpsurv(Surv(time, status) ~ x1 + x2, data = d, df = 4))
    st[["elapsed"]]
  })

  t_rstpm2 <- if (n <= 100000) {
    replicate(reps, {
      st <- system.time(stpm2(Surv(time, status) ~ x1 + x2, data = d, df = 4))
      st[["elapsed"]]
    })
  } else {
    NA_real_
  }

  data.frame(
    n = n,
    rpsurv = median(t_rpsurv),
    rstpm2 = if (all(is.na(t_rstpm2))) NA_real_ else median(t_rstpm2)
  )
}))

print(results)
saveRDS(results, "data-raw/benchmark_results.rds")

png("vignettes/figures/benchmark_speed.png", width = 1400, height = 950, res = 150)
op <- par(mar = c(4.5, 4.5, 2, 1))
plot(results$n, results$rpsurv, type = "b", pch = 19, col = "#1b6ca8", log = "xy",
     ylim = range(c(results$rpsurv, results$rstpm2), na.rm = TRUE),
     xlab = "Number of observations (n)", ylab = "Fit time in seconds (log scale)",
     main = "rpsurv vs rstpm2::stpm2 fit time")
lines(results$n, results$rstpm2, type = "b", pch = 17, col = "#c0392b")
legend("topleft", legend = c("rpsurv (Rcpp + RcppParallel)", "rstpm2::stpm2"),
       col = c("#1b6ca8", "#c0392b"), pch = c(19, 17), lty = 1, bty = "n")
par(op)
dev.off()

cat("speedup at largest common n:\n")
print(within(results, speedup <- rstpm2 / rpsurv))
