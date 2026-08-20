# Benchmarks rpsurv against the two canonical Royston-Parmar/flexible
# parametric implementations, rstpm2::stpm2 and flexsurv::flexsurvspline,
# across sample size. Not run as part of R CMD check; re-run manually to
# refresh the vignette figures and cached results in
# data-raw/benchmark_results.rds.

suppressMessages({
  library(survival)
  library(rstpm2)
  library(flexsurv)
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
# flexsurvspline() is much slower per-fit; cap it to keep the benchmark tractable
flexsurv_max_n <- 20000
rstpm2_max_n <- 200000
reps <- 3

time_median <- function(expr_fun, reps) {
  ts <- replicate(reps, system.time(expr_fun())[["elapsed"]])
  median(ts)
}

results <- do.call(rbind, lapply(ns, function(n) {
  d <- simulate_data(n)

  t_rpsurv <- time_median(function() rpsurv(Surv(time, status) ~ x1 + x2, data = d, df = 4), reps)

  t_rstpm2 <- if (n <= rstpm2_max_n) {
    time_median(function() stpm2(Surv(time, status) ~ x1 + x2, data = d, df = 4), reps)
  } else {
    NA_real_
  }

  t_flexsurv <- if (n <= flexsurv_max_n) {
    time_median(function() flexsurvspline(Surv(time, status) ~ x1 + x2, data = d, k = 3), reps)
  } else {
    NA_real_
  }

  data.frame(n = n, rpsurv = t_rpsurv, rstpm2 = t_rstpm2, flexsurv = t_flexsurv)
}))

print(results)
saveRDS(results, "data-raw/benchmark_results.rds")
file.copy("data-raw/benchmark_results.rds", "inst/extdata/benchmark_results.rds", overwrite = TRUE)

png("vignettes/figures/benchmark_speed.png", width = 1400, height = 950, res = 150)
op <- par(mar = c(4.5, 4.5, 2, 1))
plot(results$n, results$rpsurv, type = "b", pch = 19, col = "#1b6ca8", log = "xy",
     ylim = range(c(results$rpsurv, results$rstpm2, results$flexsurv), na.rm = TRUE),
     xlab = "Number of observations (n)", ylab = "Fit time in seconds (log scale)",
     main = "rpsurv vs rstpm2 vs flexsurv fit time")
lines(results$n, results$rstpm2, type = "b", pch = 17, col = "#c0392b")
lines(results$n, results$flexsurv, type = "b", pch = 15, col = "#27ae60")
legend("topleft", legend = c("rpsurv (Rcpp + RcppParallel)", "rstpm2::stpm2", "flexsurv::flexsurvspline"),
       col = c("#1b6ca8", "#c0392b", "#27ae60"), pch = c(19, 17, 15), lty = 1, bty = "n")
par(op)
dev.off()

cat("speedup vs rstpm2 and flexsurv:\n")
print(within(results, {
  speedup_rstpm2 <- rstpm2 / rpsurv
  speedup_flexsurv <- flexsurv / rpsurv
}))
