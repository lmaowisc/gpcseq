# The stored mock trial must reproduce the monitoring path documented in
# the README: three continuations, a month-72 upper rejection, and a
# refusal to update after stopping. Reference values come from an
# independent implementation of the shared-window statistic, its
# influence-based covariance, and the boundary recursion.

topcat_snapshot <- function(df, u) {
  p <- df[df$entry <= u, ]
  reach <- u - p$entry
  p$followup <- pmin(p$death_time, reach, na.rm = TRUE)
  p$death <- !is.na(p$death_time) & p$death_time <= reach
  seen <- !is.na(p$hosp_time) & p$hosp_time <= p$followup
  list(patients = p[, c("id", "arm", "entry", "followup", "death")],
       events = data.frame(id = p$id[seen], time = p$hosp_time[seen],
                           type = "hospitalization",
                           stringsAsFactors = FALSE))
}

test_that("topcat_mock is the documented final study database", {
  data(topcat_mock, envir = environment())
  expect_identical(dim(topcat_mock), c(500L, 5L))
  expect_identical(names(topcat_mock),
                   c("id", "arm", "entry", "death_time", "hosp_time"))
  expect_identical(sum(topcat_mock$arm), 250L)
  expect_true(all(topcat_mock$entry >= 0 & topcat_mock$entry <= 66))
  expect_identical(sum(!is.na(topcat_mock$death_time)), 205L)
  expect_identical(sum(!is.na(topcat_mock$hosp_time)), 218L)
  # events are recorded only through each patient's end-of-study reach
  reach <- 84 - topcat_mock$entry
  expect_true(all(topcat_mock$death_time <= reach, na.rm = TRUE))
  fu <- pmin(topcat_mock$death_time, reach, na.rm = TRUE)
  expect_true(all(topcat_mock$hosp_time <= fu, na.rm = TRUE))
})

test_that("topcat_mock monitoring reproduces the README path", {
  data(topcat_mock, envir = environment())
  trial <- gpc_design(looks = c(18, 36, 54, 72, 84), method = "shared")
  cut1 <- topcat_snapshot(topcat_mock, 18)
  trial <- gpc_update(trial, cut1$patients, cut1$events)
  cut2 <- topcat_snapshot(topcat_mock, 36)
  trial <- gpc_update(trial, cut2$patients, cut2$events)
  cut3 <- topcat_snapshot(topcat_mock, 54)
  trial <- gpc_update(trial, cut3$patients, cut3$events)
  cut4 <- topcat_snapshot(topcat_mock, 72)
  trial <- gpc_update(trial, cut4$patients, cut4$events)
  res <- gpc_results(trial)
  expect_identical(res$action,
                   c("continue", "continue", "continue", "reject_upper"))
  expect_identical(res$n0, c(66L, 129L, 198L, 250L))
  expect_identical(res$n1, c(66L, 124L, 196L, 250L))
  expect_equal(res$estimate,
               c(0.073003, 0.092461, 0.056895, 0.136064), tolerance = 1e-4)
  expect_equal(res$se,
               c(0.050234, 0.047205, 0.039852, 0.042688), tolerance = 1e-4)
  expect_equal(res$z,
               c(1.45326, 1.95868, 1.42766, 3.18743), tolerance = 1e-4)
  expect_equal(res$critical,
               c(4.23401, 2.99557, 2.48119, 2.23444), tolerance = 2e-3)
  # the trial is closed: the final look at month 84 must not be run
  cut5 <- topcat_snapshot(topcat_mock, 84)
  expect_error(gpc_update(trial, cut5$patients, cut5$events),
               "after stopping")
})
