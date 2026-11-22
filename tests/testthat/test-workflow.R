test_that("the public four-patient example follows the hand-scored trial", {
  one <- gpc_example(1)
  two <- gpc_example(2)
  fit <- gpc_design(looks = c(1, 2), spending = c(0.01, 0.05), method = "shared")
  fit <- gpc_update(fit, one$patients, one$events)
  fit <- gpc_update(fit, two$patients, two$events)
  ans <- gpc_results(fit)
  expect_equal(ans$action, c("deferred", "complete"))
  expect_equal(ans$spent, c(0, 0.05))
  expect_equal(ans$estimate[2], 0.5)
  expect_equal(ans$se[2], 0.5)
  expect_equal(ans$critical[2], qnorm(0.975), tolerance = 1e-6)
  expect_equal(unname(as.matrix(gpc_covariance(fit))), matrix(0.25, 1, 1))
  expect_equal(sort(gpc_influence(fit, 2)$phi), c(-0.5, -0.5, 0.5, 0.5))
})

test_that("natural and log scales agree with complete-data reductions", {
  d <- gpc_example(2)
  for (effect in c("net_benefit", "win_odds", "win_ratio")) {
    obj <- gpc_design(looks = 2, spending = 0.05, method = "ipcw",
                      effect = effect, horizons = 2)
    res <- gpc_results(gpc_update(obj, d$patients, d$events))
    if (effect == "net_benefit") {
      expect_equal(res$estimate, 0.5)
      expect_equal(res$se, 0.5)
    } else {
      expect_equal(res$estimate, log(3))
      expect_equal(res$natural, 3)
      expect_equal(res$se, 4 / 3)
    }
  }
  obj <- gpc_design(looks = 2, spending = 0.05, method = "aj",
                    effect = "net_win_time", horizons = 2,
                    state_order = c("well", "hospitalized", "dead"),
                    absorbing = "dead", death_state = "dead")
  aj <- gpc_update(obj, d$patients, transitions = d$transitions)
  expect_equal(gpc_results(aj)$estimate, 0.875)
})

test_that("matching ID yields the hand-derived covariance and hybrid boundary", {
  # All subjects live; each admission tier is a 0/1 indicator, and the two arm
  # distributions coincide at each look. Different people change status.
  ids <- c("C1", "C2", "C3", "T1", "T2", "T3")
  patients <- data.frame(id = ids, arm = rep(0:1, each = 3),
                         entry = 0, followup = 1, death = FALSE)
  events <- data.frame(id = c("C1", "T1", "C2", "T2"), time = c(.5,.5,1.5,1.5),
                       type = "hospitalization")
  obj <- gpc_design(looks = c(1, 2), spending = c(0.01, 0.05),
                    method = "shared", direction = "upper-final-two-sided",
                    hierarchy = gpc_hierarchy(admission_time = "none"))
  obj <- gpc_update(obj, patients, events[1:2, ])
  early <- gpc_results(obj)
  expect_equal(early$action, "continue")
  prior <- gpc_influence(obj, 1)
  patients$followup <- 2
  # Do not rely on row order for the covariance between updates.
  obj <- gpc_update(obj, patients[6:1, ], events[4:1, ])
  expect_equal(gpc_influence(obj, 1), prior)
  ans <- gpc_results(obj)
  expect_equal(ans$critical[1], early$critical)
  expect_equal(ans$action, c("continue", "complete"))
  expect_equal(ans$estimate, c(0, 0))
  expect_equal(unname(as.matrix(gpc_covariance(obj))),
               matrix(c(4, 2, 2, 4), ncol = 2) / 27)
  # Integrating the bivariate normal, with correlation 1/2 and upper-only
  # first rejection probability .01, gives this two-sided final boundary.
  expect_equal(ans$critical, c(2.32634787404, 2.03222310506), tolerance = 1e-5)
})
