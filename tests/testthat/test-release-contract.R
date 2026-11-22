test_that("trivariate corner budget reflects the backend tolerance floor", {
  R <- matrix(0.5, 3, 3); diag(R) <- 1
  # Eight inclusion-exclusion corners: a 2e-13 request implies per-corner
  # epsilon below the native 1e-14 floor, so the reported budget must carry
  # the floor, not the request. The probability itself is unaffected.
  tight <- .gpc_dependent_rect(rep(-2, 3), rep(2, 3), R, 2e-13)
  expect_gte(unname(tight["e"]), 8 * 1e-14)
  expect_gt(unname(tight["p"]), 0.88)
  expect_lt(unname(tight["p"]), 0.89)
  # An attainable request keeps a budget within its tolerance.
  loose <- .gpc_dependent_rect(rep(-2, 3), rep(2, 3), R, 1e-7)
  expect_lte(unname(loose["e"]), 1e-7)
  expect_equal(unname(tight["p"]), unname(loose["p"]), tolerance = 1e-9)
})

test_that("public argument errors name the requirement", {
  expect_error(gpc_design(2, 0.05, "ipcw", horizons = -1),
               "horizons must be positive finite follow-up times")
  expect_error(gpc_design(c(1, 2), c(0.01, 0.05), "ipcw", horizons = c(1, 2, 3)),
               "one common value, or one per look")
  expect_error(gpc_design(2, 0.05, "aj", horizons = 2,
                          state_order = c("well", "dead"), absorbing = "gone",
                          death_state = "dead"),
               "absorbing must be distinct state names")
  plan <- gpc_design(c(1, 2), c(0.01, 0.05))
  expect_error(gpc_fit(plan, 0),
               "look must be one scheduled look index between 1 and 2")
  expect_error(gpc_fit(plan, 3),
               "look must be one scheduled look index between 1 and 2")
})
