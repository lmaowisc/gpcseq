# Tests call the one-look engine before any design-specific sample-size gates.
aj_four <- function() {
  patients <- data.frame(id = c("T1", "T2", "C1", "C2"), arm = c(1L, 1L, 0L, 0L),
                         entry = 0, followup = c(2, 1.5, 1, 2),
                         death = c(FALSE, TRUE, TRUE, FALSE), state0 = "H")
  events <- data.frame(id = c("T2", "C1", "C2"), time = c(1.5, 1, 0.5),
                       from = "H", to = c("D", "D", "U"))
  list(p = patients, e = events)
}

test_that("AJ reproduces exact net time and pair-mean influences", {
  d <- aj_four()
  fit <- .aj_fit(d$p, d$e, 2, c("H", "U", "D"), "D", "D")
  # Direct times (treatment patients in rows) are (1, 1.5), (0.5, 0.5).
  expect_equal(fit$estimate, 0.875)
  expect_equal(fit$influence, c(0.375, -0.375, -0.125, 0.125))
  p <- d$p[c(4, 1, 3, 2), ]; p$arm <- 1L - p$arm
  reversed <- .aj_fit(p, d$e[3:1, ], 2, c("H", "U", "D"), "D", "D")
  expect_equal(reversed$estimate, -fit$estimate)
  expect_equal(reversed$influence, -fit$influence[c(4, 1, 3, 2)])
})

aj_mixed <- function() {
  p <- data.frame(id = c(paste0("C", 1:6), paste0("T", 1:6)),
                  arm = rep(0:1, each = 6), entry = 0,
                  state0 = c("H","H","U","H","U","H","H","U","H","H","H","U"),
                  followup = c(4,1.4,3.6,2.2,4,4,4,1.7,4,2.5,4,4),
                  death = seq_len(12) %in% c(3,10))
  e <- data.frame(id = c("C1","C1","C1","C2","C3","C3","C4","C4",
                         "T1","T1","T2","T3","T4","T4"),
                  time = c(0.8,1.5,2.6,0.5,1.0,3.6,0.8,1.8,0.8,1.5,1.2,1.7,0.8,2.5),
                  from = c("H","U","H","H","U","H","H","U",
                           "H","U","U","H","H","U"),
                  to = c("U","H","U","U","H","D","U","H",
                         "U","H","H","U","U","D"))
  list(p = p, e = e)
}

# Independent weighted point estimator: rebuild every state, risk set and
# transition matrix at every event. No calls to any AJ production helper.
aj_weighted_reference <- function(p, e, horizon, w) {
  states <- c("H", "U", "D"); ids <- match(e$id, p$id)
  grid <- sort(unique(c(0, e$time[e$time < horizon], horizon)))
  probs <- lapply(0:1, function(a) {
    in_arm <- p$arm == a
    v <- vapply(states, function(y) sum(w[in_arm & p$state0 == y]), 0.0) / sum(w[in_arm])
    out <- matrix(0, length(grid) - 1L, length(states))
    for (r in seq_len(nrow(out))) {
      t <- grid[r]; before <- as.character(p$state0)
      for (j in order(e$time)) if (e$time[j] < t) before[ids[j]] <- e$to[j]
      mat <- diag(length(states))
      for (f in seq_along(states)) {
        risk <- sum(w[in_arm & before == states[f] & p$followup >= t])
        if (!risk) next
        for (z in setdiff(seq_along(states), f)) {
          moving <- in_arm[ids] & e$time == t & e$from == states[f] & e$to == states[z]
          inc <- sum(w[ids[moving]]) / risk
          mat[f, z] <- inc; mat[f, f] <- mat[f, f] - inc
        }
      }
      v <- drop(v %*% mat); out[r, ] <- v
    }
    out
  })
  B <- outer(seq_along(states), seq_along(states), function(y, z) sign(z - y))
  sum(diff(grid) * rowSums((probs[[2L]] %*% B) * probs[[1L]]))
}

test_that("AJ patient influences equal independent case-weight refits with reversals, ties and censoring", {
  d <- aj_mixed(); p <- d$p
  for (tau in c(1.25, 2.5, 3.5)) {
    fit <- .aj_fit(p, d$e, tau, c("H", "U", "D"), "D", "D")
    ref <- function(w) aj_weighted_reference(p, d$e, tau, w)
    w <- rep(1, nrow(p)); eps <- 1e-5
    expect_equal(fit$estimate, ref(w), tolerance = 1e-12)
    derivative <- vapply(seq_len(nrow(p)), function(i) {
      up <- down <- w; up[i] <- 1 + eps; down[i] <- 1 - eps
      sum(p$arm == p$arm[i]) * (ref(up) - ref(down)) / (2 * eps)
    }, 0.0)
    expect_equal(fit$influence, derivative, tolerance = 2e-7)
    expect_lt(max(abs(tapply(fit$influence, p$arm, sum))), 1e-11)
  }
})

test_that("AJ enforces support on intervals and preserves absorbing states", {
  p <- data.frame(id = c("T","C"), arm = 1:0, entry = 0, followup = c(1,3),
                  death = FALSE, state0 = "H")
  fit <- function(...) .aj_fit(..., state_order = c("H", "U", "D"), absorbing = "D", death_state = "D")
  # A change at the horizon has zero duration; followup through tau suffices.
  expect_equal(fit(p, NULL, 1)$estimate, 0)
  expect_error(fit(p, NULL, 2), "occupied transient state H")
  p$death[1L] <- TRUE
  e <- data.frame(id = "T", time = 1, from = "H", to = "D")
  expect_equal(fit(p, e, 3)$estimate, -2)
  expect_equal(fit(p, e, 3)$influence, c(0, 0))
  expect_error(.aj_fit(p, e, 3, c("H", "U", "D"), "D"), "declared absorbing death_state")
  expect_error(fit(p, NULL, 3), "disagrees")
  p$death[1L] <- FALSE
  expect_error(fit(p, e, 3), "disagrees")
})

test_that("AJ detects loss of state support without requiring an observed transition", {
  p <- data.frame(id = c("T1","T2","C1"), arm = c(1L,1L,0L), entry = 0,
                  state0 = c("U","H","H"), followup = c(0,3,3), death = FALSE)
  expect_error(.aj_fit(p, NULL, 2, c("H", "U"), character()), "transient state U")
  fit <- .aj_fit(p, NULL, 2, c("H", "U"), "U")
  # U is declared clinically absorbing, independently of observed transitions.
  expect_equal(fit$estimate, -1)
  expect_equal(fit$influence, c(-1,1,0))
})

test_that("AJ rejects ambiguous times and invalid state paths", {
  d <- aj_four()
  fit <- function(e) .aj_fit(d$p, e, 2, c("H","U","D"), "D", "D")
  expect_error(fit(rbind(d$e, d$e[1L, ])), "multiple transitions")
  wrong <- d$e; wrong$from[1L] <- "U"
  expect_error(fit(wrong), "previous state")
  wrong <- d$e; wrong$time[1L] <- 1.51
  expect_error(fit(wrong), "after the end")
  wrong <- d$e; wrong$time[1L] <- 0
  expect_error(fit(wrong), "positive time")
})

test_that("an invalid AJ path is a data error even when final enrollment is insufficient", {
  d <- aj_four()
  plan <- gpc_design(looks = 2, spending = 0.05, method = "aj", horizons = 2,
                     state_order = c("H","U","D"), absorbing = "D", death_state = "D",
                     min_enrolled = 3L)
  # Each arm has only two patients. A malformed path must not silently close
  # the trial as an untested final look; the analyst can fix and retry the cut.
  wrong <- d$e; wrong$from[1L] <- "U"
  expect_error(gpc_update(plan, d$p, transitions = wrong), "previous state")
  expect_error(gpc_update(plan, d$p, transitions = d$e[-1L, ]), "disagrees")
  expect_false(plan$closed)
  expect_length(plan$log, 0L)
  retry <- gpc_update(plan, d$p, transitions = d$e)
  result <- gpc_results(retry)
  expect_identical(result$action, "final_untested")
  expect_equal(result$spent, 0)
  expect_match(result$reason, "Too few enrolled")
})

test_that("the cached AJ event table reproduces the direct fit", {
  d <- aj_mixed(); states <- c("H","U","D")
  indexed <- .aj_events(d$p, d$e, states, "D", "D")
  a <- .aj_fit(d$p, d$e, 2.5, states, "D", "D")
  b <- .aj_fit(d$p, d$e, 2.5, states, "D", "D", prepared = indexed)
  expect_identical(b, a)
})
