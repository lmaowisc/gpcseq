# All fixtures use the public longitudinal-data constructor. Its input event
# times and the requested common horizon are measured from participant entry.
test_that("wide records preserve subject order, earliest admissions and exact horizon events", {
  p <- data.frame(id = c("late", "alive", "dead", "never", "early"),
                  arm = c(1, 0, 0, 1, 0), entry = 0,
                  followup = c(3, 2, 2, 3, 1.5),
                  death = c(TRUE, FALSE, TRUE, FALSE, TRUE))
  h <- data.frame(id = c("dead", "late", "dead", "never", "late", "alive"),
                  time = c(1.25, 2.25, 0.25, 2.1, 2, 0.5))
  d <- gpc_trial_data(p, events = list(admission = h))
  z <- gpc_transport_outcome(d, 2)
  expect_named(z, c("id", "death_time", "death_event", "hospitalization_time",
                    "hospitalization_event", "followup"))
  expect_identical(z$id, p$id)
  expect_identical(z$death_event, c(0L, 0L, 1L, 0L, 1L))
  expect_equal(z$death_time, c(2, 2, 2, 2, 1.5))
  expect_equal(z$followup, c(2, 2, 2, 2, 1.5))
  expect_equal(z$hospitalization_time, c(2, 0.5, 0.25, 2, NA_real_))
  expect_identical(z$hospitalization_event, c(1L, 1L, 1L, 0L, NA_integer_))
  # A non-event time of 2 and a first admission at 2 have distinct flags.
  expect_equal(z$hospitalization_time[1], z$hospitalization_time[4])
  expect_false(identical(z$hospitalization_event[1], z$hospitalization_event[4]))
})

test_that("complete non-events are representable and early loss cannot be recoded", {
  p <- data.frame(id = c("D", "S"), arm = 0:1, entry = 0,
                  followup = c(0.8, 2.5), death = c(TRUE, FALSE))
  empty <- gpc_trial_data(p)
  z <- gpc_transport_outcome(empty, 2)
  expect_identical(z$hospitalization_event, c(NA_integer_, 0L))
  expect_equal(z$hospitalization_time, c(NA_real_, 2))
  incomplete <- empty
  incomplete$patients$followup[2] <- 1.999
  expect_error(gpc_transport_outcome(incomplete, 2), "Incomplete follow-up.*S")
  bad <- empty
  bad$patients$death[2] <- NA
  expect_error(gpc_transport_outcome(bad, 2), "death")
  bad <- empty
  bad$events <- data.frame(id = "unknown", time = 0.5, type = "hospitalization")
  expect_error(gpc_transport_outcome(bad, 2), "unenrolled ID")
  bad <- empty
  bad$events <- NULL
  expect_error(gpc_transport_outcome(bad, 2), "prepared")
  expect_error(gpc_transport_outcome(empty, Inf), "positive finite")
  expect_error(gpc_transport_outcome(empty, 0), "positive finite")
  expect_error(gpc_transport_outcome(p, 2), "prepared")
})

test_that("native pairtransport compares complete bridged outcomes by first admission", {
  # Optional check: attach the submitted pairtransport package in the same R
  # session before testthat::test_local(). It is not a package dependency.
  if (!("package:pairtransport" %in% search()))
    skip("Attach pairtransport to enable its native comparison oracle.")
  native <- as.environment("package:pairtransport")
  p <- data.frame(id = c("A","B","C","Z","Z2"), arm = c(1,1,0,0,1), entry = 0,
                  followup = 2, death = c(FALSE,FALSE,FALSE,TRUE,TRUE))
  # B has more admissions than C, but B's first admission is later.
  h <- data.frame(id = c("B","B","C"), time = c(1.75,1,0.5))
  z <- gpc_transport_outcome(gpc_trial_data(p, events = list(admission = h)), 2)
  actual <- native$compare_pairs(z, z, native$hierarchical_survival(2))
  rank <- c(4,3,2,1,1)
  expected <- sign(outer(rank, rank, "-"))
  expect_equal(unname(actual$win - actual$loss), expected)
  expect_equal(actual$win[2,3], 1)
  expect_equal(actual$win[4,5] + actual$loss[4,5], 0)
})
