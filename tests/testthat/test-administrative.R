admin_subjects <- function() {
  data.frame(id = c("T1", "T2", "C1", "C2"), arm = c(1, 1, 0, 0),
             entry = 0, followup = 1 + .Machine$double.eps,
             death = FALSE, state0 = "H")
}

test_that("administrative preparation is explicit, audited, and leaves other rows intact", {
  p <- admin_subjects()
  p$followup[2:3] <- c(.5, 1); p$death[3] <- TRUE
  legacy <- gpc_trial_data(p)
  expect_identical(gpc_trial_data(p, administrative = NULL, look = NULL), legacy)
  flags <- c(TRUE, FALSE, FALSE, TRUE)
  cut <- gpc_trial_data(p, administrative = flags, look = 1)
  expect_identical(cut$patients$followup, c(1, .5, 1, 1))
  expect_identical(cut$patients[2:3, ], legacy$patients[2:3, ])
  expect_identical(cut$administrative$audit$id, c("T1", "C2"))
  expect_identical(cut$administrative$audit$original, p$followup[c(1, 4)])
  expect_identical(cut$administrative$audit$delta, rep(-.Machine$double.eps, 2))
  expect_identical(unserialize(serialize(cut, NULL)), cut)
  expect_identical(gpc_trial_data(p, administrative = setNames(flags, p$id), look = 1), cut)
  expect_identical(gpc_trial_data(p, administrative = rep(FALSE, 4), look = 1)$patients,
                   legacy$patients)
  expect_error(gpc_trial_data(p, look = 1), "together")
  expect_error(gpc_trial_data(p, administrative = flags), "together")
  for (bad in list(TRUE, 1:4, c(TRUE, NA, FALSE, TRUE), matrix(flags, 2),
                   setNames(flags, rev(p$id))))
    expect_error(gpc_trial_data(p, administrative = bad, look = 1), "administrative")
  for (bad in list(NA_real_, Inf, -1, numeric(), c(1, 2), as.Date("2020-01-01")))
    expect_error(gpc_trial_data(p, administrative = flags, look = bad), "look")
  expect_error(gpc_trial_data(p, administrative = rep(TRUE, 4), look = 1), "death")
  q <- p; q$entry[1] <- 2
  expect_error(gpc_trial_data(q, administrative = flags, look = 1), "entry.*T1")
  q <- p; q$followup[1] <- .9
  expect_error(gpc_trial_data(q, administrative = flags, look = 1), "inconsistent.*T1")
})

test_that("all channels retain exact clinical times and enforce declared cut", {
  p <- admin_subjects(); before <- serialize(p, NULL)
  times <- c(1 - .Machine$double.eps/2, 1, 1 + .Machine$double.eps)
  for (kind in c("events", "measurements", "transitions")) {
    for (i in seq_along(times)) {
      tab <- data.frame(id = "T1", time = times[i])
      if (kind == "measurements") tab$value <- 7
      if (kind == "transitions") {tab$from <- "H"; tab$to <- "U"}
      args <- list(subjects = p)
      args[[kind]] <- if (kind == "transitions") tab else list(unused = tab)
      if (kind == "transitions") args$state0 <- "state0"
      if (kind == "events") args$event <- FALSE
      if (kind == "measurements") args$measurement <- FALSE
      legacy <- do.call(gpc_trial_data, args)
      args$administrative <- rep(TRUE, 4); args$look <- 1
      if (i == 3L) {
        expect_error(do.call(gpc_trial_data, args), "canonical administrative followup.*T1")
      } else {
        got <- do.call(gpc_trial_data, args)
        expect_identical(got$channels, legacy$channels)
        expect_identical(got$transitions, legacy$transitions)
        expect_identical(got$events, legacy$events)
      }
      if (i == 1L) {
        tab$available <- 1 + .Machine$double.eps
        args[[kind]] <- if (kind == "transitions") tab else list(unused = tab)
        expect_error(do.call(gpc_trial_data, args), "availability.*look.*T1")
      }
    }
  }
  expect_identical(serialize(p, NULL), before)
})

test_that("conversion uses the design clock and refuses unresolved loss of precision", {
  p <- admin_subjects(); p$entry <- c(.1, .2, .1, .2)
  p$followup <- 2 - p$entry
  for (a in c(1, 1/365.25, 365.25)) for (b in c(0, 1000)) {
    q <- p; q$entry <- p$entry*a+b; q$followup <- p$followup*a
    cut <- gpc_trial_data(q, administrative = rep(TRUE, 4), look = 2*a+b)
    expect_identical(cut$patients$followup, 2*a+b-q$entry)
  }
  q <- p; q$entry <- rep(1e15, 4); q$followup <- rep(.9, 4)
  expect_error(gpc_trial_data(q, administrative = rep(TRUE, 4), look = 1e15+1),
               "clock origin")
  q <- admin_subjects(); q$entry <- rep(1, 4); q$followup <- rep(0, 4)
  expect_identical(gpc_trial_data(q, administrative = rep(TRUE, 4), look = 1)$patients$followup,
                   rep(0, 4))
})
