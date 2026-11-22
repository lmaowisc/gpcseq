# The public fit is tested with hand-scored pairs, then with full test/IF
# workflows. Fractions for the latter are independently reproduced by
# python3 dev/evidence/hierarchy-oracle.py (case-weighted full refits).
hier_admissions <- function(id, times) {
  data.frame(id = rep(id, length(times)), time = as.numeric(times),
             type = rep("hospitalization", length(times)), value = rep(NA_real_, length(times)))
}

hier_pair <- function(rule = NULL, ht = numeric(), hc = numeric(),
                      bt = NA_real_, bc = NA_real_, followup = c(4, 4),
                      death = c(FALSE, FALSE), legacy = FALSE) {
  p <- data.frame(id = c("T", "C"), arm = c(1L, 0L), entry = 0,
                  followup = followup, death = death)
  e <- rbind(hier_admissions("T", ht), hier_admissions("C", hc))
  vals <- c(bt, bc)
  m <- !is.na(vals)
  if (any(m)) e <- rbind(e, data.frame(id = p$id[m], time = 0,
                                       type = "biomarker", value = vals[m]))
  d <- if (is.null(rule)) gpc_design(4, .05, "shared") else
    gpc_design(4, .05, "shared", hierarchy = rule)
  if (legacy) d$design$hierarchy <- NULL  # field absent in older saved designs
  gpc_fit(d, 1, p, e)$estimate
}

test_that("count precedes timing; later first and last need not favor the same patient", {
  last <- gpc_hierarchy(admission_time = "last")
  first <- gpc_hierarchy(admission_time = "first")
  count <- gpc_hierarchy(admission_time = "none")
  t <- c(1, 4); c0 <- c(2, 3)
  expect_equal(hier_pair(NULL, t, c0), 1)
  expect_equal(hier_pair(gpc_hierarchy(), t, c0), 1)
  expect_equal(hier_pair(last, t, c0), 1)
  expect_equal(hier_pair(first, t, c0), -1)
  expect_equal(hier_pair(count, t, c0), 0)
  expect_equal(hier_pair(NULL, t, c0, legacy = TRUE), 0)
  # A single earlier admission is still better than two later admissions.
  expect_equal(hier_pair(last, 1, c(2, 3)), 1)
  expect_equal(hier_pair(first, 1, c(2, 3)), 1)
})

test_that("an admission exactly at the window counts; absence in both arms falls through", {
  for (which_time in c("first", "last")) {
    rule <- gpc_hierarchy(admission_count = FALSE, admission_time = which_time,
                          biomarker_margin = 0)
    expect_equal(hier_pair(rule, hc = 4), 1)
    expect_equal(hier_pair(rule, ht = 4), -1)
    expect_equal(hier_pair(rule, bt = 2, bc = 3), 1)
    expect_equal(hier_pair(rule, bt = 2, bc = 2), 0)
  }
})

test_that("the first death settles even a hostile living comparison", {
  never <- function(treatment, control, time) stop("A settled pair reached the callback")
  rule <- gpc_hierarchy(compare = never)
  expect_equal(hier_pair(rule, ht = c(.2, .3), followup = c(4, 2),
                         death = c(FALSE, TRUE)), 1)
  expect_equal(hier_pair(rule, hc = c(.2, .3), followup = c(2, 4),
                         death = c(TRUE, FALSE)), -1)
  expect_equal(hier_pair(rule, ht = 1, followup = c(2, 2),
                         death = c(TRUE, TRUE)), 0)
})

test_that("marker direction, clinical cutoff and pairwise margin have different meanings", {
  for (direction in c("lower", "higher")) {
    s <- if (direction == "lower") 1 else -1
    rule <- gpc_hierarchy(biomarker_margin = 1, biomarker_direction = direction)
    expect_equal(hier_pair(rule, bt = 5, bc = 7), s)
    expect_equal(hier_pair(rule, bt = 5, bc = 6), 0)
    expect_error(hier_pair(rule, bt = 7), class = "gpc_deferred")
    expect_equal(hier_pair(rule, hc = 3), 1) # missing measurements cannot override a decision
    cutoff <- gpc_hierarchy(biomarker_threshold = 5, biomarker_direction = direction)
    expect_equal(hier_pair(cutoff, bt = 5, bc = 5.01), s) # == is in the lower category
    expect_equal(hier_pair(cutoff, bt = 6, bc = 12), 0)
    expect_equal(hier_pair(cutoff, bt = 1, bc = 5), 0)
    expect_error(hier_pair(cutoff, bc = 7), class = "gpc_deferred")
  }
  expect_error(gpc_hierarchy(biomarker_threshold = 5, biomarker_margin = 1))
})

test_that("a living callback receives both histories truncated to the pair window", {
  seen <- new.env(parent = emptyenv())
  seen$n <- 0L
  check <- function(treatment, control, time) {
    stopifnot(is.numeric(treatment$admissions), is.numeric(control$admissions),
              is.data.frame(treatment$biomarker), is.data.frame(control$biomarker),
              all(c("time","value") %in% names(treatment$biomarker)),
              all(c("time","value") %in% names(control$biomarker)),
              time == 2, identical(as.numeric(treatment$admissions), 1),
              identical(as.numeric(control$admissions), 2),
              identical(as.numeric(treatment$biomarker$value), 3),
              identical(as.numeric(control$biomarker$value), 8))
    seen$n <- seen$n + 1L
    sign(control$biomarker$value - treatment$biomarker$value)
  }
  p <- data.frame(id = c("T", "C"), arm = c(1, 0), entry = 0,
                  followup = c(4, 2), death = FALSE)
  e <- rbind(hier_admissions("T", c(1, 3)), hier_admissions("C", 2),
             data.frame(id = c("T", "T", "C"), time = c(1, 3, 2),
                        type = "biomarker", value = c(3, 18, 8)))
  plan <- gpc_design(4, .05, "shared", hierarchy = gpc_hierarchy(compare = check))
  expect_equal(gpc_fit(plan, 1, p, e)$estimate, 1)
  expect_gt(seen$n, 0L)
})

hier_count_last <- function(treatment, control, time) {
  t <- treatment$admissions; c0 <- control$admissions
  if (length(t) != length(c0)) return(sign(length(c0) - length(t)))
  if (!length(t)) return(0)
  sign(max(t) - max(c0))
}

test_that("equivalent hierarchies give the same two-look influences, covariance and actions", {
  h <- list(c(.6, 2.5), c(1.2, 2.1), numeric(), c(.8, 3.8),
            c(.7, 2.3), c(1.3, 2.7), numeric(), c(1.5, 2.9))
  ids <- paste0(rep(c("T", "C"), each = 4), rep(1:4, times = 2))
  history <- do.call(rbind, Map(hier_admissions, ids, h))
  phi <- list(c(-11, -3, 17, -3, 13, -3, -7, -3)/16,
              c(-7, -15, 13, 9, 9, 1, -11, 1)/16)
  baseline <- NULL
  for (method in c("shared", "ipcw")) {
    for (rule in list(gpc_hierarchy(admission_time = "last"),
                      gpc_hierarchy(compare = hier_count_last))) {
      d <- gpc_design(c(2, 4), c(.01, .05), method,
                      horizons = if (method == "ipcw") c(2, 4) else NULL, hierarchy = rule)
      for (k in 1:2) {
        u <- 2*k
        p <- data.frame(id = ids, arm = rep(c(1, 0), each = 4), entry = 0,
                        followup = u, death = FALSE)
        d <- gpc_update(d, p, history[history$time <= u, , drop = FALSE])
        expect_equal(gpc_influence(d, k)$phi, phi[[k]], tolerance = 1e-11)
      }
      log <- gpc_results(d)
      expect_equal(log$estimate, c(-5, -1)/16, tolerance = 1e-12)
      expect_equal(unname(gpc_covariance(d)), matrix(c(83, 63, 63, 91), 2)/512,
                   tolerance = 1e-11)
      expect_equal(log$action, c("continue", "complete"))
      expect_equal(log$spent, c(.01, .05))
      if (is.null(baseline)) baseline <- log$critical else
        expect_equal(log$critical, baseline, tolerance = 1e-7)
    }
  }
})

test_that("restricted timing uses first-death ascertainment and both estimated KM weights", {
  id <- c("T1", "T2", "T3", "C1", "C2", "C3", "C4")
  p <- data.frame(id = id, arm = c(1, 1, 1, 0, 0, 0, 0), entry = 0,
                  followup = c(1, 4, 4, 1, 4, 4, 3),
                  death = c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE))
  h <- list(numeric(), c(1, 3), c(2, 4), numeric(), c(1, 4), c(2, 3), c(1, 2))
  e <- do.call(rbind, Map(hier_admissions, id, h))
  for (rule in list(gpc_hierarchy(admission_time = "last"),
                    gpc_hierarchy(compare = hier_count_last))) {
    plan <- gpc_design(4, .05, "ipcw", horizons = 4, hierarchy = rule)
    f <- gpc_fit(plan, 1, p, e)
    joint <- gpc_ipcw_effects(plan, 1, p, e)
    expect_equal(f$diagnostics$n_ascertained, 7)
    expect_equal(f$diagnostics$win, 7/16)
    expect_equal(f$diagnostics$loss, 3/16)
    expect_equal(f$estimate, 1/4)
    expect_equal(gpc_influence(f)$phi, c(0, -9/16, 9/16, 3/4, -1, 1/2, -1/4))
    expect_equal(joint$influences$net_benefit, gpc_influence(f)$phi)
    expect_equal(joint$summary$natural, c(1/4, 5/3, 7/3))
  }
})
