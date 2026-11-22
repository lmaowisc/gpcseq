# Tiny trial with checks worked independently in dev/evidence/user-endpoints-oracle.py.
ue_rows <- function(ids = character(), time = numeric(), value = NULL) {
  out <- data.frame(id = ids, time = as.numeric(time), stringsAsFactors = FALSE)
  if (!is.null(value)) out$value <- as.numeric(value)
  out
}

ue_patients <- function(id = c("T", "C"), followup = 2, death = FALSE) {
  data.frame(id = id, arm = rep(c(1L,0L), each = length(id)/2), entry = 0,
             followup = followup, death = death, stringsAsFactors = FALSE)
}

ue_fit <- function(cut, tiers, method = "shared", u = 2, tau = u) {
  gpc_fit(gpc_design(u, .05, method, hierarchy = gpc_hierarchy(endpoints = tiers),
                     horizons = if (method == "ipcw") tau else NULL), 1, data = cut)
}

test_that("named endpoints determine the winner in the declared clinical order", {
  # Treatment is stroke-free; control is admission-free. All other histories
  # remain fixed when a clinical protocol switches which endpoint it prioritizes.
  cut <- gpc_trial_data(ue_patients(), events = list(stroke = ue_rows("C", 1.5),
                                                      admission = ue_rows("T", 1)))
  tiers <- list(brain = gpc_endpoint("stroke","first"),
                hospital = gpc_endpoint("admission","count"))
  for (method in c("shared","ipcw")) {
    expect_equal(ue_fit(cut, tiers, method)$estimate, 1)
    expect_equal(ue_fit(cut, rev(tiers), method)$estimate, -1)
  }
})

test_that("a configured empty event channel is observed absence; an unknown channel is an error", {
  cut <- gpc_trial_data(ue_patients(), events = list(stroke = ue_rows(),
                                                      admission = ue_rows("T", 1)))
  known <- list(strokes = gpc_endpoint("stroke","count"),
                admission = gpc_endpoint("admission","first"))
  expect_equal(ue_fit(cut, known)$estimate, -1)
  expect_equal(ue_fit(cut, known[1])$estimate, 0)
  expect_error(ue_fit(cut, list(misspelled = gpc_endpoint("strkoe","count"))))
  expect_error(ue_fit(cut, list(wrong_kind = gpc_endpoint("stroke","value"))))
})

test_that("several measurement channels retain separate units, thresholds and missingness", {
  p <- ue_patients()
  walk <- ue_rows(p$id, c(1,1), c(5,6))
  kidney <- ue_rows(p$id, c(1,1), c(4,2))
  cut <- gpc_trial_data(p, measurements = list(walk = walk, creatinine = kidney))
  walking <- gpc_endpoint("walk", "value", direction = "higher", threshold = 5)
  kidney1 <- gpc_endpoint("creatinine", "value", margin = 1)
  expect_equal(ue_fit(cut, list(walk = walking, renal = kidney1))$estimate, -1)
  # Exactly on the cutoff belongs to the lower category. Both later walk
  # measurements below the cutoff tie, so the patient's actual creatinine
  # value, in its own clinical units, determines the comparison.
  walk$value <- c(5,3)
  cut <- gpc_trial_data(p, measurements = list(walk = walk, creatinine = kidney))
  expect_equal(ue_fit(cut, list(walk = walking))$estimate, 0)
  expect_equal(ue_fit(cut, list(walk = walking, renal = kidney1))$estimate, -1)
  expect_equal(ue_fit(cut, list(renal = gpc_endpoint("creatinine","value", margin = 2)))$estimate, 0)
  # One person's creatinine is unknown. A deciding tier needs both readings.
  partial <- gpc_trial_data(p, measurements = list(walk = walk, creatinine = kidney[1,,drop=FALSE]))
  expect_error(ue_fit(partial, list(walk = walking, renal = kidney1)), class = "gpc_deferred")
  walk$value <- c(7,3)
  decisive <- gpc_trial_data(p, measurements = list(walk = walk, creatinine = kidney[1,,drop=FALSE]))
  expect_equal(ue_fit(decisive, list(walk = walking, renal = kidney1))$estimate, 1)
})

test_that("custom tiers see every named channel, clipped at the pair cutoff", {
  p <- ue_patients(followup = c(4,2))
  e <- list(stroke = ue_rows(c("T","T"), c(1,3)), admission = ue_rows(c("T","C"), c(3,2)))
  m <- list(walk = ue_rows(c("T","T","C"), c(1,3,1), c(4,25,5)),
            creatinine = ue_rows(c("T","C"), c(1,1), c(2,4)))
  cut <- gpc_trial_data(p, events = e, measurements = m)
  composite <- function(treatment, control, time) {
    # The control history ends at 2; all four named channels must be clipped
    # to 2 for each patient, including the longer-observed treatment patient.
    score <- function(x) {
      stopifnot(time == 2, all(c("stroke","admission") %in% names(x$events)),
                all(c("walk","creatinine") %in% names(x$measurements)),
                all(unlist(x$events, use.names = FALSE) <= time),
                all(vapply(x$measurements,function(b) all(b$time <= time),logical(1))))
      -2*length(x$events$stroke) - length(x$events$admission) +
        tail(x$measurements$walk$value,1) - tail(x$measurements$creatinine$value,1)
    }
    sign(score(treatment) - score(control))
  }
  # Treatment: -2 + 4 - 2 = 0; control: -1 + 5 - 4 = 0.
  expect_equal(ue_fit(cut, list(composite = composite), u = 4)$estimate, 0)
  # IPCW sees the same pair when the protocol horizon is 2.
  expect_equal(ue_fit(cut, list(composite = composite), "ipcw", u = 4, tau = 2)$estimate, 0)
  # This look-specific table has a delayed stroke report. The earlier
  # calendar look cannot consume that report, even if a different tier leads.
  early <- ue_patients()
  late <- ue_rows("C", 1); late$available <- 3
  expect_error(ue_fit(gpc_trial_data(early, events=list(stroke=late)),
                      list(s = gpc_endpoint("stroke","first"))))
  expect_error(ue_fit(gpc_trial_data(early, measurements=list(walk=ue_rows("T",3,22))),
                      list(m = gpc_endpoint("walk","value"))))
})

test_that("settled deaths never reach a user-supplied clinical tier", {
  never <- function(treatment, control, time) stop("A resolved comparison reached this callback")
  tiers <- list(outcome = never)
  p <- ue_patients(followup = c(2,2), death = c(FALSE,TRUE))
  d <- gpc_trial_data(p, events = list(stroke = ue_rows("T",1)))
  for (m in c("shared","ipcw")) expect_equal(ue_fit(d, tiers, m)$estimate, 1)
  p$death <- TRUE
  d <- gpc_trial_data(p, events = list(stroke = ue_rows("T",1)))
  for (m in c("shared","ipcw")) expect_equal(ue_fit(d, tiers, m)$estimate, 0)
})

test_that("one named-data cut yields matching monitoring results under equivalent clinical rules", {
  ids <- c("A","B","C","D","E","F")
  stroke <- ue_rows(c("A","B","B","D","E"), c(3,1,4,2,4))
  admit <- ue_rows(c("A","A","B","C","D","E","F"), c(1.5,4,2,1.5,1,1.5,1.5))
  walk <- ue_rows(c(ids,"C","F"), c(rep(0,6),3,3), c(8,1,5,9,7,4,9,9))
  creatinine <- ue_rows(ids, rep(0,6), c(1,1,4,2,5,1))
  clinical <- list(strokes = gpc_endpoint("stroke","count"),
                   hospitalization = gpc_endpoint("admission","first"),
                   mobility = gpc_endpoint("walk","value",direction="higher",threshold=5),
                   renal = gpc_endpoint("creatinine","value",margin=1))
  my_comparison <- function(treatment, control, time) {
    x <- treatment; y <- control
    d <- length(y$events$stroke) - length(x$events$stroke)
    if (d != 0) return(sign(d))
    first <- function(h) if (length(h$events$admission)) min(h$events$admission) else Inf
    a <- first(x); b <- first(y)
    if (a != b) return(if (a > b) 1 else -1)
    above <- function(h) as.integer(tail(h$measurements$walk$value,1) > 5)
    d <- above(x) - above(y)
    if (d != 0) return(sign(d))
    gap <- tail(y$measurements$creatinine$value,1) - tail(x$measurements$creatinine$value,1)
    if (abs(gap) > 1) sign(gap) else 0
  }
  expected_phi <- list(c(8,-4,-4,8,-4,-4)/9, c(4,-8,4,4,4,-8)/9)
  cuts <- lapply(c(2,4), function(u) {
    known <- function(x) x[x$time <= u,,drop=FALSE]
    gpc_trial_data(ue_patients(ids, followup = u),
                   events=list(stroke=known(stroke), admission=known(admit)),
                   measurements=list(walk=known(walk), creatinine=known(creatinine)))
  })
  cvals <- NULL
  for (method in c("shared","ipcw")) for (rules in list(clinical, list(research = my_comparison))) {
    d <- gpc_design(c(2,4), c(.01,.05), method, hierarchy = gpc_hierarchy(endpoints = rules),
                    horizons = if (method == "ipcw") c(2,4) else NULL)
    for (k in 1:2) {
      d <- gpc_update(d, data = cuts[[k]])
      expect_equal(gpc_influence(d,k)$phi, expected_phi[[k]], tolerance = 1e-12)
    }
    r <- gpc_results(d)
    expect_equal(r$estimate,c(1,-1)/9)
    expect_equal(unname(gpc_covariance(d)), matrix(c(64,32,32,64),2)/243)
    expect_equal(r$action,c("continue","complete"))
    expect_equal(r$spent,c(.01,.05))
    if (is.null(cvals)) cvals <- r$critical else expect_equal(r$critical, cvals, tolerance=1e-7)
  }
})

test_that("event and measurement channels can share a clinical label without collision", {
  p <- ue_patients()
  cut <- gpc_trial_data(p, events = list(renal = ue_rows("C",1)),
                        measurements = list(renal = ue_rows(p$id, c(1,1), c(4,1))))
  count <- gpc_endpoint("renal", "count")
  value <- gpc_endpoint("renal", "value")
  expect_equal(ue_fit(cut, list(hospitalizations = count))$estimate, 1)
  expect_equal(ue_fit(cut, list(concentration = value))$estimate, -1)
  sparse <- gpc_design(2, .05, "shared", min_enrolled = 2,
                       hierarchy = gpc_hierarchy(endpoints=list(fail = gpc_endpoint("typo","count"))))
  # One patient/arm is too few to perform this test. An invalid endpoint
  # is still a configuration error, rather than a scientific deferral.
  expect_error(gpc_update(sparse, data = cut))
  held <- gpc_design(2, .05, "shared", min_enrolled = 2,
                     hierarchy = gpc_hierarchy(endpoints=list(ok = count)))
  expect_equal(gpc_results(gpc_update(held,data=cut))$action,"final_untested")
})

test_that("unselected channels respect the cut and transport selection does not select a trial rule", {
  p <- ue_patients()
  adm <- ue_rows("C",1); adm$available <- 1
  stroke <- ue_rows("T",1); stroke$available <- 3
  tiers <- list(admit = gpc_endpoint("admission","first"))
  expect_error(ue_fit(gpc_trial_data(p, events=list(admission=adm,stroke=stroke), event="admission"),
                      tiers))
  stroke$available <- 1
  ambiguous <- gpc_trial_data(p, events=list(admission=adm,stroke=stroke))
  expect_equal(ue_fit(ambiguous, tiers)$estimate, 1)
  expect_error(gpc_transport_outcome(ambiguous, 2))
  selected <- gpc_trial_data(p, events=list(admission=adm,stroke=stroke), event="admission")
  out <- gpc_transport_outcome(selected, 2)
  expect_equal(out$hospitalization_event, c(0,1))
  expect_equal(out$hospitalization_time, c(2,1))
  # Explicitly choosing which channel populates the transport outcome must
  # not silently remove the other named channel from a sequential protocol.
  expect_equal(ue_fit(selected, list(s=gpc_endpoint("stroke","first")))$estimate,-1)
})
