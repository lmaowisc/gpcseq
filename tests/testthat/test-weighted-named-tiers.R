test_that("named last-event tiers and callbacks preserve all fitted-KM influences", {
  id <- c("T1","T2","T3","C1","C2","C3","C4")
  p <- data.frame(id = id, arm = c(1,1,1,0,0,0,0), entry = 0,
                  followup = c(1,4,4,1,4,4,3),
                  death = c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE))
  times <- list(numeric(), c(1,3), c(2,4), numeric(), c(1,4), c(2,3), c(1,2))
  strokes <- data.frame(id = rep(id, lengths(times)), time = unlist(times))
  # Reverse event order and interleave the two arms. Compare influences by ID.
  order <- c(7,1,4,3,6,2,5)
  cut <- gpc_trial_data(p[order, ], events = list(stroke = strokes[rev(seq_len(nrow(strokes))), ]),
                        event = "stroke")
  custom <- function(treatment, control, time) {
    a <- treatment$events$stroke
    b <- control$events$stroke
    if (length(a) != length(b)) return(sign(length(b) - length(a)))
    if (length(a)) sign(max(a) - max(b)) else 0
  }
  rules <- list(gpc_hierarchy(),
                gpc_hierarchy(endpoints = list(count = gpc_endpoint("stroke","count"),
                                               last = gpc_endpoint("stroke","last"))),
                gpc_hierarchy(endpoints = list(stroke = custom)))
  # Exact derivatives from dev/evidence/named-weighted-refits.py.
  phi <- list(net_benefit = c(0,-9/16,9/16,3/4,-1,1/2,-1/4),
              win_odds = c(0,-6/5,6/5,8/5,-32/15,16/15,-8/15),
              win_ratio = c(0,-15/7,15/7,16/7,-76/21,44/21,-16/21))
  for (rule in rules) {
    plan <- gpc_design(4, .05, "ipcw", horizons = 4, hierarchy = rule)
    joint <- gpc_ipcw_effects(plan, 1, data = cut)
    expect_equal(joint$diagnostics$n_ascertained, 7)
    expect_equal(unname(joint$diagnostics$G_at_horizon), c(2/3,2/3))
    expect_equal(joint$summary$effect, names(phi))
    expect_equal(joint$summary$natural, c(1/4,5/3,7/3))
    expect_equal(joint$summary$se^2, c(3/16,64/75,52/21))
    expect_equal(joint$influences$id, id[order])
    for (effect in names(phi)) {
      expect_equal(joint$influences[[effect]], phi[[effect]][order])
      single <- gpc_design(4, .05, "ipcw", effect = effect, horizons = 4, hierarchy = rule)
      fit <- gpc_fit(single, 1, data = cut)
      expect_equal(gpc_influence(fit)$phi, phi[[effect]][order])
    }
  }
})
