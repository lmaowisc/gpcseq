test_that("the longest supported horizon includes exact fractional thresholds in either arm", {
  # In the limiting arm 14/100 enrollees have reached time 87, but
  # multiplying the double .14 by 100 gives slightly more than 14.
  for (limiting in 0:1) {
    opportunity <- c(seq_len(100L), rep(100L, 73L))
    x <- data.frame(id = sprintf("p%03d", seq_along(opportunity)),
                    arm = rep(c(limiting, 1L - limiting), c(100L, 73L)),
                    entry = 100 - opportunity, followup = opportunity,
                    death = FALSE)
    for (p in c(.14 - .Machine$double.eps, .14, .14 + .Machine$double.eps)) {
      candidates <- sort(unique(opportunity))
      meets_rule <- vapply(candidates, function(t) all(vapply(0:1, function(a)
        mean(opportunity[x$arm == a] >= t) >= p, logical(1L))), logical(1L))
      plan <- gpc_design(100, .05, "ipcw", maturity = list(p = p, cap = 100))
      expect_equal(gpc_fit(plan, 1, x)$horizon, max(candidates[meets_rule]))
    }
  }
})
