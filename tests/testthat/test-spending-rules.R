test_that("the four named rules reproduce independently tabulated cumulative spending", {
  # Fractions are 1/4, 1/2, 1. Expected values use a separate normal CDF
  # implementation or high-precision arithmetic for the non-normal families.
  u <- c(3, 6, 12)
  obf <- gpc_design(u)$design
  expect_identical(obf$spending_rule, "obrien-fleming")
  expect_identical(obf$spending_fractions, c(.25, .5, 1))
  expect_equal(obf$spending, c(.00008857543832140365, .005574596680784402, .05),
               tolerance = 1e-12)
  pc <- gpc_design(u, spending = "pocock")$design
  expect_equal(pc$spending, c(.017868700975439425, .031005725347913876, .05),
               tolerance = 1e-12)
  pow <- gpc_design(u, "power", spending_param = 2)$design
  expect_equal(pow$spending, c(.003125, .0125, .05), tolerance = 1e-14)
  hsd <- gpc_design(u, "hsd")$design
  expect_identical(hsd$spending_param, -4)
  expect_equal(hsd$spending, c(.0016029301640042494, .005960146101105879, .05),
               tolerance = 1e-12)
  front <- gpc_design(u, "hsd", spending_param = 3.5)$design
  expect_equal(front$spending, c(.030064776602651187, .04259764009841553, .05),
               tolerance = 1e-12)
  linear <- c(.0125, .025, .05)
  for (rule in c("power", "linear")) {
    fit <- gpc_design(u, rule)$design
    expect_equal(fit$spending, linear)
    expect_identical(fit$spending_rule, "power")
    expect_identical(fit$spending_param, 1)
  }
  for (g in c(0, 1e-20, -1e-20))
    expect_equal(gpc_design(u, "hsd", spending_param=g)$design$spending, linear)
})

test_that("sidedness sets the total level, and the horizon does not set the spending clock", {
  u <- c(1, 2, 4); f <- c(.25, .5, 1)
  upper <- gpc_design(u, direction = "upper")$design
  expect_identical(upper$level, .025)
  expect_equal(upper$spending, c(.00000736680843586949, .0015253227579889112, .025),
               tolerance = 1e-12)
  for (direction in c("two-sided", "upper-final-two-sided"))
    expect_equal(gpc_design(u, direction=direction, level=.025)$design$spending, upper$spending)
  expect_identical(gpc_design(u, direction="upper-final-two-sided")$design$level, .05)
  restricted <- gpc_design(c(1,3,4), "obrien-fleming", "ipcw",
                           maturity = list(p=.25, cap=4), spending_fractions=f)
  expect_equal(restricted$design$spending, gpc_design(u)$design$spending)
  expect_equal(gpc_design(4, "pocock", level=.04)$design$spending, .04)
  # An early one-sided look uses its entire scheduled probability in the upper tail.
  x <- data.frame(id=c("T1","T2","C1","C2"), arm=c(1,1,0,0), entry=0,
                  followup=c(.25,1,.25,1), death=c(TRUE,FALSE,TRUE,FALSE))
  z <- gpc_results(gpc_update(gpc_design(u, direction="upper"), x))
  expect_identical(z$action, "continue")
  expect_equal(z$critical, qnorm(upper$spending[1L], lower.tail=FALSE))
})

test_that("explicit spending stays authoritative and deferred looks do not spend", {
  v <- c(.01, .05)
  z <- gpc_design(1:2, v, "shared")
  expect_identical(z$design$spending, v)
  expect_identical(z$design$spending_rule, "custom")
  expect_null(z$design$spending_fractions)
  expect_identical(z$design$level, .05)
  for (arg in list(list(level=.05), list(spending_fractions=c(.5,1)), list(spending_param=1)))
    expect_error(do.call(gpc_design, c(list(1:2, v), arg)), "numeric spending")
  named <- gpc_design(1:2, "pocock")
  custom <- gpc_design(1:2, named$design$spending)
  for (k in 1:2) {
    d <- gpc_example(k)
    named <- gpc_update(named, d$patients, d$events)
    custom <- gpc_update(custom, d$patients, d$events)
  }
  expect_equal(gpc_results(named), gpc_results(custom))
  expect_equal(gpc_results(named)$spent, c(0,.05))
})

test_that("invalid or numerically unresolved schedules are diagnosed at design", {
  for (f in list(NA_real_, c(0,1), c(.4,.8), c(1,1), c(.2,Inf), matrix(c(.5,1))))
    expect_error(gpc_design(1:2, spending_fractions=f), "spending_fractions")
  for (a in list(NA_real_, Inf, c(.01,.05), -1, 1))
    expect_error(gpc_design(1:2, level=a), "level")
  for (r in list("p", "unknown", NA_character_, c("power","hsd")))
    expect_error(gpc_design(1:2, r), "spending")
  for (r in c("pocock", "obrien-fleming", "linear"))
    expect_error(gpc_design(1:2, r, spending_param=1), "spending_param")
  for (s in list(0, -1, NA_real_, c(1,2)))
    expect_error(gpc_design(1:2, "power", spending_param=s), "spending_param")
  expect_error(gpc_design(1:2, "hsd", spending_param=41), "spending_param")
  # The planned 10% look is well below the default tolerance but usable.
  x <- data.frame(id=c("T1","T2","C1","C2"), arm=c(1,1,0,0), entry=0,
                  followup=c(.25,1,.25,1), death=c(TRUE,FALSE,TRUE,FALSE))
  early <- gpc_design(c(1,10))
  expect_equal(early$design$spending[1] / 5.720318790425403e-10, 1, tolerance=1e-11)
  first <- gpc_results(gpc_update(early,x))
  expect_equal(first$critical, 6.1979503230456, tolerance=1e-10)
  expect_identical(first$action, "continue")
  # An actual underflow or two numerically indistinguishable joint targets.
  expect_error(gpc_design(c(1,1000)), "look 1.*numerical precision")
  expect_error(gpc_design(c(1,2,100)), "look 2.*numerical precision")
})

test_that("a positive tiny first spend uses its tail quantile through the public workflow", {
  x <- data.frame(id=c("T1","T2","C1","C2"), arm=c(1,1,0,0), entry=0,
                  followup=c(.25,1,.25,1), death=c(TRUE,FALSE,TRUE,FALSE))
  for (plan in list(gpc_design(1:20), gpc_design(c(1,100)),
                    gpc_design(c(1,100), direction="upper-final-two-sided"))) {
    alpha <- plan$design$spending[1L]
    result <- gpc_results(gpc_update(plan, x))
    expect_gt(alpha, 0)
    expect_identical(1-alpha, 1)
    expect_identical(result$action, "continue")
    expect_identical(result$spent, alpha)
    sides <- if (result$direction == "two-sided") 2 else 1
    expect_equal(sides * pnorm(-result$critical) / alpha, 1, tolerance=2e-12)
  }
  expect_equal(gpc_results(gpc_update(gpc_design(1:20), x))$critical,
               8.765225405765818, tolerance=2e-12)
})
