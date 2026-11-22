test_that("calendar and reporting dates are visible and enrollment is immutable", {
  plan <- gpc_design(c(1,2),c(.01,.05),"shared")
  cut <- gpc_example(1); err <- cut$events; err$available <- 2
  expect_error(gpc_update(plan,cut$patients,err),"availability")
  pp <- cut$patients; pp$followup[1] <- 2
  expect_error(gpc_update(plan,pp,cut$events),"followup")
  early <- gpc_update(plan,cut$patients,cut$events)
  late <- gpc_example(2); late$patients$arm[1] <- 0L
  expect_error(gpc_update(early,late$patients,late$events),"arm changed")
  expect_equal(names(gpc_results(plan)),names(gpc_results(early)))
  expect_equal(nrow(gpc_results(plan)),0L)
})

test_that("maturity includes all enrolled and counts empirical quantile ties", {
  design <- gpc_design(c(10,12),c(.01,.05),"ipcw",maturity=list(p=.5,cap=30))
  p <- data.frame(arm=c(0,0,1,1,1,1,1),entry=c(1,3,0,2,11,11,12))
  expect_equal(.gpc_horizon(design$design,1,p[p$entry<=10,]),9)
  expect_equal(.gpc_horizon(design$design,2,p),1)
})

test_that("a biomarker is not needed for a pair decided by survival", {
  des <- gpc_design(2,.05,"shared",biomarker_margin=0)
  pp <- data.frame(id=c("T","C"),arm=c(1,0),entry=0,followup=c(2,1),death=c(FALSE,TRUE))
  ev <- data.frame(id="T",time=0,type="biomarker",value=NA_real_)
  f <- gpc_fit(des,1,pp,ev)
  expect_equal(f$estimate,1)
  expect_equal(gpc_influence(f)$phi,c(0,0))
  pp$death[2] <- FALSE
  expect_error(gpc_fit(des,1,pp,ev),"Missing biomarker")
})

test_that("a deferred look spends nothing and the final closes monitoring", {
  d <- gpc_design(c(1,2),c(.01,.05),"ipcw","win_ratio",horizons=c(1,2))
  a <- gpc_example(1); b <- gpc_example(2)
  x <- gpc_update(d,a$patients,a$events)
  expect_equal(gpc_results(x)$spent,0)
  expect_error(gpc_influence(x,1),"No completed")
  x <- gpc_update(x,b$patients,b$events)
  expect_equal(gpc_results(x)$natural[2],3)
  expect_error(gpc_update(x,b$patients,b$events),"after stopping")
})
