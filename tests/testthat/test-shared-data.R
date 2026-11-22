test_that("one mapped cut reproduces native fits, influences, and decisions", {
  # Exercise different names, storage types, and independent row order in
  # the clinical table. Extra baseline predictors do not enter the endpoint.
  build <- function(k) {
    raw <- gpc_example(k)
    p <- raw$patients; p$age <- c(54,65,47,70)
    names(p)[match(c("id","arm","death"), names(p))] <- c("pid","a","d")
    p$d <- as.integer(p$d)
    gpc_trial_data(p, id="pid", treatment="a", death="d", state0="state0",
                   events=list(admit=raw$events[,c("id","time")]),
                   transitions=raw$transitions)
  }
  x <- gpc_example(2); cut <- build(2)
  for (m in c("shared","ipcw","aj")) {
    plan <- switch(m, shared=gpc_design(2,.05,"shared"),
                   ipcw=gpc_design(2,.05,"ipcw",horizons=2),
                   aj=gpc_design(2,.05,"aj",horizons=2,
                                 state_order=c("well","hospitalized","dead"),
                                 absorbing="dead",death_state="dead"))
    native <- gpc_fit(plan,1,x$patients,x$events,x$transitions)
    got <- gpc_fit(plan,1,data=cut)
    expect_equal(got$estimate,native$estimate)
    expect_equal(got$subjects,native$subjects)
    expect_equal(gpc_fit(plan,1,data=x$patients,events=x$events,
                         transitions=x$transitions)$estimate, native$estimate)
  }
  plan <- gpc_design(2,.05,"ipcw",horizons=2)
  expect_equal(gpc_ipcw_effects(plan,1,data=cut),
               gpc_ipcw_effects(plan,1,x$patients,x$events))
  a <- b <- gpc_design(c(1,2),c(.01,.05),"shared")
  for (k in 1:2) {
    raw <- gpc_example(k)
    a <- gpc_update(a,data=build(k))
    b <- gpc_update(b,raw$patients,raw$events)
  }
  expect_equal(gpc_results(a),gpc_results(b))
  expect_equal(gpc_covariance(a),gpc_covariance(b))
})

test_that("data adapters preserve the cut and never guess the clinical channel", {
  raw <- gpc_example(2); p <- raw$patients; e <- raw$events[,c("id","time")]
  plan <- gpc_design(c(1,2),c(.01,.05))
  cut <- gpc_trial_data(p, events=list(admit=e))
  expect_error(gpc_update(plan, data=cut),"beyond the look")
  expect_error(gpc_update(plan,p,data=cut),"both")
  expect_error(gpc_update(plan,data=cut,events=e),"histories")
  ambiguous <- gpc_trial_data(p,events=list(admit=e,fall=e))
  expect_error(gpc_fit(gpc_design(2,.05),1,data=ambiguous),"Select event")
  expect_equal(gpc_trial_data(p,events=list(admit=e,fall=e),event="admit")$events, cut$events)
  expect_error(gpc_trial_data(transform(p,arm=factor(arm))),"0/1")
  early <- gpc_example(1); late <- transform(e,available=1.5)
  delayed <- gpc_trial_data(early$patients, events=list(admit=late))
  expect_error(gpc_update(plan,data=delayed), "availability")
  e$id <- "not-enrolled"
  expect_error(gpc_trial_data(p,events=list(admit=e)))
})

test_that("a gpclearn outcome representation has the same keyed clinical data", {
  x <- gpc_example(2)
  p <- x$patients[c(3,1,4,2),]; p$pid <- p$id; p$id <- NULL
  e <- list(admit=x$events[,c("id","time")])
  # A contract fixture; dev/check-interop.R constructs the actual package class.
  y <- structure(list(subjects=p, id_column="pid", ids=p$pid,
                      end=p$followup,death=as.integer(p$death),
                      events=e,measurements=list(),states=list()),class="gpc_outcome")
  expect_identical(gpc_trial_data(y),gpc_trial_data(p,id="pid",events=e))
  y$ids <- rev(y$ids)
  expect_error(gpc_trial_data(y),"order")
})

test_that("printed scientific notation from a numeric subject key remains joined", {
  p <- data.frame(key=c(1e6,2e6), followup=2, death=0L, arm=0:1,entry=0)
  e <- data.frame(id="1e+06",time=.5)
  z <- structure(list(subjects=p,id_column="key",ids=c("1e+06","2e+06"),
                      end=p$followup,death=p$death,events=list(admit=e),
                      measurements=list(),states=list()),class="gpc_outcome")
  d <- gpc_trial_data(z)
  expect_identical(d$patients$id, c("1000000","2000000"))
  expect_identical(d$events$id, "1000000")
})

test_that("a later named cut cannot switch the endpoint under a preset", {
  one <- gpc_example(1); two <- gpc_example(2)
  a <- list(admission=one$events[,c("id","time")])
  b <- list(stroke=two$events[,c("id","time")])
  design <- gpc_update(gpc_design(c(1,2),c(.01,.05)),
                       data=gpc_trial_data(one$patients,events=a))
  expect_error(gpc_update(design,data=gpc_trial_data(two$patients,events=b)),
               "previous look")
  names(b) <- "admission"
  expect_equal(gpc_results(gpc_update(design,data=gpc_trial_data(two$patients,events=b)))$look,1:2)
})

test_that("custom callback aliases retain their clinical meaning between looks", {
  one <- gpc_example(1); two <- gpc_example(2)
  my_rule <- gpc_hierarchy(compare=function(treatment, control, time)
    sign(length(control$admissions)-length(treatment$admissions)))
  d <- gpc_design(c(1,2), c(.01,.05), hierarchy=my_rule)
  a <- list(admission=one$events[,c("id","time")])
  b <- list(stroke=two$events[,c("id","time")])
  d <- gpc_update(d, data=gpc_trial_data(one$patients,events=a))
  expect_error(gpc_update(d,data=gpc_trial_data(two$patients,events=b)), "previous look")
})

test_that("known reporting dates survive combining dated and undated channels", {
  raw <- gpc_example(1)
  adm <- transform(raw$events[,c("id","time")],available=1.5)
  marker <- data.frame(id="T1",time=.25,value=1)
  cut <- gpc_trial_data(raw$patients,events=list(admission=adm),measurements=list(marker=marker))
  expect_equal(cut$events$available,c(1.5,NA_real_))
  d <- gpc_design(1,.05)
  expect_error(gpc_fit(d,1,data=cut),"availability")
  expect_error(gpc_fit(d,1,cut$patients,cut$events))
})
