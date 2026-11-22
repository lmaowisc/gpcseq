test_that("risk counts and ranks reproduce weighted patient and grid totals", {
  # A death and administrative censoring can coincide. Multiple deaths and
  # recurrent events at the horizon exercise the inclusive time convention.
  p <- data.frame(id=letters[1:12], arm=rep(1:0,each=6), entry=0,
                  followup=c(0,1,2,2,3,3,0,1,2,2,3,3),
                  death=as.logical(c(1,0,1,0,0,1,0,1,1,0,1,0)))
  e <- data.frame(id=c("d","e","j","k","f","f"),time=c(2,0,2,0,0,2),
                  type="hospitalization")
  p <- p[c(7,1,9,2,10,3,11,4,12,5,8,6), ]
  ti <- which(p$arm==1)
  ci <- which(p$arm==0)
  rule <- gpc_hierarchy()
  for (tau in c(0,1,2,2.5)) {
    km <- lapply(list(ti,ci), function(i) .pair_km(p$followup[i],p$death[i],tau))
    actual <- .pair_summaries(p,e,ti,ci,tau,TRUE,km,NULL,rule)
    expected <- .pair_summaries_reference(p,e,ti,ci,tau,TRUE,km,NULL,rule)
    expect_equal(actual,expected,tolerance=1e-12)
    for (contrast in list(c(1,0),c(0,1),c(1,-1))) {
      expect_equal(.if_ipcw_linear(actual,contrast,p,list(ti,ci),km),
                   .if_ipcw_linear(expected,contrast,p,list(ti,ci),km),tolerance=1e-12)
    }
  }
})

test_that("only the built-in restricted count and last-event rule uses ranking", {
  h <- gpc_hierarchy()
  expect_true(.rank_eligible(TRUE,NULL,h))
  expect_false(.rank_eligible(FALSE,NULL,h))
  expect_false(.rank_eligible(TRUE,1,h))
  for (rule in list(NULL,gpc_hierarchy(admission_count=FALSE),
                    gpc_hierarchy(admission_time="none"),gpc_hierarchy(admission_time="first"),
                    gpc_hierarchy(biomarker_margin=0),gpc_hierarchy(biomarker_threshold=5),
                    gpc_hierarchy(endpoints=list(count=gpc_endpoint("admission","count"))),
                    gpc_hierarchy(compare=function(treatment,control,time) 0))) {
    expect_false(.rank_eligible(TRUE,NULL,rule))
  }
})

test_that("the ranked calculation retains pair counts beyond 32-bit integers", {
  # All controls die at the horizon. Every treated patient is observed alive.
  n <- 50000L
  p <- data.frame(arm=rep(1:0,each=n),followup=2,death=rep(c(FALSE,TRUE),each=n))
  ti <- seq_len(n)
  ci <- n+ti
  km <- lapply(list(ti,ci),function(i) .pair_km(p$followup[i],p$death[i],2))
  z <- .rank_ipcw(p,NULL,ti,ci,2,km)
  expect_equal(z$n_pairs,2.5e9)
  expect_equal(z$n_ascertained,2.5e9)
  expect_equal(z$mean,c(1,0))
  expect_true(all(z$average[,1L]==1 & z$average[,2L]==0))
  expect_equal(as.numeric(z$mass),c(1,0))
})
