test_that("custom scores see cut-specific histories at repeated and varying cutoffs", {
  p <- data.frame(id=c("a","b","c","d"),arm=c(1,1,0,0),entry=0,
                  followup=c(4,3,4,2),death=FALSE)
  visits <- data.frame(id=rep(p$id,each=3), time=c(1,2,4,1,2,3,1,2,4,0,1,2))
  values <- transform(visits,value=time)
  cut <- gpc_trial_data(p,events=list(visit=visits),measurements=list(marker=values))
  observed <- numeric()
  scorer <- function(treatment,control,time) {
    observed <<- c(observed,time)
    for (x in list(treatment,control)) {
      expect_true(all(x$events$visit <= time))
      expect_equal(x$measurements$marker$value,x$events$visit)
    }
    sign(length(control$events$visit)-length(treatment$events$visit))
  }
  rule <- gpc_hierarchy(endpoints=list(custom=scorer))
  for (h in c(2,3)) {
    observed <- numeric()
    plan <- gpc_design(4,.05,"ipcw",horizons=h,hierarchy=rule)
    gpc_fit(plan,1,data=cut)
    expect_equal(observed,rep(h,if (h==2) 4 else 2))
  }
  observed <- numeric()
  gpc_fit(gpc_design(4,.05,"shared",hierarchy=rule),1,data=cut)
  expect_equal(sort(observed),c(2,2,3,4))
})
