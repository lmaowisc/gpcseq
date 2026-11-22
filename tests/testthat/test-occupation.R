test_that("occupation curves retain tied clinical jumps at the horizon and respect each arm", {
  p <- data.frame(id=c(paste0("t",1:4),paste0("c",1:4)), arm=rep(1:0,each=4L),
                  entry=0, followup=c(2,2,1,2,2,2,2,.5), death=FALSE,
                  state0=c("U","H","H","H",rep("H",4)))
  p$death[4] <- TRUE
  e <- data.frame(id=c("t1","t2","t4","t2","c3"), time=c(1,1,2,2,1.5),
                  from=c("U","H","H","U","H"),to=c("H","U","D","H","U"))
  plan <- gpc_design(2,.05,"aj",horizons=2,state_order=c("H","U","D"),
                     absorbing="D",death_state="D")
  fit <- gpc_fit(plan,1,p[c(7,2,8,6,4,3,5,1),],transitions=e[5:1,],keep_occupation=TRUE)
  plain <- gpc_fit(plan,1,p[c(7,2,8,6,4,3,5,1),],transitions=e[5:1,])
  stripped <- fit; stripped$occupation <- NULL
  expect_identical(stripped,plain)
  expect_error(gpc_occupation(plain),"No curves stored")
  # At 1, two treatment patients exchange states before censoring in H.
  # At 2, one of two H patients dies and the sole U patient returns to H.
  # In control, censoring at .5 leaves three risks at the jump to U at 1.5.
  times <- c(2,0,1,1.75,2)
  x <- gpc_occupation(fit, times=times)
  expect_equal(x$time,rep(rep(times,each=3),2))
  expect_equal(x$state,rep(c("H","U","D"),10))
  # Named and ordinary elapsed-time grids follow the same clinical comparisons.
  named <- gpc_occupation(fit, times = setNames(times, letters[seq_along(times)]))
  expect_identical(x, named)
  expect_identical(gpc_occupation(fit, times = integer())$time, integer())
  want1 <- rbind(c(.625,0,.375),c(.75,.25,0),c(.75,.25,0),c(.75,.25,0),c(.625,0,.375))
  want0 <- rbind(c(2/3,1/3,0),c(1,0,0),c(1,0,0),c(2/3,1/3,0),c(2/3,1/3,0))
  expect_equal(x$probability[x$arm==1],as.vector(t(want1)))
  expect_equal(x$probability[x$arm==0],as.vector(t(want0)))
  raw <- gpc_occupation(fit)
  expect_equal(unique(raw$time[raw$arm==0]),c(0,1.5,2))
  expect_equal(unique(raw$time[raw$arm==1]),c(0,1,2))
  done <- gpc_update(plan,p,transitions=e,keep_occupation=TRUE)
  expect_equal(gpc_occupation(done,look=1,times=times),x)
  expect_error(gpc_occupation(fit,look=2),"does not match")
  expect_error(gpc_occupation(done,look=2),"No completed test")
  expect_error(gpc_occupation(plan),"No completed test")
  expect_error(gpc_occupation(fit,times=2.1),"fitted horizon")
  expect_equal(nrow(gpc_occupation(fit,times=numeric())),0)
})

test_that("occupation curves require an AJ fit and an explicit request for retention", {
  p <- gpc_example(2)
  plan <- gpc_design(2,.05,"shared")
  fit <- gpc_fit(plan,1,p$patients,p$events)
  expect_error(gpc_occupation(fit),"Aalen")
  expect_error(gpc_fit(plan,1,p$patients,p$events,keep_occupation=TRUE),"keep_occupation")
})
