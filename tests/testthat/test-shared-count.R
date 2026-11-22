test_that("count-only shared windows preserve both patients' chronological projections", {
  p <- data.frame(id=letters[1:8], arm=rep(1:0,4), entry=0,
                  followup=c(2,2,0,3,2,1,3,3),
                  death=c(TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE))
  e <- data.frame(id=c("a","a","b","b","c","e","f","g","h","h"),
                  time=c(0,2,1,1,0,2,1,2,0,3),type="hospitalization")
  count <- gpc_hierarchy(admission_time="none")
  ti <- which(p$arm==1); ci <- which(p$arm==0)
  q <- .pair_summaries(p,e,ti,ci,NA_real_,FALSE,NULL,NULL,count)
  pair <- outer(ti,ci,Vectorize(function(i,j) {
    u <- min(p$followup[c(i,j)])
    died <- p$death[c(i,j)] & p$followup[c(i,j)]==u
    if (any(died)) return(as.integer(died[2L])-as.integer(died[1L]))
    sign(sum(e$id==p$id[j] & e$time<=u)-sum(e$id==p$id[i] & e$time<=u))
  }))
  for (col in 1:2) {
    scored <- if (col==1L) pair>0 else pair<0
    expect_equal(q$mean[col],mean(scored))
    expect_equal(q$average[ti,col],rowMeans(scored))
    expect_equal(q$average[ci,col],colMeans(scored))
  }
  expect_equal(q$n_pairs,length(pair))
  expect_equal(q$n_ascertained,length(pair))
})

test_that("roundoff accepted at a cut does not extend the shared comparison window", {
  p <- data.frame(id=c("T","C"),arm=1:0,entry=0,followup=1,death=FALSE)
  d <- gpc_design(2,.05,"shared",hierarchy=gpc_hierarchy(admission_time="none"))
  at <- function(time) gpc_fit(d,1,p,data.frame(id="T",time=time,type="hospitalization"))$estimate
  expect_equal(at(1),-1)
  expect_equal(at(1+5e-13),0)
  expect_error(at(1+2e-12),"outside observed followup|follow-up|followup")
})
