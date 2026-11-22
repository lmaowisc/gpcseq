test_that("custom histories retain data-frame construction semantics", {
  p <- data.frame(id=c("T","C","E"),arm=c(1L,0L,0L),followup=3,death=FALSE)
  empty <- data.frame(id=character(),time=numeric(),value=numeric())
  for (value in list(c(2,3,4),c(2L,3L,4L),c(NA_real_,NaN,4),I(c(2,3,4)))) {
    tab <- data.frame(id=c("T","T","C"),time=c(0,2,1))
    tab$value <- value
    seen <- list()
    callback <- function(treatment,control,time) {
      seen[[length(seen)+1L]] <<- list(treatment,control)
      0
    }
    channels <- list(events=list(),measurements=list(marker=tab,empty=empty))
    provider <- gpcseq:::.pair_tier_provider(p,NULL,1L,2:3,
      gpc_hierarchy(compare=callback),channels,restricted=TRUE)
    expect_identical(provider(1L,1L,1:2,c(3,3)),c(0,0))
    for (k in 1:2) for (side in 1:2) {
      id <- if(side==1) "T" else p$id[k+1L]
      rows <- which(tab$id==id)
      expected <- data.frame(time=tab$time[rows],value=tab$value[rows])
      expect_identical(seen[[k]][[side]]$measurements$marker,expected)
      expect_identical(seen[[k]][[side]]$measurements$empty,empty[c("time","value")])
    }
  }
})
