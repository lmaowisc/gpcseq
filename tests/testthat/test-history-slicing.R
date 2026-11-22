test_that("pair-specific measurement slices preserve table semantics", {
p <- data.frame(id=c('T','C','E'),arm=c(1L,0L,0L),followup=3,death=FALSE)
checked <- 0L
for(value in list(c(2,3,4),c(2L,3L,4L),c(NA_real_,NaN,4),I(c(2,3,4)),setNames(c(2,3,4),letters[1:3])))
for(tt in list(c(0,2,1),c(0L,2L,1L),I(c(0,2,1))))
for(at in c(0,.5,1,1.5,2,3)) for(restricted in c(FALSE,TRUE)) {
 tab <- data.frame(id=c('T','T','C'));tab$time<-tt;tab$value<-value
 empty <- tab[FALSE,]
 seen <- list()
 callback <- function(treatment,control,time){seen[[length(seen)+1L]]<<-list(treatment,control);0}
 provider <- gpcseq:::.pair_tier_provider(p,NULL,1L,2:3,gpc_hierarchy(compare=callback),
   list(events=list(),measurements=list(marker=tab,empty=empty)),restricted=restricted)
 stopifnot(identical(provider(1L,1L,1:2,rep(at,2)),c(0,0)))
 for(k in 1:2)for(side in 1:2){
  id <- if(side==1)'T' else p$id[k+1L]
  rows<-which(tab$id==id)
  expected<-data.frame(time=tab$time[rows],value=tab$value[rows])
  expected<-expected[which(expected$time<=at),,drop=FALSE]
  expect_identical(seen[[k]][[side]]$measurements$marker,expected)
  checked<-checked+1L
 }
}

})
