test_that("ordered exact tiers ignore missing values only in resolved groups", {
  p <- data.frame(id=c("t","c"), arm=1:0, entry=0, followup=2, death=FALSE)
  h <- gpc_hierarchy(endpoints=list(first=gpc_endpoint("stroke","first"),
    walk=gpc_endpoint("walk","value",margin=.1)))
  d <- gpc_design(2, .05, "ipcw", horizons=2, hierarchy=h)
  walk <- data.frame(id="c",time=0,value=1)
  stroke <- data.frame(id="t",time=2)
  make <- function(event) gpc_trial_data(p,events=list(stroke=event),measurements=list(walk=walk))
  expect_equal(gpc_fit(d,1,data=make(stroke))$estimate,-1)
  expect_error(gpc_fit(d,1,data=make(stroke[0,])),class="gpc_deferred")
})

test_that("terminal margin respects strict pair differences including rounded boundaries", {
  score <- c(-Inf,-.3,-.2,-.1,0,.1,.2,.3,Inf)
  key <- rep(score,2)
  tr <- seq_along(score); co <- length(tr)+tr
  group <- rep(1L,length(key))
  for (margin in c(.1,.2,1)) {
    result <- .rank_margin_counts(group,key,tr,co,margin)
    gap <- outer(score,score,"-")
    gap[outer(score,score,"==")] <- 0
    expect_equal(result$better,rowSums(gap>margin))
    expect_equal(result$worse,rowSums(gap < -margin))
  }
  positive <- gpc_endpoint("walk","value",margin=.1)
  exact <- gpc_endpoint("relapse","count")
  expect_null(.rank_tiers(gpc_hierarchy(endpoints=list(a=positive,b=exact))))
  expect_length(.rank_tiers(gpc_hierarchy(endpoints=list(a=exact,b=positive))),2)
  expect_null(.rank_tiers(gpc_hierarchy(compare=function(x,y,t) 0)))
})

test_that("general restricted ranks reproduce patient projections and weighted grid mass", {
  p <- data.frame(id=paste0("p",1:8), arm=rep(1:0,each=4),entry=0,
                  followup=c(1,2,3,3,1,2,3,3),death=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE))
  ev <- data.frame(id=p$id[c(3,4,7)],time=c(1,3,2))
  v <- data.frame(id=p$id,time=0,value=c(.2,.3,.5,1.6,1.8,2.5,.3,1.5))
  cut <- gpc_trial_data(p,events=list(relapse=ev),measurements=list(walk=v))
  h <- gpc_hierarchy(endpoints=list(count=gpc_endpoint("relapse","count"),
                    walk=gpc_endpoint("walk","value",direction="higher",margin=.1)))
  km <- lapply(list(1:4,5:8),function(ix) .pair_km(p$followup[ix],p$death[ix],3))
  a <- list(patients=p, events=cut$events, ti=1:4, ci=5:8, tau=3,
            restricted=TRUE, km=km, biomarker_margin=NULL,hierarchy=h,channels=cut$channels)
  expect_equal(do.call(.pair_summaries,a),do.call(.pair_summaries_reference,a),tolerance=1e-13)
})

test_that("callbacks with copied endpoint attributes use the supplied comparison", {
  p <- data.frame(id=c("t","c"), arm=1:0, entry=0, followup=2, death=FALSE)
  cut <- gpc_trial_data(p, events=list(stroke=data.frame(id="c",time=1)))
  original <- gpc_endpoint("stroke", "count")
  retained <- unserialize(serialize(original, NULL))
  expect_identical(.pair_endpoint_spec(retained),attr(original,"gpcseq_endpoint"))
  reversed <- function(treatment, control, time)
    sign(length(treatment$events$stroke)-length(control$events$stroke))
  attr(reversed,"gpcseq_endpoint") <- attr(original,"gpcseq_endpoint")
  changed <- original
  body(changed) <- body(reversed)
  for (f in list(reversed, changed)) {
    expect_null(.pair_endpoint_spec(f))
    for (m in c("ipcw", "shared")) {
      plan <- gpc_design(2,.05,m,horizons=if (m=="ipcw") 2 else NULL,
                         hierarchy=gpc_hierarchy(endpoints=list(stroke=f)))
      expect_equal(gpc_fit(plan,1,data=cut)$estimate,-1)
    }
  }
})

test_that("missing deciding marker reports the original patients and time", {
  p <- data.frame(id=c("alpha","beta"), arm=1:0,entry=0,followup=2,death=FALSE)
  m <- data.frame(id="alpha",time=1,type="biomarker",value=1)
  rule <- gpc_hierarchy(biomarker_margin=1)
  plan <- gpc_design(2,.05,"ipcw",horizons=2,hierarchy=rule)
  expect_error(gpc_fit(plan,1,p,events=m),
    "Missing biomarker when needed for patients alpha and beta at follow-up 2")
})

test_that("first-event timing alone is rankable for either input form", {
  expect_length(.rank_tiers(gpc_hierarchy(admission_count=FALSE,admission_time="first")),1)
  expect_length(.rank_tiers(gpc_hierarchy(endpoints=list(stroke=gpc_endpoint("stroke","first")))),1)
})
