# Slow, deliberately independent case-weight refit oracle. No optimized
# scorer, reverse-KM function, or influence-function helper is called here.
.p_oracle <- function(p, e, tau, method, effect, w, margin = NULL) {
  arm <- lapply(c(1L, 0L), function(a) which(p$arm == a)); D <- p$followup
  D[!p$death] <- Inf
  g <- function(indices, t) {
    j <- indices
    s <- unique(p$followup[j][!p$death[j] & p$followup[j] < t])
    prod(vapply(s, function(z) {
      at_risk <- p$followup[j] > z | (p$followup[j] == z & !p$death[j])
      censored <- p$followup[j] == z & !p$death[j]
      1 - sum(w[j][censored]) / sum(w[j][at_risk])
    }, numeric(1L)))
  }
  probs <- c(0, 0)
  for (i in arm[[1L]]) for (j in arm[[2L]]) {
    end <- min(p$followup[i], p$followup[j])
    if (method == "ipcw" && end < min(D[i], D[j], tau)) next
    tt <- if (method == "shared") end else min(D[i], D[j], tau)
    if (min(D[i], D[j]) <= tt) {
      h <- as.integer(D[i] > D[j]) - as.integer(D[i] < D[j])
    } else {
      hist <- e[e$time <= tt, , drop = FALSE]
      n <- function(a) sum(hist$id == a & hist$type == "hospitalization")
      h <- sign(n(p$id[j]) - n(p$id[i]))
      if (!is.null(margin) && h == 0) {
        bm <- function(a) {
          z <- hist[hist$id == a & hist$type == "biomarker", , drop = FALSE]
          if (!nrow(z)) return(NA_real_)
          z$value[which.max(z$time)]
        }
        gap <- bm(p$id[j]) - bm(p$id[i])
        if (is.na(gap)) stop("Unascertained biomarker in oracle fixture.")
        h <- (gap > margin) - (gap < -margin)
      }
    }
    censor_weight <- if (method == "shared") 1 else g(arm[[1L]], tt) * g(arm[[2L]], tt)
    probs <- probs + w[i] * w[j] * c(h == 1, h == -1) / censor_weight
  }
  probs <- probs / prod(vapply(arm, function(j) sum(w[j]), numeric(1L)))
  net <- probs[1] - probs[2]
  switch(effect, net_benefit = net, win_odds = log1p(net) - log1p(-net),
         win_ratio = log(probs[1]) - log(probs[2]))
}

.p_data <- function() {
  p <- data.frame(id = c(paste0("T", 1:6), paste0("C", 1:6)),
                  arm = rep(c(1L, 0L), each = 6L), entry = 0,
                  followup = c(.8,.8,2.5,4.5,5.5,5, .6,1.3,1.7,3.5,4.5,5),
                  death = c(TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,
                           FALSE,TRUE,FALSE,TRUE,TRUE,FALSE))
  e <- data.frame(id = c("T1","T4","T6","T6","C2","C4","C4","C5","C5","C5"),
                  time = c(.2, 2.9, .7,4, .1, 1,1.1, 1.1,1.8,2.8),
                  type = "hospitalization", value = NA_real_)
  marks <- data.frame(id = c(p$id, "T4","C5","C6"), time = c(rep(0,12), .8, 2, 2),
                      type = "biomarker", value = c(2,1,2,2,1,2, 1,2,3,1,2,0, 0,2,1))
  list(p = p, e = rbind(e, marks))
}

test_that("all pairwise influence curves match case-weight refits of the full estimator", {
  z <- .p_data(); p <- z$p; e <- z$e; n <- nrow(p); base <- rep(1, n)
  for (method in c("shared","ipcw")) {
    effects <- if (method == "shared") "net_benefit" else c("net_benefit", "win_odds","win_ratio")
    for (effect in effects) {
      f <- .pair_fit(p, e, method, effect, 4, biomarker_margin = 0.5)
      est <- function(w) .p_oracle(p, e, 4, method, effect, w, .5)
      expect_equal(f$estimate, est(base), tolerance = 1e-12)
      h <- 1e-5
      fd <- vapply(seq_len(n), function(i) {
        w1 <- w2 <- base; w1[i] <- 1 + h; w2[i] <- 1 - h
        sum(p$arm == p$arm[i]) * (est(w1) - est(w2)) / (2 * h)
      }, numeric(1L))
      expect_equal(f$influence, fd, tolerance = 2e-7)
      expect_equal(sum(f$influence[p$arm == 1]), 0, tolerance = 1e-10)
      expect_equal(sum(f$influence[p$arm == 0]), 0, tolerance = 1e-10)
      # Input row order, cohort order and clinical-event order are incidental.
      ix <- c(12:8, 1:7)
      f2 <- .pair_fit(p[ix, ], e[nrow(e):1, ], method, effect, 4, .5)
      expect_equal(f2$estimate, f$estimate, tolerance = 1e-12)
      expect_equal(f2$influence, f$influence[ix], tolerance = 1e-12)
    }
  }
})

test_that("the first death settles a pair; unsupported pairs retain the full denominator", {
  p <- data.frame(id = c("T1","T2","C1","C2"), arm = c(1,1,0,0),
                  entry = 0, followup = c(2,4,1,4), death = c(FALSE, FALSE, TRUE, FALSE))
  # T1 is censored before its own horizon, but beats C1 because C1 dies first.
  # T1 vs C2 is unascertained; T2 vs C2 is a tie. G(1)=1.
  f <- .pair_fit(p, NULL, "ipcw","net_benefit", 4)
  expect_equal(f$diagnostics$n_ascertained, 3)
  expect_equal(f$natural, 2/4)
  expect_error(.pair_fit(p, NULL, "ipcw","win_ratio", 4), "positive estimated win and loss")
  expect_error(.pair_fit(p, NULL, "ipcw","net_benefit", 4, min_risk = 2L), "[Ii]nsufficient")
})

test_that("a biomarker margin is strict and a missing deciding measurement is not a tie", {
  p <- data.frame(id=c("t","c"), arm=1:0, entry=0, followup=2, death=FALSE)
  e <- data.frame(id=c("t","c"), time=0, type="biomarker", value=c(1,2))
  expect_equal(.pair_fit(p,e,"shared","net_benefit",2,1)$natural, 0)
  expect_equal(.pair_fit(p,e,"shared","net_benefit",2,.9)$natural, 1)
  incomplete <- e[-1, ]
  expect_error(.pair_fit(p,incomplete,"shared","net_benefit",2,1), "[Mm]issing biomarker")
  expect_error(.pair_fit(p,rbind(e,e[1, ]),"ipcw","net_benefit",2,0), "same time")
  p$death[2] <- TRUE
  expect_equal(.pair_fit(p,incomplete,"ipcw","net_benefit",2,1)$natural, 1)
  p$death[1] <- TRUE
  expect_equal(.pair_fit(p,incomplete,"ipcw","net_benefit",2,1)$natural, 0)
  # A biomarker is irrelevant to an IPCW pair if first censored before settlement.
  p$followup <- c(1,2); p$death <- c(FALSE,FALSE)
  p <- rbind(p, data.frame(id="t2", arm=1,entry=0,followup=2,death=FALSE))
  p$death[p$id == "t2"] <- TRUE
  expect_equal(.pair_fit(p,NULL,"ipcw","net_benefit",2,0)$diagnostics$n_ascertained, 1)
})

test_that("reverse KM handles same-date clinical deaths before censorings and excludes the jump at t", {
  km <- .pair_km(c(2,2,3,4), c(TRUE,FALSE,TRUE,TRUE), 4)
  expect_equal(km$risk, 3)
  expect_equal(km$jump, 1/3)
  expect_equal(.pair_G(km, c(2,3)), c(1,2/3))
  # For mass one above s=2, this is minus the relative KM influence.
  v <- .pair_weight_if(c(2,2,3,4), c(TRUE,FALSE,TRUE,TRUE), km,
                       c(2,3), c(1,1,0))
  expect_equal(v, c(0,4/3,-2/3,-2/3))
})

test_that("separate win/loss components yield distinct ratio and odds with clinical ties", {
  z <- .p_data()
  get <- function(method,effect) .pair_fit(z$p,z$e,method,effect,4,.5)
  f <- get("ipcw","net_benefit"); odds <- get("ipcw","win_odds"); ratio <- get("ipcw","win_ratio")
  # Independent exact enumeration of all 36 pairs on the mixed-tie fixture.
  expect_equal(f$diagnostics$win, 1/2)
  expect_equal(f$diagnostics$loss, 4/9)
  expect_equal(f$natural, 1/18)
  expect_equal(odds$natural, 19/17)
  expect_equal(ratio$natural, 9/8)
  expect_equal(odds$influence, 2 / (1 - (1/18)^2) * f$influence)
  expect_equal(get("shared","net_benefit")$natural, 0)
})

test_that("descriptive IPCW effects reuse the clinical fit and agree with single-effect fits", {
  z <- .p_data()
  ds <- gpc_design(6,.05,"ipcw",horizons=4,biomarker_margin=.5)
  old <- serialize(ds,NULL)
  many <- gpc_ipcw_effects(ds,1,z$p,z$e)
  expect_identical(old, serialize(ds,NULL))
  expect_equal(many$horizon, 4)
  expect_equal(many$summary$effect, c("net_benefit","win_odds","win_ratio"))
  expect_equal(many$influences$id, z$p$id)
  expect_true(all(many$summary$valid))
  for (i in 1:3) {
    eff <- many$summary$effect[i]
    single <- gpc_fit(gpc_design(6,.05,"ipcw",effect=eff,horizons=4,
                                 biomarker_margin=.5),1,z$p,z$e)
    phi <- many$influences[[eff]]
    expect_equal(many$summary$estimate[i], single$estimate)
    expect_equal(many$summary$natural[i], single$natural)
    expect_equal(phi,single$subjects$phi)
    expect_equal(many$summary$se[i], sqrt(sum((phi/6)^2)))
  }
})

test_that("undefined descriptive transforms leave valid individual effects available", {
  # Two treatment wins, two clinical ties: win odds is finite; ratio is not.
  p <- data.frame(id=c("t1","t2","c1","c2"),arm=c(1,1,0,0),entry=0,
                  followup=2, death=FALSE)
  e <- data.frame(id=c("t2","c1","c2"),time=.5,type="hospitalization",value=NA_real_)
  plan <- gpc_design(2,.05,"ipcw",horizons=2)
  out <- gpc_ipcw_effects(plan,1,p,e)
  expect_equal(out$summary$valid,c(TRUE,TRUE,FALSE))
  expect_equal(out$summary$natural[1:2],c(.5,3))
  expect_true(all(is.na(out$influences$win_ratio)))
  expect_match(out$summary$reason[3],"positive estimated win and loss")
  # Unanimous wins give a descriptive net benefit of 1 with SE 0; neither
  # log transform exists, even if floating point pair means approach 1 from below.
  unanimous <- gpc_ipcw_effects(plan,1,p,e[e$id != "t2", ])
  expect_equal(unanimous$summary$valid,c(TRUE,FALSE,FALSE))
  expect_equal(unanimous$summary$se[1],0)
  expect_error(.if_ipcw_benefit(list(mean = c(1-.Machine$double.eps,0)),
               NULL,NULL,NULL,TRUE,base_if=0), "working precision")
  expect_error(gpc_ipcw_effects(gpc_design(2,.05,"shared"),1,p,e),"requires a design")
})

test_that("the descriptive extractor runs one score pass, two KM fits, and reuses net influence for odds", {
  # Count calls through private copies with an isolated lexical parent; never
  # unlock or modify a binding in the installed package's namespace.
  ns <- environment(gpc_ipcw_effects)
  env <- new.env(parent = ns)
  called <- new.env(parent = emptyenv())
  for (nm in c("gpc_ipcw_effects", ".pair_core", ".pair_effect",
               ".if_ipcw_benefit", ".if_ipcw_win_ratio", ".pair_fit")) {
    fun <- get(nm, envir = ns); environment(fun) <- env; assign(nm,fun,envir=env)
  }
  tap <- function(nm) {
    force(nm); orig <- get(nm,envir=ns); called[[nm]] <- 0L
    assign(nm, function(...) {
      called[[nm]] <- called[[nm]] + 1L
      orig(...)
    }, envir = env)
  }
  invisible(lapply(c(".pair_km",".pair_summaries",".if_ipcw_linear"), tap))
  z <- .p_data(); design <- gpc_design(6,.05,"ipcw",horizons=4,biomarker_margin=.5)
  ans <- env$gpc_ipcw_effects(design,1,z$p,z$e)
  expect_true(all(ans$summary$valid))
  expect_identical(called$.pair_km,2L)
  expect_identical(called$.pair_summaries,1L)
  expect_identical(called$.if_ipcw_linear,2L)  # signed and ratio bases; odds uses signed
  for (nm in ls(called,all.names=TRUE)) called[[nm]] <- 0L
  only <- env$.pair_fit(z$p,z$e,"ipcw","win_odds",4,.5, hierarchy=design$design$hierarchy)
  expect_equal(only$estimate,ans$summary$estimate[2L])
  expect_identical(called$.pair_km,2L)
  expect_identical(called$.pair_summaries,1L)
  expect_identical(called$.if_ipcw_linear,1L)  # the monitored log odds does not fit the ratio
})

test_that("unanimity under censoring is not turned into finite descriptive log odds by roundoff", {
  # Random stress fixture 9: only T3 supplies ascertained pairs; its four
  # known comparisons are all survival wins. Exact net benefit equals one.
  dat <- data.frame(id=paste0("p",1:9),arm=c(rep(1L,3),rep(0L,6)),entry=0,
                    followup=c(0,0,7,2,5,1,1,4,1),
                    death=as.logical(c(0,0,0,0,1,1,1,1,0)))
  out <- gpc_ipcw_effects(gpc_design(7,.05,"ipcw",horizons=5),1,dat)
  expect_equal(out$summary$natural[1L],1,tolerance=1e-14)
  expect_equal(out$diagnostics$n_ascertained,4)
  expect_equal(out$summary$valid,c(TRUE,FALSE,FALSE))
  expect_match(out$summary$reason[2],"working precision")
  expect_true(all(is.na(out$influences$win_odds)))
})
