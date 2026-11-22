test_that("covariance uses common IDs and the enrolled size in each arm at each look", {
  first <- data.frame(id = c("t1","t2","c1","c2","c3"),
                      arm = c(1L,1L,0L,0L,0L), phi = c(1,-1,1,-1,0))
  second <- data.frame(id = c("t1","t2","t3","c1","c2","c3","c4"),
                       arm = c(1L,1L,1L,0L,0L,0L,0L),
                       phi = c(.5,.5,-1,1,-1,0,0))
  history <- list(list(subjects = first), list(subjects = second))
  oracle <- matrix(c(2/4 + 2/9, 2/12, 2/12, 1.5/9 + 2/16), 2L, 2L)
  expect_equal(.gpc_cov(history), oracle)
  history[[2L]]$subjects <- second[c(5,2,6,1,4,7,3), ]
  expect_equal(.gpc_cov(history), oracle)
  # A previously recorded block is not re-estimated.
  saved <- matrix(oracle[1L,1L] + 1e-10, 1L, 1L)
  update <- .gpc_cov(history, previous = saved)
  expect_identical(update[1L,1L], saved[1L,1L])
  expect_equal(update[2L, ], oracle[2L, ])
  history[[2L]]$subjects$arm[match("t1", history[[2L]]$subjects$id)] <- 0L
  expect_error(.gpc_cov(history), "changed treatment arm")
})

test_that("first usable test uses the scheduled sidedness and cumulative alpha", {
  # A skipped first look has no row in the covariance.
  a <- .gpc_boundary(diag(1L), numeric(), "upper", .025)
  expect_equal(a$critical, qnorm(.975))
  expect_equal(a$probability, .975)
  b <- .gpc_boundary(diag(1L), numeric(), "two-sided", .05)
  expect_equal(b$critical, qnorm(.975))
  expect_equal(b$probability, .95)
})

test_that("independent hybrid rule spends remaining error after a skipped interim", {
  # Three scheduled looks: upper at .01 and .025, two-sided at .05.
  # Stage 2 skipped; probability of reaching the final is .99.
  a <- .gpc_boundary(diag(2L), qnorm(.99), c("upper", "two-sided"), .05)
  expect_equal(a$critical, qnorm(.04 / (2 * .99), lower.tail = FALSE))
  expect_equal(.99 * 2 * pnorm(-a$critical), .04)
  expect_equal(a$probability, .95)
  # If spent error exceeds the putative cumulative target, no valid boundary.
  expect_error(.gpc_boundary(diag(2L), qnorm(.9), rep("upper",2L), .05),
               "spending target")
})

test_that("correlated, upper-only interim and two-sided final use the original regions", {
  corr <- matrix(c(1,.5,.5,1), 2L)
  a <- .gpc_boundary(corr, qnorm(.99), c("upper","two-sided"), .05, tol = 1e-9)
  expect_equal(a$critical, 2.032223105059194, tolerance = 5e-9)
  expect_lt(abs(a$probability - .95), 1e-9)
  # Independent derivation: integrate over the first test, not over rho.
  s <- sqrt(1-.5^2)
  cont <- integrate(function(z) dnorm(z) *
                    (pnorm((a$critical - .5*z)/s) -
                     pnorm((-a$critical - .5*z)/s)),
                    -Inf, qnorm(.99), rel.tol = 1e-10)$value
  expect_equal(cont, .95, tolerance = 5e-9)
  expect_lte(a$error, 1e-9)
})

test_that("multilook correlation is not replaced by canonical information fractions", {
  # Represent Z_j=(X_0-X_j)/sqrt(2) using iid normal X's.
  # All k Z_j <= 0 exactly when X_0 is the smallest; P=1/(k+1).
  for (k in 3:4) {
    R <- matrix(.5, k, k); diag(R) <- 1
    ans <- .gpc_boundary(R, rep(0,k-1L), rep("upper",k), k/(k+1))
    expect_equal(ans$critical, 0, tolerance = 4e-7)
    expect_lte(ans$error, 1e-7)
  }
  # Negative correlation is allowed. With zero upper boundaries the exact
  # three-dimensional orthant probability is known for general correlations.
  R <- matrix(c(1,-.4,.3,-.4,1,.2,.3,.2,1), 3L)
  p <- 1/8 + sum(asin(R[upper.tri(R)]))/(4*pi)
  a <- .gpc_boundary(R, c(0,0), rep("upper",3L), 1-p)
  expect_equal(a$critical, 0, tolerance = 4e-7)
})

test_that("singular correlations and invalid look directions do not yield a boundary", {
  expect_error(.gpc_boundary(matrix(1,2L,2L), 2, c("upper","upper"), .05),
               "singular")
  expect_error(.gpc_boundary(diag(2L), 2, c("sideways","upper"), .05),
               "same looks")
})

test_that("rounded duplicate and linearly dependent looks are consistently deferred", {
  # Repeated patient influences at looks two and three in a staggered trial.
  # The off-diagonal correlation rounded below 1 when covariance was scaled.
  r <- 0.71980084072973427
  d <- 0.99999999999999967
  repeated <- matrix(c(1,r,r, r,1,d, r,d,1), 3L)
  # No two components coincide when Z3 = (Z1 - Z2)/sqrt(2).
  s <- sqrt(.5)
  dependent <- matrix(c(1,0,s, 0,1,-s, s,-s,1), 3L)
  for (R in list(repeated, dependent)) {
    # Monitoring effects may have very different test scales.
    scales <- c(1e-4, 2, 1e3)
    scaled <- cov2cor(R * outer(scales,scales))
    for (idx in list(1:3, c(3,1,2), c(2,3,1))) {
      expect_error(.gpc_boundary(scaled[idx,idx], c(3,2.5), rep("upper",3L), .05),
                   "numerically singular", class="gpc_deferred")
    }
  }
})

test_that("almost matching looks with resolvable variation still give a boundary", {
  # At this separation the probability of crossing the upper interim
  # boundary but missing the looser final upper boundary is negligible.
  rho <- 1 - 1e-8
  R <- matrix(c(1,rho,rho,1),2L)
  fit <- .gpc_boundary(R, qnorm(.99), c("upper","upper"), .05)
  expect_equal(fit$critical, qnorm(.95), tolerance=1e-6)
  expect_equal(fit$probability, .95, tolerance=1e-7)
})

test_that("high-dimensional calibration preserves the caller's random stream", {
  skip_if_not_installed("mvtnorm")
  k <- 5L; R <- matrix(.5, k, k); diag(R) <- 1
  set.seed(123)
  before <- .Random.seed
  a <- .gpc_boundary(R, rep(0,k-1L), rep("upper",k), k/(k+1), tol = 2e-4)
  expect_identical(.Random.seed, before)
  expect_lt(abs(a$critical), .005)
  expect_lte(a$error, 2e-4)
})

test_that("successive four-look cutoffs match an independent one-factor oracle", {
  # Equicorrelation 1/2: Z_j=(U+eps_j)/sqrt(2), independent U, eps_j.
  # Conditional on U only, all look rectangles factor. The reference constants
  # were computed independently by scalar integration over U, not by the
  # product quadrature used in the package.
  alpha <- c(.002,.010,.025,.05)
  expected <- list(c(3.090232306168,2.640661145539,2.397165357172,2.176478542609),
                   c(2.878161739095,2.392956806716,2.120019704718,2.189975180518))
  patterns <- list(rep("two-sided",4L), c(rep("upper",3L),"two-sided"))
  for (i in seq_along(patterns)) {
    cut <- numeric()
    for (k in 1:4) {
      R <- matrix(.5,k,k); diag(R) <- 1
      d <- patterns[[i]][seq_len(k)]
      ans <- .gpc_boundary(R, cut, d, alpha[k])
      cut <- c(cut, ans$critical)
      expect_lt(abs(ans$critical - expected[[i]][k]), 2e-6)
      # Recheck probability at the produced, rounded-to-machine boundary;
      # an independent integral also tests the joint accumulation of earlier errors.
      reference <- integrate(function(u) {
        p <- rep(1,length(u))
        for (j in seq_len(k)) {
          upper <- pnorm(sqrt(2)*cut[j]-u)
          lower <- if (d[j] == "two-sided") pnorm(-sqrt(2)*cut[j]-u) else 0
          p <- p * (upper-lower)
        }
        dnorm(u)*p
      }, -9, 9, rel.tol = 1e-11)$value
      expect_lt(abs(reference - (1-alpha[k])), 1e-7)
    }
  }
})

test_that("independent looks allocate small crossing probabilities on the tail scale", {
  for (dirs in list(c("upper", "upper"), c("upper", "two-sided"),
                    c("two-sided", "two-sided"))) {
    p <- 1e-300
    sides <- ifelse(dirs == "two-sided", 2, 1)
    old <- qnorm(p / sides[1], lower.tail = FALSE)
    ans <- .gpc_boundary(diag(2), old, dirs, 10*p, tol = p/32)
    added <- sides[2] * pnorm(ans$critical, lower.tail = FALSE)
    expect_lt(abs(added/(9*p) - 1), 1e-12)
  }
  # A few ulps of the previously spent probability cannot be certified
  # to this stricter tolerance even if a root numerically rounds to its target.
  expect_error(.gpc_boundary(diag(2), qnorm(.975), rep("two-sided",2),
                             .05 + .Machine$double.eps, tol=1e-17),
               class="gpc_deferred")
})

test_that("nearly repeated or sign-reversed looks preserve cumulative spending", {
  for (k in 3:5) {
    old <- c(2.326348, 2.15, 2.07, 1.99)[seq_len(k-1)]
    for (sgn in list(rep(1, k), rep(c(1,-1),length.out=k))) {
      R <- (1 - 1e-8) * tcrossprod(sgn)
      diag(R) <- 1
      dirs <- c(rep("upper",k-1), "two-sided")
      # Spacing between competing same-tail cutoffs exceeds eight SDs of
      # the difference. The active bounds agree with the common-factor limit
      # well within the requested absolute probability tolerance.
      upper <- min(old[sgn[-k] == 1])
      left <- old[sgn[-k] == -1]
      lower <- if (length(left)) -min(left) else -Inf
      alpha <- pnorm(upper, lower.tail=FALSE) + pnorm(lower) + .012
      ans <- .gpc_boundary(R, old, dirs, alpha)
      actual <- pnorm(min(upper, ans$critical),lower.tail=FALSE) +
        pnorm(max(-ans$critical, lower))
      expect_lt(abs(actual - alpha), 1e-7)
      expect_lte(ans$error, 1e-7)
    }
  }
  # Four nearly identical, equally restrictive constraints remain distinct.
  R <- matrix(1 - 1e-8, 5, 5); diag(R) <- 1
  expect_error(.gpc_boundary(R, rep(2,4), rep("upper",5), .05),
               "could not be resolved", class="gpc_deferred")
})

test_that("near linear dependence without duplicate looks uses the Gaussian orthant law", {
  s <- sqrt(.5 - 1e-9)
  R <- matrix(c(1,0,s, 0,1,-s, s,-s,1), 3L)
  # P(all <= 0) = 1/8 + sum(asin(rho_ij))/(4*pi) = 1/8.
  set.seed(135)
  before <- .Random.seed
  ans <- .gpc_boundary(R, c(0,0), rep("upper",3), .875)
  expect_lt(abs(ans$critical), 1e-6)
  expect_identical(.Random.seed, before)
})


test_that("root evaluation reuses exact cutoffs and preserves integration error", {
  seen <- numeric()
  evaluate <- function(cut) {seen <<- c(seen,cut); c(p=pnorm(cut),e=2e-9)}
  z <- .gpc_boundary_root(evaluate,.975,FALSE,1e-7)
  expect_equal(z$critical,qnorm(.975),tolerance=1e-7)
  expect_identical(z$integration_error,2e-9)
  expect_equal(z$error,abs(z$probability-.975)+2e-9)
  expect_false(anyDuplicated(seen)>0)
  other <- .gpc_boundary_root(function(cut)c(p=pnorm(cut-1),e=3e-9),.975,FALSE,1e-7)
  expect_equal(other$critical,qnorm(.975)+1,tolerance=1e-7)
  expect_identical(other$integration_error,3e-9)
})
