# All helpers below are private. Package documentation in inst/technical/boundaries.md
# gives their contracts, integration tolerances, failure policy and reference examples.

# history: completed fits with one (id, arm, phi) row per patient; previous: optional
# (K-1)-square saved covariance. Returns a K-square covariance of estimates.
# Same-patient covariance. Divide each influence by its own arm's enrolled
# count *at that look*. An absent ID contributes zero; never re-fit old data.
.gpc_cov <- function(history, previous = NULL) {
  k <- length(history)
  if (!k) return(matrix(numeric(), 0L, 0L))
  rows <- lapply(history, function(look) {
    x <- look$subjects
    if (!is.data.frame(x) || !all(c("id", "arm", "phi") %in% names(x)) ||
        !is.character(x$id) || anyNA(x$id) || any(!nzchar(x$id)) ||
        anyDuplicated(x$id) || !is.numeric(x$arm) || anyNA(x$arm) || !all(x$arm %in% 0:1) ||
        !is.numeric(x$phi) || any(!is.finite(x$phi)))
      stop("Each tested look needs unique patient IDs, arms 0/1 and finite influences.",
           call. = FALSE)
    n <- tabulate(as.integer(x$arm) + 1L, nbins = 2L)
    if (any(n == 0L)) .gpc_defer("No enrolled patient in one arm.")
    list(id = x$id, arm = x$arm, score = x$phi / n[x$arm + 1L])
  })
  inner <- function(left, right) {
    hit <- match(left$id, right$id)
    both <- !is.na(hit)
    if (any(left$arm[both] != right$arm[hit[both]]))
      stop("A patient ID changed treatment arm between looks.", call. = FALSE)
    sum(left$score[both] * right$score[hit[both]])
  }
  V <- matrix(0, k, k)
  if (!is.null(previous)) {
    if (!is.numeric(previous) || !is.matrix(previous) ||
        !all(dim(previous) == k - 1L) || any(!is.finite(previous)) ||
        !isTRUE(isSymmetric(unname(previous))))
      stop("The saved covariance block must match the earlier completed looks.", call. = FALSE)
    earlier <- seq_len(k - 1L)
    V[earlier, earlier] <- previous
    start <- k
  } else start <- 1L
  for (j in seq.int(start, k)) {
    # Existing rows are read only. A new look adds one column and its transpose.
    d <- vapply(seq_len(j - 1L), function(i) inner(rows[[i]], rows[[j]]), numeric(1))
    V[seq_len(j - 1L), j] <- d
    V[j, seq_len(j - 1L)] <- d
    V[j, j] <- sum(rows[[j]]$score ^ 2)
  }
  V
}

# Equally sized numeric vectors a, b (possibly infinite); one probability per interval.
# P(a < N(0,1) <= b); retain right-tail accuracy when both limits are positive.
.gpc_normal_interval <- function(a, b) {
  right <- a > 0
  p <- stats::pnorm(b) - stats::pnorm(a)
  if (any(right)) p[right] <- stats::pnorm(a[right], lower.tail = FALSE) -
    stats::pnorm(b[right], lower.tail = FALSE)
  p
}

# lo/hi: two limits each; rho in (-1,1); tol: allowed absolute probability error.
# Returns named numeric p/e (continuation probability / quadrature error estimate).
# Bivariate normal rectangle: integrate the derivative with respect to the
# correlation angle on a finite interval (Plackett identity). Infinite limits
# have zero derivative. No random integration or subtraction of four large CDFs.
.gpc_rect2 <- function(lo, hi, rho, tol) {
  base <- prod(.gpc_normal_interval(lo, hi))
  if (rho == 0) return(c(p = base, e = 0))
  kernel <- function(x, y, th) {
    if (!is.finite(x) || !is.finite(y)) return(numeric(length(th)))
    exp(-0.5 * (y * y + ((x - sin(th) * y) / cos(th)) ^ 2))
  }
  f <- function(th) kernel(hi[1L], hi[2L], th) - kernel(lo[1L], hi[2L], th) -
    kernel(hi[1L], lo[2L], th) + kernel(lo[1L], lo[2L], th)
  angle <- asin(rho)
  z <- tryCatch(stats::integrate(f, min(0, angle), max(0, angle),
                                  abs.tol = tol * pi / 4, rel.tol = 0,
                                  subdivisions = 250L), error = function(e) NULL)
  if (is.null(z)) .gpc_defer("Bivariate normal integration did not converge.")
  c(p = base + sign(angle) * z$value / (2 * pi), e = z$abs.error / (2 * pi))
}

# Positive integer n; returns list(x, w) of nodes and weights, of length n.
# Fixed Gauss-Legendre nodes from the symmetric Jacobi matrix, on [-1,1].
.gpc_legendre <- local({
  cache <- new.env(parent = emptyenv())
  function(n) {
    key <- as.character(n)
    if (!is.null(cache[[key]])) return(cache[[key]])
    i <- seq_len(n - 1L)
    j <- matrix(0, n, n)
    j[cbind(i, i + 1L)] <- i / sqrt(4 * i ^ 2 - 1)
    ev <- eigen(j + t(j), symmetric = TRUE)
    out <- list(x = ev$values, w = 2 * ev$vectors[1L, ] ^ 2)
    cache[[key]] <- out
    out
  }
})

# L: lower Cholesky factor of the current K-square correlation, lo/hi: K-1 old
# continuation limits. ngrid: rule order; tol: absolute probability tolerance.
# Returns weights, conditional means and common SD of the last test.
# Conditional normal quadrature for 3-4 looks: integrate earlier continuation
# limits once per grid, then vary only the final normal limits during bisection.
# The integration dimensions are independent standard normal innovations.
.gpc_prefix <- function(L, lo, hi, ngrid, tol) {
  k <- nrow(L)
  z <- .gpc_legendre(ngrid)
  tail <- stats::qnorm(tol / (64 * k), lower.tail = FALSE)
  mass <- 1
  mu <- matrix(0, 1L, k)
  for (j in seq_len(k - 1L)) {
    a <- pmax(-tail, (lo[j] - mu[, 1L]) / L[j, j])
    b <- pmin(tail, (hi[j] - mu[, 1L]) / L[j, j])
    good <- b > a
    if (!any(good)) return(list(weight = 0, mean = 0, sd = L[k, k]))
    a <- a[good]; b <- b[good]; mass <- mass[good]; mu <- mu[good, , drop = FALSE]
    half <- (b - a) / 2
    ii <- rep(seq_along(a), each = ngrid)
    draws <- ((b + a) / 2)[ii] + half[ii] * rep(z$x, length(a))
    mass <- (mass * half)[ii] * rep(z$w, length(a)) * stats::dnorm(draws)
    # Drop the consumed coordinate. At the largest grid only one column
    # remains, rather than materializing all K conditional means.
    future <- seq.int(j + 1L, k)
    mu <- mu[ii, -1L, drop = FALSE] + tcrossprod(draws, L[future, j])
  }
  list(weight = mass, mean = mu[, 1L], sd = L[k, k])
}

# prefix: output of .gpc_prefix; lo/hi: candidate lower and upper Z for the last
# test. Returns one unconditional continuation probability (summing prefix mass).
.gpc_prefix_eval <- function(prefix, lo, hi) {
  sum(prefix$weight * .gpc_normal_interval((lo - prefix$mean) / prefix$sd,
                                            (hi - prefix$mean) / prefix$sd))
}

# lo/hi: K-vectors, corr: K-square positive-definite correlation, tol: absolute
# probability accuracy. Returns named numeric p/e, or defers if no finite result.
# Use the same random integration points for every candidate in a root search.
# Restore the caller's RNG state even with mvtnorm versions predating seed=.
.gpc_rect_qmc <- function(lo, hi, corr, tol) {
  if (any(lo == hi)) return(c(p = 0, e = 0))
  old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", .GlobalEnv)
  kind <- RNGkind()
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(kind)))
    if (is.null(old)) {
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(120819L, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  val <- tryCatch(mvtnorm::pmvnorm(lower = lo, upper = hi, corr = corr,
    algorithm = mvtnorm::GenzBretz(maxpts = 2000000L, abseps = tol / 4,
                                    releps = 0)),
    error = function(e) NULL)
  if (is.null(val) || !is.finite(val) || length(attr(val, "error")) != 1L ||
      !is.finite(attr(val, "error")))
    .gpc_defer("Multivariate normal integration failed at this look.")
  c(p = as.numeric(val), e = as.numeric(attr(val, "error")))
}

# In a nearly dependent rectangle the random-integration error estimate can be
# too small when a narrow transition lies between integration points. If two
# signed constraints are well separated, bound the probability that the weaker
# side binds while the stronger side holds using their Gaussian difference.
.gpc_dependent_rect <- function(lo, hi, R, tol) {
  if (any(lo == hi)) return(c(p = 0, e = 0))
  k <- length(hi)
  omission <- 0
  if (k > 3L) {
    finite_lower <- which(is.finite(lo))
    index <- c(seq_len(k), finite_lower)
    sign <- c(rep(1, k), rep(-1, length(finite_lower)))
    limit <- c(hi, -lo[finite_lower])
    kept <- integer()
    for (face in order(limit)) {
      rho <- sign[face] * sign[kept] * R[index[face], index[kept]]
      sd <- sqrt(2 * pmax(0, 1 - rho))
      risk <- stats::pnorm((limit[face] - limit[kept]) / sd, lower.tail = FALSE)
      if (length(risk) && min(risk) <= tol / (16 * length(limit))) {
        omission <- omission + min(risk)
      } else kept <- c(kept, face)
    }
    lower <- rep(-Inf, k)
    upper <- rep(Inf, k)
    for (face in kept) {
      if (sign[face] > 0) upper[index[face]] <- limit[face]
      else lower[index[face]] <- -limit[face]
    }
    use <- is.finite(lower) | is.finite(upper)
    lo <- lower[use]; hi <- upper[use]; R <- R[use, use, drop = FALSE]
    k <- length(hi)
  }
  if (k == 1L) val <- c(p = .gpc_normal_interval(lo, hi), e = 0)
  else if (k == 2L) val <- .gpc_rect2(lo, hi, R[1L, 2L], tol / 2)
  else if (k == 3L) {
    # TVPACK evaluates one lower orthant at a time. Flip coordinates that
    # have only a lower limit so that every upper limit is finite.
    flip <- !is.finite(hi)
    if (any(flip)) {
      hi[flip] <- -lo[flip]
      lo[flip] <- -Inf
      signs <- ifelse(flip, -1, 1)
      R <- R * outer(signs, signs)
    }
    # Some mvtnorm versions initialize the caller's RNG even for TVPACK.
    if (!exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      on.exit(if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    varying <- which(is.finite(lo))
    corners <- 2 ^ length(varying)
    eps <- tol / (4 * corners)
    p <- 0
    for (mask in seq_len(corners) - 1L) {
      use <- varying[as.logical(intToBits(mask))[seq_along(varying)]]
      upper <- hi
      upper[use] <- lo[use]
      term <- mvtnorm::pmvnorm(upper = upper, corr = R,
                               algorithm = mvtnorm::TVPACK(abseps = eps))
      p <- p + (-1) ^ length(use) * as.numeric(term)
    }
    # The native TVPACK code floors each trivariate CDF request at 1e-14 and
    # reports the requested epsilon rather than the achieved error. Budget
    # the larger of the request and that floor for every corner, so an
    # unattainably tight request defers instead of appearing to succeed.
    # This does not certify convergence at the backend's subdivision cap.
    val <- c(p = p, e = corners * (max(eps, 1e-14) + 4 * .Machine$double.eps))
  } else {
    spectrum <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    if (min(spectrum) < 1e-3 * max(spectrum))
      .gpc_defer("The nearly dependent look boundaries could not be resolved at this tolerance.")
    val <- .gpc_rect_qmc(lo, hi, R, tol / 2)
  }
  val["e"] <- val["e"] + omission
  val
}

# eval(cut): named numeric p/e; target: unconditional continuation; twosided:
# logical; tol: absolute probability tolerance. Returns fitted cut, p and errors.
# Bracket an increasing continuation probability. A benefit-only look can
# legitimately have a negative critical value for a sufficiently large alpha.
.gpc_boundary_root <- function(eval, target, twosided, tol) {
  # Evaluators use fixed integration points within a root search. Retain the
  # full probability/error pair at exact repeated cutoffs; never round keys.
  points <- numeric(); values <- list()
  evaluate <- function(cut) {
    i <- match(cut, points, nomatch = 0L)
    if (i) return(values[[i]])
    value <- eval(cut)
    points <<- c(points, cut)
    values[[length(points)]] <<- value
    value
  }
  lower <- if (twosided) 0 else -1
  upper <- 2
  f <- function(cut) unname(evaluate(cut)["p"]) - target
  fl <- f(lower); fu <- f(upper)
  while (!twosided && fl > 0 && lower > -40) { lower <- lower * 2; fl <- f(lower) }
  while (fu < 0 && upper < 40) { upper <- upper * 2; fu <- f(upper) }
  if (!is.finite(fl) || !is.finite(fu) || fl > 0 || fu < 0)
    .gpc_defer("Earlier boundaries leave insufficient probability for this spending target.")
  ans <- tryCatch(stats::uniroot(f, c(lower, upper), f.lower = fl, f.upper = fu,
                        tol = max(.Machine$double.eps, tol / 4), check.conv = TRUE),
                  error = function(e) .gpc_defer(paste("Boundary root search failed:",
                                                       conditionMessage(e))))
  cut <- ans$root
  val <- evaluate(cut)
  list(critical = cut, probability = unname(val["p"]),
       integration_error = unname(val["e"]),
       error = abs(unname(val["p"]) - target) + unname(val["e"]))
}

# correlation: K-square correlation; previous: K-1 fixed critical values;
# directions: K test directions; alpha: cumulative spend; tol: absolute probability
# error budget. Returns the critical value, continuation probability, error and method.
# Probability and error are on the *continuation* scale. Previous contains
# only performed tests, each with its original critical value/direction;
# alpha is the cumulative spending at this scheduled look (not a conditional
# level). A skipped interim therefore needs no special numerical branch.
.gpc_boundary <- function(correlation, previous, directions, alpha, tol = 1e-7) {
  k <- length(directions)
  if (!k || !is.numeric(correlation) || !is.matrix(correlation) ||
      !all(dim(correlation) == k) || length(previous) != k - 1L ||
      !is.numeric(previous) || any(!is.finite(previous)) ||
      anyNA(directions) || !all(directions %in% c("upper", "two-sided")) ||
      any(previous[directions[seq_along(previous)] == "two-sided"] < 0))
    stop("The correlation, directions and earlier critical values must describe the same looks.",
         call. = FALSE)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= 1)
    stop("alpha and the probability tolerance must lie strictly between zero and one.",
         call. = FALSE)
  R <- unname(correlation)
  if (any(!is.finite(R)) || any(abs(diag(R) - 1) > 1e-12) ||
      any(abs(R) > 1 + 1e-12) || max(abs(R - t(R))) > 1e-12)
    .gpc_defer("The look correlation is not a finite symmetric correlation matrix.")
  # Correlation standardization can give a tiny positive eigenvalue even
  # for identical patient influences. Test numerical rank on the correlation
  # scale: individual Cholesky pivots depend on the ordering of the looks.
  if (k > 1L) {
    spectrum <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    if (min(spectrum) <= k * .Machine$double.eps * max(spectrum))
      .gpc_defer("The look correlation matrix is numerically singular.")
  }
  L <- tryCatch(t(chol(R)), error = function(e) NULL)
  if (is.null(L))
    .gpc_defer("The look correlation matrix is numerically singular.")
  two <- directions[k] == "two-sided"
  if (k == 1L) {
    # No earlier test conditions the first boundary. Work on the tail scale
    # so a positive scheduled crossing probability can be smaller than 1 ulp at 1.
    log_sides <- if (two) log(2) else 0
    log_tail <- log(alpha) - log_sides
    cut <- stats::qnorm(log_tail, lower.tail = FALSE, log.p = TRUE)
    fitted_log_tail <- stats::pnorm(cut, lower.tail = FALSE, log.p = TRUE)
    err <- alpha * abs(expm1(fitted_log_tail - log_tail))
    if (!is.finite(cut) || !is.finite(err) || err > min(tol, alpha / 8))
      .gpc_defer("The normal calibration is below numerical precision.")
    return(list(critical = cut, probability = 1 - exp(fitted_log_tail + log_sides),
                integration_error = 0, error = err, method = "normal"))
  }
  lo <- ifelse(directions[seq_along(previous)] == "two-sided", -previous, -Inf)
  tol <- min(tol, alpha / 8)
  if (all(R[upper.tri(R)] == 0)) {
    # Accumulate old continuation in logs. Differences of probabilities near
    # one can otherwise round positive spending increments to an unrelated ulp.
    log_reach <- stats::pnorm(previous, log.p = TRUE)
    symmetric <- is.finite(lo)
    log_reach[symmetric] <- log1p(-2 * stats::pnorm(previous[symmetric], lower.tail = FALSE))
    log_reach <- sum(log_reach)
    inc <- alpha + expm1(log_reach)
    if (inc <= 0) .gpc_defer("The spending target must increase beyond previous spending.")
    log_tail <- log(inc) - log_reach - if (two) log(2) else 0
    cut <- stats::qnorm(log_tail, lower.tail = FALSE, log.p = TRUE)
    fitted_log_tail <- stats::pnorm(cut, lower.tail = FALSE, log.p = TRUE)
    difference <- inc * expm1(fitted_log_tail - log_tail)
    err <- abs(difference) + 8 * .Machine$double.eps * alpha
    if (!is.finite(cut) || !is.finite(err) || err > min(tol, inc / 8))
      .gpc_defer("The normal calibration is below numerical precision.")
    return(list(critical = cut, probability = 1 - (alpha + difference),
                integration_error = 0, error = err, method = "normal"))
  }
  target <- 1 - alpha
  if (tol < 8 * .Machine$double.eps)
    .gpc_defer("The joint continuation probability cannot resolve the requested tolerance.")
  if (k >= 3L && min(spectrum) < 1e-3 * max(spectrum)) {
    prior <- .gpc_dependent_rect(lo, previous, R[-k, -k, drop = FALSE], tol)
    room <- prior["p"] - target - prior["e"]
    if (room <= 0) .gpc_defer("A positive new spending increment could not be resolved.")
    tol <- min(tol, unname(room) / 8)
    if (tol < 8 * .Machine$double.eps)
      .gpc_defer("The joint continuation probability cannot resolve the spending increment.")
    eval <- function(cut) .gpc_dependent_rect(c(lo, if (two) -cut else -Inf),
                                              c(previous, cut), R, tol)
    ans <- .gpc_boundary_root(eval, target, two, tol)
    if (!is.finite(ans$error) || ans$error > tol)
      .gpc_defer("The nearly dependent look boundaries did not meet the requested tolerance.")
    ans$method <- "reduced-rectangle"
    return(ans)
  }
  if (k == 2L) {
    pold <- .gpc_normal_interval(lo, previous)
    if (pold <= target) .gpc_defer("The spending target must increase beyond previous spending.")
    tol <- min(tol, (pold - target) / 8)
    if (tol < 8 * .Machine$double.eps)
      .gpc_defer("The joint continuation probability cannot resolve the spending increment.")
    eval <- function(cut) .gpc_rect2(c(lo, if (two) -cut else -Inf),
                                     c(previous, cut), R[1L, 2L], tol)
    ans <- .gpc_boundary_root(eval, target, two, tol); ans$method <- "bivariate"
    if (ans$error <= tol) return(ans)
    .gpc_defer("The bivariate boundary did not meet its probability tolerance.")
  }
  if (k <= 4L) {
    coarse <- .gpc_prefix(L, lo, previous, 12L, tol)
    for (n in c(24L, 48L, 80L)) {
      fine <- .gpc_prefix(L, lo, previous, n, tol)
      pold <- sum(fine$weight)
      room <- pold - target
      if (room <= 0) { coarse <- fine; next }
      allowed <- min(tol, room / 8)
      eval <- function(cut) c(p = .gpc_prefix_eval(fine, if (two) -cut else -Inf,
                                                    cut), e = 0)
      ans <- .gpc_boundary_root(eval, target, two, allowed)
      pf <- .gpc_prefix_eval(coarse, if (two) -ans$critical else -Inf, ans$critical)
      # Grid disagreement estimates quadrature error; the discarded normal
      # tails contribute at most tol/32 by the union bound.
      err <- max(abs(ans$probability - pf),
                 abs(pold - sum(coarse$weight))) + tol / 32
      ans$integration_error <- err; ans$error <- err + abs(ans$probability - target)
      ans$method <- paste0("conditional-", n)
      if (ans$error <= allowed) return(ans)
      coarse <- fine
    }
  }
  # More than four coordinates or unresolved small-dimensional quadrature.
  # A failed stochastic integration is reported rather than replacing R.
  prior <- R[-k, -k, drop = FALSE]
  prev <- .gpc_rect_qmc(lo, previous, prior, tol)
  room <- prev["p"] - target - prev["e"]
  if (room <= 0) .gpc_defer("A positive new spending increment could not be resolved.")
  tol <- min(tol, unname(room) / 8)
  eval <- function(cut) .gpc_rect_qmc(c(lo, if (two) -cut else -Inf),
                                      c(previous, cut), R, tol)
  ans <- .gpc_boundary_root(eval, target, two, tol); ans$method <- "Genz-Bretz"
  if (!is.finite(ans$error) || ans$error > tol)
    .gpc_defer("Multivariate normal boundary integration did not meet its probability tolerance.")
  ans
}
