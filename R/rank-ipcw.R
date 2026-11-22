#' Identify the restricted IPCW hierarchy with an exact rank calculation
#' @param restricted Whether comparisons use a common follow-up horizon.
#' @param biomarker_margin Legacy marker margin, NULL for no marker.
#' @param hierarchy Validated clinical hierarchy or NULL for an older plan.
#' @return TRUE for the built-in survival, count and last-event rule with no
#'   biomarker tier. Old plans and user-specified rules use their own scorer.
#' @noRd
.rank_eligible <- function(restricted, biomarker_margin, hierarchy) {
  # Older saved designs have NULL hierarchy: their lower-tier rule was count
  # only. Explicit ordered lists or user callbacks are never inferred equivalent.
  isTRUE(restricted) && inherits(hierarchy, "gpcseq_hierarchy") &&
    isTRUE(hierarchy$admission_count) && identical(hierarchy$admission_time, "last") &&
    is.null(biomarker_margin) && is.null(hierarchy$biomarker_margin) &&
    is.null(hierarchy$biomarker_threshold) && is.null(hierarchy$compare) &&
    is.null(hierarchy$endpoints)
}

#' Compute restricted pair summaries from risk counts and clinical ranks
#' @param patients,events Validated look-specific records. events is the flat
#'   view of the selected recurrent-event channel (or NULL). Its times and
#'   the follow-up field have a common time-since-entry scale.
#' @param ti,ci Treatment and control row indices in patients.
#' @param tau A supported common follow-up horizon.
#' @param specs Ordered declarative tier specifications after survival, or NULL
#'   for the specialized default count/last calculation.
#' @param channels Named event and measurement tables for declarative tiers.
#'   Preset tiers use the selected native event view.
#' @param km Reverse Kaplan--Meier fits in treatment, control order.
#' @return A list with mean (win, loss), average (one row in original patient
#'   order), n_pairs, n_ascertained, grid (sorted truncated follow-up times),
#'   and mass (weighted win and loss contributions by settlement time).
#'   Its scale and observation convention are those of .pair_summaries.
#'   NULL requests direct scoring when a deciding measurement is missing.
#' @details Inclusive risk counts handle a death and opposite-arm censoring
#'   on the same day. Tied deaths count toward ascertainment. Among patients
#'   observed alive at tau, tied admission counts are resolved by last time.
#'   Grid mass retains the contribution of estimating both censoring curves.
#' @noRd
.rank_ipcw <- function(patients, events, ti, ci, tau, km, specs = NULL, channels = NULL) {
  # One endpoint grid for both arms, as in .pair_summaries(). Follow-up at t
  # is inclusive. A death at tau counts before the clinical horizon ranking.
  x <- patients$followup
  s <- pmin(x, tau)
  grid <- sort(unique(s))
  m <- length(grid)
  ix <- match(s, grid)
  iw <- 1 / (.pair_G(km[[1L]], grid) * .pair_G(km[[2L]], grid))
  died <- patients$death & x <= tau
  nt <- length(ti); nc <- length(ci)
  N <- as.double(nt) * nc
  d1 <- as.double(tabulate(ix[ti[died[ti]]], m))
  d0 <- as.double(tabulate(ix[ci[died[ci]]], m))
  revsum <- function(z) rev(cumsum(rev(as.double(z))))
  r1 <- revsum(tabulate(ix[ti], m))
  r0 <- revsum(tabulate(ix[ci], m))
  win <- d0 * (r1 - d1)
  loss <- d1 * (r0 - d0)
  # Weighted counts through the focal patient's observed end. For two
  # simultaneous deaths the earlier-death contribution must be excluded.
  # Use a strict prefix directly rather than subtracting possibly large
  # nearly equal values from inclusive sums.
  prior1 <- c(0, cumsum(iw * d1))
  prior0 <- c(0, cumsum(iw * d0))
  tset <- ix[ti]; cset <- ix[ci]
  tr_dead <- died[ti]; co_dead <- died[ci]
  raw <- matrix(0.0, nrow(patients), 2L)
  raw[ti, 1L] <- prior0[tset + as.integer(!tr_dead)]
  raw[ti, 2L] <- tr_dead * iw[tset] * (r0[tset] - d0[tset])
  raw[ci, 1L] <- co_dead * iw[cset] * (r1[cset] - d1[cset])
  raw[ci, 2L] <- prior1[cset + as.integer(!co_dead)]
  live1 <- ti[x[ti] >= tau & !tr_dead]
  live0 <- ci[x[ci] >= tau & !co_dead]
  nqual <- sum(win + loss + d0 * d1) + as.double(length(live1)) * length(live0)
  if (length(live1) && length(live0)) {
    # These are the only pairs requiring an admission comparison. Ranking
    # uses exact numeric times; no concatenated string key or time rounding.
    live <- c(live1, live0)
    l1 <- seq_along(live1); l0 <- length(live1) + seq_along(live0)
    if (!is.null(specs)) {
      totals <- .rank_tier_counts(specs, patients$id[live], length(live1), tau, events, channels)
      if (is.null(totals)) return(NULL)
      win1 <- totals$win1; loss1 <- totals$loss1
      loss0 <- totals$loss0; win0 <- totals$win0
    } else {
      nlive <- length(live)
      count <- integer(nlive); last <- numeric(nlive)
      if (!is.null(events) && nrow(events)) {
        use <- events$type == "hospitalization" & events$time <= tau
        k <- match(events$id[use], patients$id[live], nomatch = 0L)
        t <- events$time[use]
        keep <- k > 0L
        if (any(keep)) {
          k <- k[keep]; t <- t[keep]
          count <- tabulate(k, nlive)
          o <- order(k, t, method = "radix")
          kk <- k[o]; end <- !duplicated(kk, fromLast = TRUE)
          last[kk[end]] <- t[o[end]]
        }
      }
      # Worst to best: more admissions first and earlier last event first.
      ord <- order(-count, last, method = "radix")
      group <- integer(nlive)
      group[ord] <- cumsum(c(TRUE, diff(count[ord]) != 0L | diff(last[ord]) != 0))
      ng <- max(group)
      freq1 <- as.double(tabulate(group[l1], ng))
      freq0 <- as.double(tabulate(group[l0], ng))
      p1 <- c(0, cumsum(freq1)); p0 <- c(0, cumsum(freq0))
      g1 <- group[l1]; g0 <- group[l0]
      win1 <- p0[g1]; loss1 <- length(live0) - p0[g1 + 1L]
      loss0 <- p1[g0]; win0 <- length(live1) - p1[g0 + 1L]
    }
    # The eligible horizon pairs all settle in the last grid slot.
    raw[live1, 1L] <- raw[live1, 1L] + iw[m] * win1
    raw[live1, 2L] <- raw[live1, 2L] + iw[m] * loss1
    raw[live0, 1L] <- raw[live0, 1L] + iw[m] * win0
    raw[live0, 2L] <- raw[live0, 2L] + iw[m] * loss0
    win[m] <- win[m] + sum(win1)
    loss[m] <- loss[m] + sum(loss1)
  }
  weighted_mass <- cbind(win, loss) * (iw/N)
  dimnames(weighted_mass) <- NULL
  means <- colSums(raw[ti, , drop = FALSE]) / N
  raw[ti, ] <- raw[ti, , drop = FALSE] / nc
  raw[ci, ] <- raw[ci, , drop = FALSE] / nt
  list(mean=means, average=raw, n_pairs=N, n_ascertained=nqual,
       grid=grid, mass=weighted_mass)
}
