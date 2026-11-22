# Reverse KM, projected pair scores and weight influences are separate so a
# ratio/odds transformation never has to refit censoring or rescore patients.

#' Fit the arm-specific reverse Kaplan--Meier jumps at one look
#'
#' @param x Nonnegative follow-up times of the enrolled patients in one arm,
#'   on the same scale as the horizon. Each is the observed time to death,
#'   censoring, or the calendar cutoff, whichever occurred first.
#' @param died Logical vector aligned with x: TRUE is a documented death;
#'   FALSE is an observed censoring time, including the calendar cutoff.
#' @param tau Nonnegative supported horizon. Only jumps strictly before tau
#'   enter P(C >= t) for t <= tau.
#' @return List with sorted censoring jump times (time), death-adjusted
#'   risk counts (risk), censoring hazard increments (jump), left endpoint
#'   survival 1 followed by post-jump products (G), and enrollment (n).
#'   When there are no jumps before tau, G contains only 1.
#' @details Tied documented deaths leave the censoring risk set before that
#'   date's censorings. The risk denominator is the number whose follow-up
#'   exceeds the jump time plus those censored exactly then. This helper fits
#'   the empirical curve even if it reaches zero; .pair_fit checks support
#'   before the curve is inverted. O(n log n) time and O(n) storage.
#' @keywords internal
#' @noRd
.pair_km <- function(x, died, tau) {
  n <- length(x)
  is_cens <- !died & x < tau
  s <- sort(unique(x[is_cens]))
  if (!length(s)) return(list(time = numeric(), risk = numeric(),
                              jump = numeric(), G = 1, n = n))
  nc <- tabulate(match(x[is_cens], s), length(s))
  # Documented clinical death precedes censoring on a shared date.
  risk <- n - findInterval(s, sort(x), left.open = TRUE) -
    tabulate(match(x[died], s, nomatch = 0L), length(s))
  h <- nc / risk
  list(time = s, risk = risk, jump = h, G = c(1, cumprod(1 - h)), n = n)
}

#' Evaluate the censoring survival immediately before each requested time
#'
#' @param km Result of .pair_km for the arm and current look.
#' @param t Vector of follow-up times from zero through the horizon used
#'   to fit km, in arbitrary order.
#' @return Numeric vector P-hat(C >= t) in the order of t. Returns 1
#'   before and at the first censoring jump, including a jump at zero.
#' @details Uses the last jump strictly earlier than t. No curve is refit;
#'   this lookup assumes the caller already selected the appropriate arm,
#'   look and supported horizon.
#' @keywords internal
#' @noRd
.pair_G <- function(km, t) {
  if (!length(km$time)) return(rep(1, length(t)))
  km$G[1L + findInterval(t, km$time, left.open = TRUE)]
}

#' Index the clinical event histories for a one-pass treatment--control comparison
#'
#' @param patients Validated patient snapshot with unique id, numeric arm,
#'   followup and death status. Records must contain only observations
#'   available through the patient cutoff at this look.
#' @param events Long-form clinical records: id, time since entry, type, and
#'   value for a biomarker. Each hospitalization row counts as one admission.
#'   NULL or no rows represents no recorded events in a complete history.
#' @param ti,ci Patient row indices in treatment and control, respectively.
#' @param biomarker Whether to index the last biomarker tier. When FALSE,
#'   biomarker data are not inspected.
#' @return th, tb and tv are lists indexed in ti order, giving sorted
#'   admission times, biomarker times and marker values. ch/chi and cb/cbi/cv
#'   give the sorted event times, within-control patient positions and (for
#'   markers) observed values. Empty histories get empty vectors.
#' @details Rejects unresolved duplicate biomarker records for the same
#'   patient and time. Multiple admission records are governed by the input
#'   contract. One time sort suffices; O(number of events + patients) storage.
#' @keywords internal
#' @noRd
.pair_events <- function(patients, events, ti, ci, biomarker) {
  nt <- length(ti)
  emp <- rep(list(numeric()), nt)
  o <- list(th = emp, tb = emp, tv = emp, ch = numeric(), chi = integer(),
            cb = numeric(), cbi = integer(), cv = numeric())
  if (is.null(events) || !nrow(events)) return(o)
  e <- events[events$type == "hospitalization" |
                (biomarker & events$type == "biomarker"), , drop = FALSE]
  e <- e[e$time <= patients$followup[match(e$id, patients$id)], , drop = FALSE]
  if (!nrow(e)) return(o)
  # One global sort is sufficient, even when the caller shuffles event rows.
  e <- e[order(e$time), , drop = FALSE]
  tr <- match(e$id, patients$id[ti])
  co <- match(e$id, patients$id[ci])
  h <- e$type == "hospitalization"
  b <- !h
  if (any(b) && anyDuplicated(e[b, c("id", "time"), drop = FALSE]))
    stop("Multiple biomarker values for one patient at the same time; resolve their order before analysis.",
         call. = FALSE)
  htr <- h & !is.na(tr); btr <- b & !is.na(tr)
  hco <- h & !is.na(co); bco <- b & !is.na(co)
  bucket <- function(idx) factor(tr[idx], levels = seq_len(nt))
  o$th <- unname(split(e$time[htr], bucket(htr)))
  if (biomarker) {
    o$tb <- unname(split(e$time[btr], bucket(btr)))
    o$tv <- unname(split(e$value[btr], bucket(btr)))
    o$cb <- e$time[bco]; o$cbi <- co[bco]; o$cv <- e$value[bco]
  }
  o$ch <- e$time[hco]; o$chi <- co[hco]
  o
}

#' Pool win and loss counts over patients sharing a settlement time
#'
#' @param counts Matrix with one row per endpoint owner and two columns,
#'   the numbers of treatment wins and losses assigned to that endpoint.
#' @param slots Integer vector mapping each owner to its position on the
#'   sorted union of endpoint times.
#' @param ngrid Number of distinct positions in that union.
#' @return An ngrid by two numeric matrix of counts in grid order, with
#'   zero rows for empty positions.
#' @details Treatment owns a pair if its endpoint equals the control endpoint.
#'   This helper groups counts once per arm, after the pair loop, without
#'   storing each pair. Its inputs must come from the same scoring pass.
#' @keywords internal
#' @noRd
.pair_endpoint_totals <- function(counts, slots, ngrid) {
  a <- rowsum(counts, slots, reorder = FALSE)
  out <- matrix(0, ngrid, ncol(counts))
  out[as.integer(rownames(a)), ] <- a
  out
}

#' Accumulate sufficient summaries from the clinical win and loss comparisons
#'
#' @param patients,events Validated patient and event snapshots. Follow-up
#'   and event times use the same time-since-entry units. See .pair_events.
#' @param ti,ci Indices of the treatment and control patients in patients.
#' @param tau Common follow-up horizon for IPCW; ignored for the shared window.
#' @param restricted Logical: FALSE compares through the shorter observed
#'   follow-up; TRUE requires ascertainment through the first death or tau.
#' @param km For a restricted analysis, list of the two supported .pair_km
#'   fits in treatment, control order. Ignored when restricted is FALSE.
#' @param biomarker_margin NULL omits the biomarker tier; otherwise a lower
#'   last-observed value must beat this nonnegative margin strictly.
#' @param hierarchy Prespecified gpc_hierarchy; NULL uses the legacy biomarker margin.
#' @param channels Optional named lists of event and measurement tables validated at the cut.
#' @return mean contains the two averages over all pairs, win then loss;
#'   average is the patient-by-two matrix of opposite-arm averages in input
#'   patient order. n_pairs and n_ascertained distinguish the denominator
#'   from the number observed. For IPCW, grid contains sorted observed pair
#'   endpoints and mass is a grid-by-two matrix of cutoff-specific weighted
#'   totals divided by the total number of pairs.
#' @details The first death by the comparison time settles the score.
#'   For supported restricted IPCW tiers, risk counts and ranks
#'   give all patient and time-specific totals. For count-only shared windows,
#'   each observed follow-up end is evaluated against the patients still
#'   under observation in the opposite arm. Other rules compare pairs
#'   directly. Both calculations return identical fields to the inference
#'   routines. A required missing measurement defers the analysis.
#' @keywords internal
#' @noRd
.pair_summaries <- function(patients, events, ti, ci, tau, restricted, km, biomarker_margin,
                            hierarchy = NULL, channels = NULL) {
  if (.shared_count_eligible(restricted, biomarker_margin, hierarchy))
    return(.shared_count_summaries(patients, events, ti, ci))
  if (.rank_eligible(restricted, biomarker_margin, hierarchy))
    return(.rank_ipcw(patients, events, ti, ci, tau, km))
  if (restricted) {
    specs <- .rank_tiers(hierarchy)
    if (!is.null(specs)) {
      .pair_check_endpoint_channels(hierarchy, channels)
      # Presets use the clinical channels selected for the native event view.
      source_channels <- if (!is.null(hierarchy$endpoints) || !is.null(hierarchy$compare)) channels
      result <- .rank_ipcw(patients, events, ti, ci, tau, km, specs, source_channels)
      if (!is.null(result)) return(result)
      # Missing deciding values require patient- or channel-specific diagnostics.
      # Delegate their ordering and wording to the same comparison as other fits.
    }
  }
  .pair_summaries_reference(patients, events, ti, ci, tau, restricted, km,
                            biomarker_margin, hierarchy, channels)
}

#' Evaluate comparisons directly for any supported clinical hierarchy
#' @inheritParams .pair_summaries
#' @return The patient, time and pair totals described in .pair_summaries.
#' @details Only ascertained pairs with two living members reach the clinical
#'   score provider. One row of pairs is held in memory. Win and loss totals
#'   share both the clinical comparisons and, when relevant, the two KM fits.
#' @noRd
.pair_summaries_reference <- function(patients, events, ti, ci, tau, restricted,
                                      km, biomarker_margin, hierarchy = NULL, channels = NULL) {
  nt <- length(ti); nc <- length(ci)
  xt <- patients$followup[ti]; xc <- patients$followup[ci]
  dt <- patients$death[ti]; dc <- patients$death[ci]
  tx <- if (restricted) pmin(xt, tau) else xt
  cx <- if (restricted) pmin(xc, tau) else xc
  txs <- matrix(0, nt, 2L); cxs <- matrix(0, nc, 2L)
  if (restricted) {
    grid <- sort(unique(c(tx, cx)))
    trslot <- match(tx, grid); cslot <- match(cx, grid)
    iw <- 1 / (.pair_G(km[[1L]], grid) * .pair_G(km[[2L]], grid))
    tct <- matrix(0, nt, 2L); cct <- matrix(0, nc, 2L)
  }
  # Compare observed survivors under the chosen clinical rule. The rest of
  # this pass uses only signed scores and observation times.
  order_t <- order(tx)
  if (is.null(hierarchy)) hierarchy <- gpc_hierarchy(admission_time = "none", biomarker_margin = biomarker_margin)
  score_alive <- .pair_score_provider(patients, events, ti, ci, tx[order_t], hierarchy, channels, restricted)
  nqual <- 0
  for (r in seq_along(order_t)) {
    i <- order_t[r]
    at <- pmin(tx[i], cx)
    tdead <- dt[i] & xt[i] <= at
    cdead <- dc & xc <= at
    qual <- if (restricted) (at >= tau | tdead | cdead) else rep(TRUE, nc)
    nqual <- nqual + sum(qual)
    # Death settles the sign; two observed deaths at the same time give 0.
    score <- as.integer(cdead) - as.integer(tdead)
    living <- which(qual & !tdead & !cdead)
    if (length(living)) score[living] <- score_alive(r, i, living, at[living])
    win <- score == 1L; loss <- score == -1L
    if (restricted) {
      by_t <- trslot[i] <= cslot
      w <- iw[pmin(trslot[i], cslot)]
      rw <- w * win; rl <- w * loss
      # Attribute each pair to one of its two endpoints (ties to treatment).
      # No n-pairs by n-patients array or full pair-score matrix is needed.
      tct[i, ] <- c(sum(win & by_t), sum(loss & by_t))
      cct[, 1L] <- cct[, 1L] + (win & !by_t)
      cct[, 2L] <- cct[, 2L] + (loss & !by_t)
    } else {
      rw <- win; rl <- loss
    }
    txs[i, ] <- c(sum(rw), sum(rl))
    cxs[, 1L] <- cxs[, 1L] + rw
    cxs[, 2L] <- cxs[, 2L] + rl
  }
  N <- as.double(nt) * nc
  averages <- matrix(0, nrow(patients), 2L)
  averages[ti, ] <- txs / nc
  averages[ci, ] <- cxs / nt
  p <- colSums(txs) / N
  out <- list(mean = p, average = averages, n_pairs = N, n_ascertained = nqual)
  if (restricted) {
    out$grid <- grid
    out$mass <- (.pair_endpoint_totals(tct, trslot, length(grid)) +
                   .pair_endpoint_totals(cct, cslot, length(grid))) * (iw / N)
  }
  out
}

#' Compute the two-sample projection for a shared-window signed win--loss score
#'
#' @param q Shared-window output from .pair_summaries.
#' @return Numeric vector in input patient order; its entries are the mean
#'   signed score against the patient's opposite arm minus the overall mean.
#'   Control scores retain the treatment-win orientation.
#' @details Averages over each arm are zero. Normalize by that arm's
#'   enrollment when assembling the covariance; no estimated weights enter.
#'   O(number of patients) arithmetic and storage.
#' @keywords internal
#' @noRd
.if_shared_win_loss <- function(q) {
  q$average[, 1L] - q$average[, 2L] - (q$mean[1L] - q$mean[2L])
}

#' Add the influence of one fitted reverse-KM curve to a scalar pair contrast
#'
#' @param x,died Follow-up times and death indicators in one arm, aligned
#'   with the arm sample used for km. Positions in the result follow x.
#' @param km Supported arm-specific fit from .pair_km. No unit jump may
#'   occur before the horizon; this is checked before weights are used.
#' @param grid Sorted cutoff times from a restricted .pair_summaries fit.
#' @param mass_after Reverse cumulative scalar weighted score mass on grid,
#'   ending with an appended zero. Its jth entry is the all-pair mean
#'   contributed by cutoffs at grid[j] or later.
#' @return One patient contribution from estimating the censoring weight,
#'   in the input arm order (zero if there are no relevant censoring jumps).
#' @details At jump s, take the score mass at cutoffs strictly greater
#'   than s, divide by mean risk times one minus the censoring hazard jump,
#'   and multiply by the individual censoring count minus its fitted risk
#'   contribution. A cumulative sum evaluates all risk terms up to each
#'   patient's observed end; its own censoring jump adds the count term.
#'   The positive sign here arises because the derivative of the inverse
#'   censoring survival cancels the negative KM survival derivative. No
#'   subject-by-jump array is formed.
#' @keywords internal
#' @noRd
.pair_weight_if <- function(x, died, km, grid, mass_after) {
  if (!length(km$time)) return(numeric(length(x)))
  after <- mass_after[1L + findInterval(km$time, grid)]
  coeff <- after / ((km$risk / km$n) * (1 - km$jump))
  before <- findInterval(x, km$time, left.open = TRUE)
  ans <- -c(0, cumsum(coeff * km$jump))[1L + before]
  ix <- match(x, km$time)
  cens <- which(!died & !is.na(ix))
  if (length(cens))
    ans[cens] <- ans[cens] + coeff[ix[cens]] * (1 - km$jump[ix[cens]])
  ans
}

#' Combine fixed-weight pair and fitted-censoring contributions for a linear contrast
#'
#' @param q Restricted win and loss summaries from .pair_summaries.
#' @param contrast Coefficients for win and loss, in that order, or the
#'   gradient of a smooth log contrast evaluated at their fitted means.
#' @param patients The same patient snapshot and order used in q.
#' @param arms List of two patient index vectors, treatment then control.
#' @param km Corresponding list of fitted censoring curves.
#' @return Numeric vector in input patient order, already on the contrast
#'   scale but not divided by enrollment or by the standard error.
#' @details Center the opposite-arm means at the overall linear contrast.
#'   Then add the appropriate arm-specific KM contribution. Linearity lets
#'   us form the transform gradient before the censoring derivative, keeping
#'   the covariance of wins and losses without refitting the two kernels.
#'   Extra storage is linear in patients and distinct cutoff/jump times.
#' @keywords internal
#' @noRd
.if_ipcw_linear <- function(q, contrast, patients, arms, km) {
  phi <- drop(q$average %*% contrast) - sum(q$mean * contrast)
  m <- drop(q$mass %*% contrast)
  after <- c(rev(cumsum(rev(m))), 0)
  for (a in 1:2) {
    j <- arms[[a]]
    phi[j] <- phi[j] + .pair_weight_if(patients$followup[j], patients$death[j],
                                        km[[a]], q$grid, after)
  }
  phi
}

#' Obtain patient contributions for net benefit or its log win-odds transform
#'
#' @param q,patients,arms,km The restricted summaries, aligned patient records,
#'   arm indices and KM fits used by .if_ipcw_linear.
#' @param log_odds If FALSE, return net benefit contributions. If TRUE,
#'   multiply those contributions by 2 / (1 - benefit squared).
#' @return Numeric vector in patient order, on the net-benefit or log-odds
#'   scale. No rescaling by enrollment is applied.
#' @details Both effects use the same signed kernel. Log win odds is not
#'   defined for a fitted benefit outside the open interval (-1,1): defer
#'   that effect rather than truncate or add pseudocounts.
#' @keywords internal
#' @noRd
.if_ipcw_benefit <- function(q, patients, arms, km, log_odds = FALSE, base_if = NULL) {
  b <- q$mean[1L] - q$mean[2L]
  if (log_odds && (!is.finite(b) || 1 - abs(b) <= 32 * .Machine$double.eps))
    .gpc_defer("Log win odds require net benefit strictly between -1 and 1 at working precision.")
  phi <- if (is.null(base_if)) .if_ipcw_linear(q, c(1, -1), patients, arms, km) else base_if
  if (log_odds) (2 / ((1 - b) * (1 + b))) * phi else phi
}

#' Obtain patient contributions for the log IPCW win ratio
#'
#' @param q,patients,arms,km The restricted summaries, aligned patient records,
#'   arm indices and KM fits used by .if_ipcw_linear.
#' @return Numeric vector in patient order on the log win-ratio scale;
#'   the influence of a win is divided by its fitted win probability,
#'   and the corresponding loss term divided by its fitted loss probability
#'   is subtracted. No rescaling by enrollment is applied.
#' @details Both component means must be positive. Otherwise defer only this
#'   transform. The two components use the same ascertainment times and
#'   censoring curves and their covariance is retained automatically.
#' @keywords internal
#' @noRd
.if_ipcw_win_ratio <- function(q, patients, arms, km) {
  if (any(q$mean <= 0)) .gpc_defer("Log win ratio requires positive estimated win and loss probabilities.")
  # The delta-method contrast combines BOTH components, including both KM
  # derivatives, before constructing the variance or cross-look covariance.
  .if_ipcw_linear(q, c(1 / q$mean[1L], -1 / q$mean[2L]), patients, arms, km)
}

#' Prepare all clinical and censoring summaries for one pairwise analysis
#'
#' @param patients,events Validated patient and clinical snapshots at this
#'   look; see .pair_fit for fields, units and data-cut conventions.
#' @param method Either "shared" or "ipcw".
#' @param tau,biomarker_margin Follow-up horizon (IPCW only) and optional
#'   strict margin on the final comparison tier.
#' @param hierarchy Prespecified gpc_hierarchy; NULL uses the legacy biomarker margin.
#' @param channels Optional named lists of event and measurement tables validated at the cut.
#' @param min_risk,min_g IPCW support limits for risk count and censoring
#'   survival at tau in each arm.
#' @return A reusable list: q is the paired sums from .pair_summaries, arms
#'   are treatment/control row indices, km are the corresponding censoring
#'   fits (NULL for shared window), and diagnostics contain ascertainment
#'   counts, component means, and (for IPCW) horizon support.
#' @details This is the only clinical-scoring pass for a look, regardless
#'   of the number of descriptive IPCW effects requested. Empty arms or
#'   insufficient support defer all contrasts using this fit.
#' @keywords internal
#' @noRd
.pair_core <- function(patients, events, method, tau, biomarker_margin, min_risk, min_g,
                       hierarchy = NULL, channels = NULL) {
  restricted <- method == "ipcw"
  arms <- list(which(patients$arm == 1L), which(patients$arm == 0L))
  if (any(lengths(arms) == 0L)) .gpc_defer("Both treatment arms must be enrolled.")
  km <- NULL
  diag <- list()
  if (restricted) {
    if (length(tau) != 1L || !is.finite(tau) || tau < 0)
      stop("A restricted comparison needs one finite, nonnegative horizon.", call. = FALSE)
    km <- lapply(arms, function(ix) .pair_km(patients$followup[ix], patients$death[ix], tau))
    risk <- vapply(arms, function(ix) sum(patients$followup[ix] >= tau), integer(1L))
    gh <- vapply(km, .pair_G, numeric(1L), t = tau)
    names(gh) <- names(risk) <- c("treatment", "control")
    diag$risk_at_horizon <- risk; diag$G_at_horizon <- gh
    if (any(risk < min_risk)) .gpc_defer("Insufficient patients at risk at the IPCW horizon.")
    if (any(!is.finite(gh) | gh <= 0 | gh < min_g))
      .gpc_defer("Insufficient estimated censoring survival at the IPCW horizon.")
  }
  q <- .pair_summaries(patients, events, arms[[1L]], arms[[2L]], tau,
                       restricted, km, biomarker_margin, hierarchy, channels)
  diag$n_pairs <- q$n_pairs; diag$n_ascertained <- q$n_ascertained
  diag$win <- q$mean[1L]; diag$loss <- q$mean[2L]
  list(q = q, arms = arms, km = km, diagnostics = diag)
}

#' Derive a single pairwise effect and its patient influences from a common fit
#'
#' @param core Common summaries and estimated censoring fits from .pair_core.
#' @param patients Input snapshot in its original row order.
#' @param effect Selected net benefit, win odds or win ratio; factors and
#'   weights are not re-estimated.
#' @param base_if Optional influence for the net benefit on these records
#'   when forming its log-odds transformation.
#' @return Test-scale estimate, natural-scale estimate, estimated influence
#'   vector on the test scale and common pair/support diagnostics.
#' @details Ordinarily only the requested influence is computed. An undefined
#'   log contrast raises gpc_deferred; a descriptive multi-effect caller may
#'   retain valid results for the other contrasts.
#' @keywords internal
#' @noRd
.pair_effect <- function(core, patients, effect, base_if = NULL) {
  q <- core$q
  b <- q$mean[1L] - q$mean[2L]
  if (is.null(core$km)) {
    estimate <- natural <- b
    phi <- .if_shared_win_loss(q)
  } else if (effect == "win_ratio") {
    phi <- .if_ipcw_win_ratio(q, patients, core$arms, core$km)
    estimate <- log(q$mean[1L]) - log(q$mean[2L]); natural <- exp(estimate)
  } else {
    use_odds <- effect == "win_odds"
    phi <- .if_ipcw_benefit(q, patients, core$arms, core$km, use_odds, base_if)
    estimate <- if (use_odds) log1p(b) - log1p(-b) else b
    natural <- if (use_odds) exp(estimate) else b
  }
  list(estimate = estimate, natural = natural, influence = unname(phi),
       diagnostics = core$diagnostics)
}

#' Fit one look of a shared-window or restricted pairwise comparison
#'
#' @param patients Validated enrolled snapshot: unique id, arm (1 treatment,
#'   0 control), followup (time since entry) and death (logical). Row order
#'   is preserved throughout. Prior and future looks are never accessed.
#' @param events Validated available event records as described in .pair_events.
#' @param method Either "shared" (shorter observed follow-up) or "ipcw"
#'   (first death or common horizon).
#' @param effect "net_benefit", "win_odds" or "win_ratio". Shared window
#'   supports only the signed score, whose mean is censoring dependent.
#' @param tau Nonnegative follow-up horizon for a restricted analysis.
#' @param biomarker_margin Optional nonnegative margin for the last biomarker
#'   comparison after survival and hospitalizations. NULL omits that tier.
#' @param hierarchy Prespecified gpc_hierarchy; NULL uses the legacy biomarker margin.
#' @param channels Optional named lists of event and measurement tables validated at the cut.
#' @param min_risk Minimum number of observations at risk at the horizon
#'   required in each arm for IPCW (inclusive of that time).
#' @param min_g Minimum fitted censoring survival at the horizon in each arm.
#'   Positive survival is required even if min_g equals zero.
#' @return List with estimate on the test scale, natural on its original
#'   scale, influence in input patient order, and diagnostics: total pairs,
#'   ascertained pairs, win and loss means, plus arm-specific horizon risk
#'   counts and fitted censoring survival for IPCW. Both ratio effects use
#'   log estimates and log-scale influences.
#' @details The public design and snapshot validators handle malformed input.
#'   An empty arm, inadequate support, missing deciding biomarker or undefined
#'   chosen transform signals a deferral. A nonpositive variance is checked
#'   by the monitoring layer. One KM fit per arm and one score pass suffice.
#' @keywords internal
#' @noRd
.pair_fit <- function(patients, events, method, effect, tau, biomarker_margin = NULL,
                      min_risk = 1L, min_g = 0, hierarchy = NULL, channels = NULL) {
  if (!(method %in% c("shared", "ipcw"))) stop("Unknown pair method.", call. = FALSE)
  if (!(effect %in% c("net_benefit", "win_odds", "win_ratio")) ||
      (method == "shared" && effect != "net_benefit"))
    stop("Unsupported effect for the selected pair method.", call. = FALSE)
  core <- .pair_core(patients, events, method, tau, biomarker_margin, min_risk, min_g, hierarchy, channels)
  .pair_effect(core, patients, effect)
}
