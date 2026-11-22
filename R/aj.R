# Exact arm-specific Aalen--Johansen fits and patient influences for the
# restricted net win time. Times below are follow-up times, not calendar dates.
# The derivatives in Supporting information C.3 are accumulated backwards:
# we need one scalar per patient, not an occupation curve for every patient.

#' Add increments by integer index
#'
#' @param x A numeric vector to update.
#' @param index Positive integer positions in `x`, possibly repeated.
#' @param value One signed increment per position.
#' @return A numeric vector of the original length and order. Repeated positions
#'   receive the sum of their increments; no positions means no change.
#' @details Uses direct indexed addition for unique positions and a grouped sum
#'   otherwise. Callers supply valid, nonmissing positions. The additions require
#'   work linear in the input increments, plus the output vector if copied.
#' @noRd
.aj_add <- function(x, index, value) {
  if (!length(index)) return(x)
  if (anyDuplicated(index)) {
    sums <- rowsum(matrix(value, ncol = 1L), index, reorder = FALSE)
    index <- as.integer(rownames(sums))
    value <- sums[, 1L]
  }
  x[index] <- x[index] + value
  x
}

#' Validate and index exact transition records
#'
#' @param p One row per patient enrolled at this look, with unique `id`,
#'   initial state `state0`, nonnegative observed `followup` and logical `death`.
#' @param transitions A data frame with `id,time,from,to`, or `NULL` for no
#'   transitions. Times are measured from entry; history at zero is in `state0`.
#' @param states Distinct state labels, ordered best to worst.
#' @param absorbing Names of states that cannot be left, possibly empty.
#' @param death_state Explicit absorbing label for death, or `NULL` in the
#'   absence of observed terminal deaths. Other absorbing labels may be nonfatal.
#' @return A data frame ordered by patient and follow-up time, with `id` now
#'   the row position in `p`, and `from,to` the positions in `states`. Every
#'   recorded transition is validated, including those after the target horizon.
#' @details Rejects unknown patients or states, changes outside observed
#'   follow-up, conflicting/simultaneous changes in one patient, and inconsistent
#'   death or absorbing-state records. With M transitions and n enrollees,
#'   workspace is O(n+M); ordering is at most O(M log M).
#' @noRd
.aj_events <- function(p, transitions, states, absorbing, death_state) {
  if (!is.character(states) || !length(states) || anyNA(states) ||
      any(states == "") || anyDuplicated(states))
    stop("state_order must list distinct state names, best to worst.", call. = FALSE)
  if (!is.character(absorbing) || anyNA(absorbing) ||
      any(!absorbing %in% states) || anyDuplicated(absorbing))
    stop("absorbing must contain distinct names in state_order.", call. = FALSE)
  if (!is.null(death_state) &&
      (!is.character(death_state) || length(death_state) != 1L ||
       is.na(death_state) || !death_state %in% absorbing))
    stop("death_state must be one of the declared absorbing states.", call. = FALSE)
  if (any(p$death) && is.null(death_state))
    stop("An observed death requires an explicitly declared absorbing death_state.", call. = FALSE)
  if (!"state0" %in% names(p)) stop("patients must include state0 for the AJ method.", call. = FALSE)
  initial <- match(p$state0, states)
  if (anyNA(initial))
    stop("Unknown or missing state0 for patient ", p$id[which(is.na(initial))[1L]], ".", call. = FALSE)

  if (is.null(transitions))
    transitions <- data.frame(id = character(), time = numeric(), from = character(), to = character())
  need <- c("id", "time", "from", "to")
  if (!is.data.frame(transitions) || !all(need %in% names(transitions)))
    stop("transitions must contain id, time, from, and to.", call. = FALSE)
  sid <- match(transitions$id, p$id)
  if (anyNA(sid))
    stop("Unknown transition patient ID: ", transitions$id[which(is.na(sid))[1L]], ".", call. = FALSE)
  t <- transitions$time
  if (!is.numeric(t) || any(!is.finite(t)))
    stop("Transition time must be a finite numeric follow-up time.", call. = FALSE)
  f <- match(transitions$from, states); z <- match(transitions$to, states)
  ord <- order(sid, t)
  sid <- sid[ord]; t <- t[ord]; f <- f[ord]; z <- z[ord]
  fail <- function(bad, why) {
    if (any(bad)) {
      j <- which(bad)[1L]
      stop("Invalid transition for patient ", p$id[sid[j]], " at time ", t[j],
           ": ", why, ".", call. = FALSE)
    }
  }
  if (length(sid)) {
    fail(t <= 0, "positive time required; record the state at entry in state0")
    fail(t > p$followup[sid], "transition after the end of observed follow-up")
    fail(is.na(f) | is.na(z), "unknown or missing from or to state")
    first <- !duplicated(sid)
    before_t <- c(0, t[-length(t)])
    fail(!first & t == before_t, "multiple transitions at one time; adjudicate their exact times")
    fail(f == z, "from and to must differ")
    expected <- initial[sid]; last <- c(0L, z[-length(z)])
    expected[!first] <- last[!first]
    fail(f != expected, "from does not match the previous state")
    fail(f %in% match(absorbing, states), "departure from an absorbing state")
  }
  if (!is.null(death_state)) {
    d <- match(death_state, states)
    flag <- initial == d & p$followup == 0
    if (any(initial == d & p$followup != 0))
      stop("A patient in death_state at entry must have followup zero.", call. = FALSE)
    at_death <- z == d
    if (any(at_death)) {
      fail(at_death & t != p$followup[sid], "entry to death_state must be at the observed death time")
      flag[sid[at_death]] <- TRUE
    }
    if (any(flag != p$death))
      stop("death disagrees with the recorded death_state at followup for patient ",
           p$id[which(flag != p$death)[1L]], ".", call. = FALSE)
  }
  data.frame(id = sid, time = t, from = f, to = z)
}

# One grid for both arms, including censoring times so that a state losing
# its last observed patient midway through an interval cannot go unnoticed.
#' Group records by grid position, retaining empty positions
#'
#' @param index Integer grid position of each input record, between 1 and L.
#' @param L Number of grid positions.
#' @return Unnamed list of length `L`; element l contains input row positions
#'   belonging to grid position l, in input order, or an empty integer vector.
#' @details Used for transition and censoring times. An empty index is allowed.
#'   Matching and filtering to valid grid positions are the caller's responsibility.
#'   Work and memory are O(L + length(index)).
#' @noRd
.aj_group <- function(index, L) {
  unname(split(seq_along(index), factor(index, levels = seq_len(L))))
}

#' Fit arm-specific occupation probabilities on a shared integration grid
#'
#' @param initial Integer state positions at entry, one per patient in arm order.
#' @param exits Observed terminal/censoring follow-up times in the same order.
#' @param e Validated transitions strictly before the horizon, in patient/time
#'   order. Patient indices are local to this arm; state positions are global.
#' @param grid Distinct interval starts, beginning at zero and including all
#'   transition and observation-exit times strictly before the horizon.
#' @param states Ordered names of the states, best to worst.
#' @param absorbing Names of clinically absorbing states.
#' @param arm Arm identifier (0 or 1) for support diagnostics.
#' @return A list containing the interval-by-state occupation matrix `p`, entry
#'   and last observed pre-horizon state per patient, exit times, a jump-to-event
#'   index, the augmented event table (`g` grid position, `risk` pretransition
#'   source-state count, `factor` equal to arm size times the old source-state
#'   probability divided by that count), and the minimum occupied transient-state
#'   risk on positive-length intervals (NA if none).
#' @details Tied transitions of different patients use the same pretransition
#'   counts; observation exits at their recorded times follow the clinical
#'   updates. A positive occupation probability in a transient state with no
#'   observed patient over an interval defers the look. For S states, L intervals,
#'   n enrollees and M transitions: O(n+M+SL) work after indexing, with O(n+M+SL)
#'   storage. No patient-by-time-by-state array is constructed.
#' @noRd
.aj_forward <- function(initial, exits, e, grid, states, absorbing, arm) {
  m <- length(initial); S <- length(states); L <- length(grid)
  current <- initial; risk <- tabulate(current, S); p <- risk / m
  occ <- matrix(0, L, S)
  e$g <- match(e$time, grid)
  jumps <- .aj_group(e$g, L)
  exit_index <- match(exits, grid)
  exit_rows <- which(!is.na(exit_index))
  departures <- .aj_group(exit_index[exit_rows], L)
  event_risk <- event_factor <- numeric(nrow(e))
  from <- e$from; to <- e$to; id <- e$id
  minrisk <- Inf
  transient <- !states %in% absorbing
  for (l in seq_len(L)) {
    j <- jumps[[l]]
    if (length(j) == 1L) {
      # One transition changes just two coordinates of the occupation and risk.
      f <- from[j]; z <- to[j]
      nr <- risk[f]
      if (nr <= 0L) stop("AJ transition exceeds its pretransition risk set.", call. = FALSE)
      step <- p[f] / nr
      event_risk[j] <- nr; event_factor[j] <- m * step
      p[f] <- p[f] * ((nr - 1L) / nr)
      p[z] <- p[z] + step
      current[id[j]] <- z
      risk[f] <- nr - 1L; risk[z] <- risk[z] + 1L
    } else if (length(j)) {
      f <- from[j]; z <- to[j]
      nr <- risk[f]; out <- tabulate(f, S)
      if (any(nr <= 0L) || any(out > risk)) stop("AJ transition exceeds its pretransition risk set.", call. = FALSE)
      # Count clinical transitions before removing patients whose observation
      # ends at this same time. All tied patients use the same old risk sets.
      step <- p[f] / nr
      event_risk[j] <- nr; event_factor[j] <- m * step
      ff <- which(out > 0)
      p[ff] <- p[ff] * ((risk[ff] - out[ff]) / risk[ff])
      p <- .aj_add(p, z, step)
      current[id[j]] <- z
      risk <- risk - out + tabulate(z, S)
    }
    ending <- exit_rows[departures[[l]]]
    if (length(ending)) risk <- risk - tabulate(current[ending], S)
    required <- transient & p > 0
    if (any(required & risk == 0L)) {
      bad <- which(required & risk == 0L)[1L]
      .gpc_defer(paste0("AJ: no patients under observation in occupied transient state ",
                        states[bad], " in arm ", arm, " after follow-up time ", grid[l], "."))
    }
    if (any(required)) minrisk <- min(minrisk, risk[required])
    occ[l, ] <- p
  }
  e$risk <- event_risk; e$factor <- event_factor
  list(p = occ, initial = initial, final = current, exits = exits, events = e,
       jumps = jumps, minrisk = if (is.finite(minrisk)) minrisk else NA_integer_)
}

# Derivative w.r.t. each state's occupation when the other arm has the
# supplied occupation probabilities. A smaller index is a better state.
#' Compute each state's instantaneous signed comparison
#'
#' @param opposite Matrix of fitted opposite-arm state probabilities, rows
#'   for integration intervals, columns ordered best to worst.
#' @return Matrix of the same dimensions: for a row and a candidate state,
#'   probability the opposite arm is worse minus probability it is better.
#'   Equality of state contributes zero. Negate for a control patient's effect on
#'   the treatment-minus-control contrast.
#' @details Uses cumulative probabilities instead of a square comparison matrix;
#'   work and storage are O(L S) for L intervals and S states.
#' @noRd
.aj_payoff <- function(opposite) {
  mass <- rowSums(opposite)
  better <- numeric(nrow(opposite)); ans <- opposite
  for (s in seq_len(ncol(opposite))) {
    ans[, s] <- mass - opposite[, s] - 2 * better
    better <- better + opposite[, s]
  }
  ans
}

# Multiply the case-weight derivative by the arm size as in C.3. Work
# backwards to compute the effect of each occupation on later net win time.
# At each jump, account separately for patients who made a transition and
# those at risk. Prefix sums apply the risk correction to observed intervals.
#' Calculate an arm's patient influences on integrated net win time
#'
#' @param fit One arm's output from `.aj_forward`.
#' @param payoff Interval-by-state derivative of the instantaneous contrast
#'   with respect to this arm's occupation probabilities, with its arm sign.
#' @param grid Increasing interval starts, beginning at zero.
#' @param widths Positive lengths through the horizon, aligned with `grid`.
#' @return Numeric vector of arm-size-scaled case-weight derivatives in the
#'   same patient order as the arm's initial states. Contributions sum to zero
#'   within arm up to rounding, including in the presence of censoring.
#' @details Reverse differentiation of the forward AJ product equals the
#'   patient recursion in Supporting information C.3. Transition residuals and
#'   risk-interval prefix sums avoid per-patient refits or per-patient propagation
#'   at every time. After grid matching, the core sums cost O(n+M+SL);
#'   finding patients' endpoint indices costs O(n log L) in the general case.
#'   Workspace is O(n+M+SL). Support and chronology are checked by the forward fit.
#' @noRd
.aj_influence <- function(fit, payoff, grid, widths) {
  n <- length(fit$initial); S <- ncol(fit$p); L <- length(grid)
  e <- fit$events
  from <- e$from; to <- e$to; factor <- e$factor; risk <- e$risk
  future <- numeric(S); charge <- matrix(0, L, S); gain <- numeric(nrow(e))
  for (l in rev(seq_len(L))) {
    h <- widths[l] * payoff[l, ] + future
    j <- fit$jumps[[l]]
    if (length(j) == 1L) {
      f <- from[j]
      diff <- h[to[j]] - h[f]
      gain[j] <- factor[j] * diff
      charge[l, f] <- gain[j] / risk[j]
      future <- h
      future[f] <- h[f] + diff / risk[j]
    } else if (length(j)) {
      f <- from[j]; diff <- h[to[j]] - h[f]
      gain[j] <- factor[j] * diff
      charge[l, ] <- .aj_add(numeric(S), f, gain[j] / risk[j])
      future <- .aj_add(h, f, diff / risk[j])
    } else future <- h
  }
  for (s in seq_len(S)) charge[, s] <- cumsum(charge[, s])
  phi <- future[fit$initial] - mean(future[fit$initial])
  last <- rep.int(1L, n)
  if (nrow(e)) {
    # Events are ordered by patient and time, with no simultaneous changes
    # for one patient. Risk intervals include their right endpoints.
    previous <- c(1L, e$g[-nrow(e)])
    previous[!duplicated(e$id)] <- 1L
    interval <- charge[cbind(e$g, e$from)] - charge[cbind(previous, e$from)]
    phi <- .aj_add(phi, e$id, gain - interval)
    last_event <- !duplicated(e$id, fromLast = TRUE)
    last[e$id[last_event]] <- e$g[last_event]
  }
  end <- findInterval(fit$exits, grid)
  phi - (charge[cbind(end, fit$final)] - charge[cbind(last, fit$final)])
}

#' Fit restricted net win time and its patient influences at one look
#'
#' @param patients Validated enrolled records with columns `id`, `arm`,
#'   `entry`, `followup`, `death` and `state0`. Influence values follow this row order.
#' @param transitions Exact clinical transitions available at this look;
#'   columns `id,time,from,to`, or `NULL` for none. Time is since entry.
#' @param tau Common nonnegative follow-up horizon in the units of `followup`.
#' @param state_order Finite vector of distinct clinical state names, best first.
#' @param absorbing Clinically absorbing states within `state_order`.
#' @param death_state Absorbing state representing documented deaths; `NULL` if
#'   there are none and no such label is used.
#' @param keep_occupation Whether to retain the two fitted right-continuous state
#'   probability curves, including transitions at the reported horizon.
#' @param prepared Optional internal table returned by `.aj_events` for exactly
#'   these patients, transitions, state ordering and terminal-state convention.
#'   The public snapshot validator supplies it before checking statistical support.
#'   If omitted, validate all submitted paths here.
#' @return A list: `estimate` and `natural` are the same signed net win time
#'   (time units); `influence` is numeric with one value per input row;
#'   `diagnostics` contains the horizon, interval count, pre-horizon transition
#'   counts, and minimum occupied transient-state risks for control/treatment.
#' @details Invalid paths are input errors. Empty arms, a zero horizon, and
#'   insufficient state support defer the look. The monitoring layer handles
#'   zero estimated variance. No Markov restriction is imposed on the clinical
#'   history. Sorting and interval matching supplement the O(n+M+SL) forward
#'   and reverse sums, for n enrollees, M changes, S states and L intervals.
#' @noRd
.aj_fit <- function(patients, transitions, tau, state_order, absorbing,
                    death_state = NULL, prepared = NULL, keep_occupation = FALSE) {
  if (length(tau) != 1L || !is.numeric(tau) || !is.finite(tau) || tau < 0)
    stop("AJ horizon must be one finite nonnegative follow-up time.", call. = FALSE)
  e <- if (is.null(prepared)) .aj_events(patients, transitions, state_order, absorbing, death_state) else prepared
  n <- nrow(patients); index <- lapply(0:1, function(a) which(patients$arm == a))
  if (any(lengths(index) == 0L)) .gpc_defer("AJ needs at least one enrolled patient in each arm.")
  if (tau == 0) .gpc_defer("AJ net win time has zero variance at horizon zero.")
  if (keep_occupation) endpoint <- e[e$time == tau, , drop = FALSE]
  e <- e[e$time < tau, , drop = FALSE]  # a change at tau has zero duration
  grid <- sort(unique(c(0, e$time, patients$followup[patients$followup < tau])))
  width <- diff(c(grid, tau))
  local_id <- integer(n)
  for (ii in index) local_id[ii] <- seq_along(ii)
  ee <- split(e, factor(patients$arm[e$id], levels = 0:1))
  initial <- match(patients$state0, state_order)
  fit <- lapply(1:2, function(b) {
    ii <- index[[b]]; eb <- ee[[b]]; eb$id <- local_id[eb$id]
    .aj_forward(initial[ii], patients$followup[ii], eb, grid, state_order, absorbing, b - 1L)
  })
  pay1 <- .aj_payoff(fit[[1L]]$p); pay0 <- -.aj_payoff(fit[[2L]]$p)
  estimate <- sum(width * rowSums(fit[[2L]]$p * pay1))
  influence <- numeric(n)
  influence[index[[1L]]] <- .aj_influence(fit[[1L]], pay0, grid, width)
  influence[index[[2L]]] <- .aj_influence(fit[[2L]], pay1, grid, width)
  result <- list(estimate = estimate, natural = estimate, influence = influence,
       diagnostics = list(horizon = tau, intervals = length(width),
                          transitions = stats::setNames(vapply(ee, nrow, integer(1)), c("control", "treatment")),
                          min_transient_risk = stats::setNames(vapply(fit, function(x) x$minrisk, numeric(1)),
                                                       c("control", "treatment"))))
  if (keep_occupation) {
    result$occupation <- lapply(1:2, function(b) {
      arm_fit <- fit[[b]]
      loc <- sort(unique(c(1L, arm_fit$events$g)))
      prob <- arm_fit$p[loc, , drop = FALSE]
      # The integral stops at tau. For a right-continuous probability at tau,
      # apply any endpoint jumps with the same pretransition risk convention.
      last <- prob[nrow(prob), ]
      j <- which(patients$arm[endpoint$id] == b - 1L)
      if (length(j)) {
        src <- endpoint$from[j]
        risk <- tabulate(arm_fit$final[arm_fit$exits >= tau], length(last))
        move <- last[src] / risk[src]
        out <- tabulate(src, length(last)); changed <- which(out > 0L)
        last[changed] <- last[changed] * ((risk[changed] - out[changed]) / risk[changed])
        last <- .aj_add(last, endpoint$to[j], move)
      }
      prob <- rbind(prob, last)
      dimnames(prob) <- list(NULL, state_order)
      list(time = c(grid[loc], tau), probability = prob)
    })
  }
  result
}
