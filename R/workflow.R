#' Fit a pairwise treatment effect at one look
#'
#' @description Fit one prespecified comparison and obtain one estimated influence per enrolled patient. A standalone fit neither spends error nor updates the monitoring history.
#' @param object A trial plan from \code{\link{gpc_design}}.
#' @param keep_occupation Logical. For an Aalen--Johansen analysis, use TRUE
#'   to retain the fitted state occupation curves for this look. Extract them
#'   with \code{\link{gpc_occupation}}. Defaults to FALSE.
#' @param look Scheduled integer index, from 1 through the number of planned looks.
#' @param data Optional replacement for \code{patients}: a cut from \code{\link{gpc_trial_data}} containing its own histories, or a native patient data frame with separate event or transition arguments. Supply only one of \code{patients} and \code{data}.
#' @param patients Data frame with one row per patient enrolled by this calendar look. Required columns: immutable, unique character \code{id}; numeric \code{arm} (1 treatment, 0 control); nonnegative numeric calendar \code{entry}; nonnegative numeric \code{followup} from entry to observed death or last complete observation at this cut; logical \code{death} indicating an observed death at \code{followup}. For \code{"aj"} also provide character \code{state0}, the state at entry.
#' @param events For win--loss comparisons, data frame of the clinical events known by this look, or \code{NULL} for none. Supply character \code{id}, numeric \code{time} since entry and character \code{type} (\code{"hospitalization"} or \code{"biomarker"}). Supply numeric \code{value} for biomarker measurements. Each hospitalization row counts as one admission. Optional numeric \code{available} is the calendar time the record became known.
#' @param transitions For \code{"aj"}, data frame of all changes of state observed by this look, or \code{NULL} if none. Supply character \code{id}, numeric \code{time} since entry, and character \code{from} and \code{to} among the declared states. Optional numeric \code{available} is the calendar time the transition became known.
#' @details Each data cut must contain all enrolled patients and their complete observed histories through \code{followup}. Use one time unit throughout. Do not include a patient or event from a later cut; \code{entry + followup} cannot exceed the scheduled calendar time. \code{available}, when supplied, must be between occurrence and the data cut. Without it the caller must still construct an as-of snapshot. An unreported interval is not evidence that no event occurred.
#'
#' For shared-window or IPCW comparisons, the first death settles the pair, including a death at the comparison time. Equal observed death times tie. While both are alive, the default favors fewer admissions, breaking a tie by later last admission. The plan can substitute time to first admission or omit the timing tier. An optional biomarker tier then compares the last available values at or before the comparison time. Absence of a required value is not a tie.
#'
#' For IPCW, ascertain a pair if both patients were observed through the horizon or until the first death in the pair, whichever occurs first. The denominator is the product of the arm-specific reverse Kaplan--Meier probabilities of observation through that time; censoring jumps at the evaluation time are excluded. At a known common date, observed deaths precede censoring in the censoring risk set; unresolved rounding must be adjudicated before analysis. The influences include the contribution from estimating the censoring distributions. Win odds is a transformation of net benefit; win ratio uses separate weighted win and loss probabilities. Ratios are fitted on the log scale. Undefined transformations are not replaced by pseudocounts.
#'
#' For Aalen--Johansen, record the initial state and all transitions: a patient's \code{from} state must agree with the previous state and there can be no departure from a declared absorbing state. The declared death state must agree with \code{death} and \code{followup}. Transitions of different patients on a common date use pretransition risk sets; adjudicate multiple changes in the same patient recorded at one time. Net win time integrates the signed difference in state favorability over the common follow-up horizon.
#' @returns A \code{gpcseq_fit} object. \code{estimate} is on the test scale; \code{natural} is on the original scale. They coincide for the signed score, net benefit and net win time. For win odds or ratio they are the log estimate and its exponential. Other components include \code{look}, calendar \code{time}, \code{horizon} (\code{NA} for a shared window), \code{method}, \code{effect}, \code{subjects} (patient \code{id}, \code{arm}, and test-scale influence \code{phi}) and method-specific \code{diagnostics}. An unsupported effect raises a \code{gpc_deferred} condition. An invalid input raises an error. A finite fit can still have zero variance and be unsuitable for monitoring.
#' @seealso \code{\link{gpc_design}}, \code{\link{gpc_update}}, \code{\link{gpc_influence}}.
#' @examples
#' cut <- gpc_example(2)
#' plan <- gpc_design(2, 0.05, "ipcw", effect = "win_ratio", horizons = 2)
#' fit <- gpc_fit(plan, 1, cut$patients, cut$events)
#' c(log_win_ratio = fit$estimate, win_ratio = fit$natural)
#' gpc_influence(fit)
#'
#' aj <- gpc_design(2, 0.05, "aj", horizons = 2,
#'                  state_order = c("well", "hospitalized", "dead"),
#'                  absorbing = "dead", death_state = "dead")
#' gpc_fit(aj, 1, cut$patients, transitions = cut$transitions)$estimate
#' @export
gpc_fit <- function(object, look, patients=NULL, events=NULL, transitions=NULL, data=NULL,
                    keep_occupation=FALSE) {
  d <- .gpc_check_look(object,look)
  x <- .gpc_data_snapshot(patients, events, transitions, data, d, look)
  horizon <- .gpc_horizon(d,look,x$p)
  .gpc_fit_cut(x,d,look,horizon,keep_occupation)
}

# Shared index and support checks for all public single-look interfaces.
.gpc_check_look <- function(object,look) {
  if (!inherits(object,"gpcseq")) stop("object must be a gpcseq design.",call.=FALSE)
  d <- object$design
  if (!is.numeric(look) || length(look)!=1L || !is.finite(look) ||
      look!=floor(look) || look<1L || look>length(d$looks))
    stop("look must be one scheduled look index between 1 and ",length(d$looks),".",call.=FALSE)
  d
}

.gpc_check_support <- function(p,d,horizon) {
  count <- tabulate(as.integer(p$arm)+1L,nbins=2L)
  if (any(count<d$min_enrolled)) .gpc_defer("Too few enrolled patients in at least one arm.")
  if (d$method!="shared" && (!is.finite(horizon) || horizon<=0))
    .gpc_defer("No positive supported horizon.")
  invisible(count)
}

.gpc_fit_cut <- function(x,d,k,horizon,keep_occupation=FALSE) {
  if (!identical(keep_occupation, FALSE) &&
      !(identical(keep_occupation, TRUE) && d$method == "aj"))
    stop("keep_occupation must be FALSE, or TRUE for an Aalen--Johansen fit.", call.=FALSE)
  .gpc_check_support(x$p,d,horizon)
  f <- if (d$method=="aj") .aj_fit(x$p,x$tr,horizon,d$state_order,d$absorbing,d$death_state,prepared=x$aj_events,keep_occupation=keep_occupation) else
    .pair_fit(x$p,x$e,d$method,d$effect,horizon,d$biomarker_margin,d$min_risk,d$min_g,hierarchy=d$hierarchy,channels=x$channels)
  if (length(f$influence)!=nrow(x$p) || any(!is.finite(f$influence)) ||
      length(f$estimate)!=1L || !is.finite(f$estimate))
    .gpc_defer("The chosen estimate or a subject influence is not finite.")
  f$subjects <- data.frame(id=x$p$id,arm=x$p$arm,phi=as.numeric(f$influence))
  f$influence <- NULL
  f$look <- k; f$time <- d$looks[k]; f$horizon <- horizon; f$method <- d$method
  f$effect <- d$effect; class(f) <- "gpcseq_fit"; f
}

#' Process the next scheduled look
#'
#' @description Update a design with the records available at its next look. The function fits the chosen effect, obtains covariance with earlier tests from the same patients' contributions, and calculates a stopping boundary.
#' @inheritParams gpc_fit
#' @param object An open \code{gpcseq} object, either the initial design or the object returned at the previous scheduled look.
#' @param data Optional replacement for \code{patients}: a cut from \code{\link{gpc_trial_data}} containing its own histories, or a native patient data frame with separate event or transition arguments. Supply only one of \code{patients} and \code{data}.
#' @param patients Complete snapshot of the enrolled patients at the next scheduled calendar time. Required fields and conventions are given in \code{\link{gpc_fit}}. Retain previously enrolled IDs; their arm, entry time, baseline state if applicable, and observed death are fixed.
#' @param events For shared-window and IPCW comparisons, the cumulative clinical observations available at this cut; see \code{\link{gpc_fit}}. Do not supply only the newly reported events.
#' @param transitions For Aalen--Johansen, the cumulative observed transitions available at this cut; see \code{\link{gpc_fit}}.
#' @details Always assign the returned object. Earlier influences and boundaries remain as first fitted. Among previously completed tests, let \eqn{j} be the last tested look; use zero previous spending if there was none. At look \eqn{k} the new boundary allocates unconditional probability \eqn{\alpha_k-\alpha_j}{alpha[k] - alpha[j]} to a first crossing at this look under the joint normal approximation. An upper crossing favors treatment; a lower crossing in a two-sided test favors control. In the hybrid design only the scheduled final look is two-sided, even if an interim was skipped.
#'
#' An unsupported look, undefined ratio, nonpositive estimated variance, singular correlation, or unsuccessful calibration spends no error. At the next scheduled look the data are assessed anew. If the scheduled final test cannot be performed and no earlier test rejected, the trial ends without rejection. Input errors, such as a changed arm or an impossible transition, are not converted into scientific deferrals. Once a boundary has been crossed, or the scheduled final look has been processed, no further update is allowed.
#' Previously supplied event and measurement records must remain in subsequent
#' cuts, including records supplied at a deferred look. Rows may be reordered
#' and new records added; removing records or changing their values or supplied
#' availability dates raises an error. Reconcile corrections against the saved
#' cuts before proceeding; past fitted decisions are not overwritten.
#' Older saved plans without retained histories issue a warning on their first
#' update: checking starts with that cut and cannot verify earlier records.
#'
#' @returns A new \code{gpcseq} object. \code{\link{gpc_results}} returns the decision at each processed look, including the reason when no test was performed. Use \code{\link{gpc_influence}} and \code{\link{gpc_covariance}} to inspect completed tests.
#' @seealso \code{\link{gpc_fit}}, \code{\link{gpc_results}}, \code{\link{gpc_example}}.
#' @examples
#' plan <- gpc_design(c(1, 2), c(0.01, 0.05), "shared")
#' for (k in 1:2) {
#'   cut <- gpc_example(k)
#'   plan <- gpc_update(plan, cut$patients, cut$events)
#' }
#' gpc_results(plan)[, c("look", "estimate", "se", "spent", "action")]
#' @export
gpc_update <- function(object, patients=NULL, events=NULL, transitions=NULL, data=NULL,
                       keep_occupation=FALSE) {
  if (!inherits(object,"gpcseq") || isTRUE(object$closed))
    stop("Provide an open gpcseq design; no update is possible after stopping.",call.=FALSE)
  k <- length(object$log)+1L; d <- object$design; K <- length(d$looks)
  if (k>K) stop("No remaining scheduled look.",call.=FALSE)
  x <- .gpc_data_snapshot(patients,events,transitions,data,d,k)
  # Bind active preset labels across named data cuts. An initially unnamed
  # empty catalog can acquire a label later; a bound label cannot change.
  for (kind in names(x$binding)) {
    before <- object$endpoint_channels[[kind]]
    now <- x$binding[[kind]]
    if (!is.null(before) && !identical(now, before))
      stop("The ", kind, " endpoint was '", before,
           "' at a previous look. Supply the same named channel, including a zero-row table if needed.", call. = FALSE)
    if (!is.null(now)) object$endpoint_channels[[kind]] <- now
  }
  object$registry <- .gpc_registry(x$p,object$registry)
  # Older saved objects lack historical rows. Preserve their numerical path,
  # but identify the first cut from which continuity can be checked.
  if (is.null(object$history_registry) && length(object$log))
    warning("Earlier event histories are unavailable in this saved plan; cumulative-history checking starts at this cut. Earlier cuts cannot be verified.", call. = FALSE)
  object$history_registry <- .gpc_history_check(object$history_registry,
    .gpc_history_records(x,d$method))
  n <- tabulate(as.integer(x$p$arm)+1L,nbins=2L)
  direction <- if (d$direction=="upper-final-two-sided") {
      if (k==K) "two-sided" else "upper"
    } else d$direction
  prior <- object$fits
  last_spent <- if (length(prior)) d$spending[tail(as.integer(names(prior)),1L)] else 0
  row <- .gpc_logrow(k,d$looks[k],n,d$spending[k],last_spent,direction)
  next_result <- tryCatch({
    row$horizon <- .gpc_horizon(d,k,x$p)
    f <- .gpc_fit_cut(x,d,k,row$horizon,keep_occupation)
    row$estimate <- f$estimate; row$natural <- f$natural
    f$direction <- direction
    fits <- prior; fits[[as.character(k)]] <- f
    V <- .gpc_cov(fits,previous=object$covariance)
    dimnames(V) <- list(names(fits),names(fits))
    j <- length(fits); vv <- V[j,j]
    if (!is.finite(vv) || vv<=0) .gpc_defer("The estimated variance is zero or nonfinite.")
    row$se <- sqrt(vv); row$z <- row$estimate/row$se
    if (!is.finite(row$z)) .gpc_defer("The standardized statistic is nonfinite.")
    R <- stats::cov2cor(V)
    prev <- unname(vapply(prior,function(i) i$critical,numeric(1)))
    dirs <- unname(vapply(fits,function(i) i$direction,character(1)))
    boundary <- .gpc_boundary(R,prev,dirs,d$spending[k],d$tol)
    row$critical <- boundary$critical; row$spent <- d$spending[k]
    row$action <- if (row$z>row$critical) "reject_upper" else if (
      direction=="two-sided" && row$z < -row$critical) "reject_lower" else if (k==K) "complete" else "continue"
    f$critical <- row$critical; f$calibration <- boundary
    fits[[as.character(k)]] <- f; list(fits=fits,V=V)
  },gpc_deferred=function(e) e)
  if (inherits(next_result,"gpc_deferred")) {
    row$reason <- conditionMessage(next_result)
    row$action <- if (k==K) "final_untested" else "deferred"
  } else {
    object$fits <- next_result$fits; object$covariance <- next_result$V
  }
  object$log[[k]] <- row
  object$closed <- k==K || row$action %in% c("reject_upper","reject_lower")
  object
}

# Construct one scheduling record; before any look its zero-row slice supplies
# the exact same field names and types to the public result accessor.
.gpc_logrow <- function(k,time,n,alpha,spent,direction) {
  data.frame(look=as.integer(k),time=as.numeric(time),horizon=NA_real_,n0=n[1],n1=n[2],
    estimate=NA_real_,natural=NA_real_,se=NA_real_,z=NA_real_,critical=NA_real_,
    alpha=alpha,spent=spent,direction=direction,action="deferred",reason="",
    stringsAsFactors=FALSE)
}

#' Inspect decisions, patient influences, and covariance across looks
#'
#' @description Extract the monitoring record, the fitted influence functions for a particular test, or the covariance matrix of the unstandardized test-scale estimates at completed looks.
#' @param object For results and covariance, a \code{gpcseq} monitoring object. For \code{gpc_influence}, either a monitoring object or a standalone \code{gpcseq_fit} returned by \code{\link{gpc_fit}}.
#' @details \code{gpc_results} distinguishes cumulative error scheduled (\code{alpha}) from the amount used (\code{spent}). A deferred look does not change the latter. When available, a deferred look may retain descriptive estimates or a standard error; take the decision from \code{action}. In a hybrid design the final level controls the total probability of crossing either direction; it need not allocate half of the overall error to each direction.
#'
#' For looks \eqn{j} and \eqn{k}, the covariance is the sum, over arms, of products of same-patient influences at the two looks divided by the product of the enrolled arm sizes. Only patients present at both looks enter that sum. The diagonal gives the variance of the corresponding estimate. Use the matrix on its returned scale: log for win odds and ratio, time units for net win time. Covariance refers to repeated looks for one specified method, not to comparisons between different methods.
#' @returns \code{gpc_results} returns a data frame with one row per processed scheduled look (empty before the first update). Its columns are:
#' \itemize{
#' \item \code{look}, \code{time}, \code{horizon}, \code{n0}, \code{n1}: scheduled index, calendar time, follow-up horizon (\code{NA} for shared-window), control and treatment enrollment.
#' \item \code{estimate}, \code{natural}, \code{se}, \code{z}, \code{critical}: test-scale estimate, original-scale estimate, test-scale standard error, standardized statistic, and boundary. Positive and negative boundaries have magnitude \code{critical} for a two-sided test.
#' \item \code{alpha}, \code{spent}, \code{direction}, \code{action}, \code{reason}: planned and used cumulative spending, sidedness at this look, decision, and any deferral explanation. Actions are \code{"continue"}, \code{"reject_upper"}, \code{"reject_lower"}, \code{"complete"} (tested final without rejection), \code{"deferred"} and \code{"final_untested"}.
#' }
#' \code{gpc_influence} returns a data frame with columns \code{id}, \code{arm}, and \code{phi}. \code{phi} is the patient's contribution on the test scale before division by the enrolled arm size. For IPCW this includes the uncertainty due to estimating the weights.
#'
#' \code{gpc_covariance} returns a square matrix for completed tests, excluding skipped looks. It is a zero-by-zero matrix before any test is performed.
#' @seealso \code{\link{gpc_update}}, \code{\link{gpc_fit}}.
#' @examples
#' x <- gpc_design(c(1, 2), c(0.01, 0.05))
#' for (k in 1:2) {
#'   cut <- gpc_example(k)
#'   x <- gpc_update(x, cut$patients, cut$events)
#' }
#' gpc_results(x)[, c("look", "spent", "action", "reason")]
#' gpc_influence(x, 2)
#' gpc_covariance(x)
#' @export
gpc_results <- function(object) {
  if (!inherits(object,"gpcseq")) stop("Expected a gpcseq object.",call.=FALSE)
  if (!length(object$log)) return(.gpc_logrow(1L,0,c(0L,0L),0,0,"")[FALSE,,drop=FALSE])
  do.call(rbind,object$log)
}

#' @rdname gpc_results
#' @param look A scheduled look number. With a monitoring object the default is the most recent tested look; requesting a deferred or absent test raises an error. A standalone fit contains one look and its contributions are returned directly.
#' @export
gpc_influence <- function(object,look=NULL) {
  if (inherits(object,"gpcseq_fit")) return(object$subjects)
  if (!inherits(object,"gpcseq")) stop("Expected a fitted look or monitoring object.",call.=FALSE)
  if (is.null(look)) look <- tail(names(object$fits),1L)
  if (length(look)!=1L || !as.character(look) %in% names(object$fits))
    stop("No completed test at this look.",call.=FALSE)
  object$fits[[as.character(look)]]$subjects
}

#' @rdname gpc_results
#' @export
gpc_covariance <- function(object) {
  if (!inherits(object,"gpcseq")) stop("Expected a gpcseq object.",call.=FALSE)
  object$covariance
}

#' Print a monitoring plan or an estimate at one look
#'
#' @description Print the planned method and effect. A monitoring object also shows the number of processed looks and the decision record, if any. A standalone fit shows the scheduled look and estimate on the test scale.
#' @param x A \code{gpcseq} monitoring object or a \code{gpcseq_fit} from a single look.
#' @param ... Additional arguments accepted for compatibility with \code{print}; currently unused.
#' @returns The input object, invisibly. For numerical output use \code{\link{gpc_results}}, \code{\link{gpc_influence}}, \code{\link{gpc_covariance}}, or the components of a standalone fit.
#' @examples
#' trial <- gpc_design(2, 0.05, "ipcw", horizons = 2)
#' print(trial)
#' x <- gpc_example(2)
#' print(gpc_fit(trial, 1, x$patients, x$events))
#' @export
print.gpcseq <- function(x,...) {
  cat("Group sequential ",x$design$method," / ",x$design$effect,": ",
      length(x$log),"/",length(x$design$looks)," looks\n",sep="")
  if (length(x$log)) print(gpc_results(x),row.names=FALSE)
  invisible(x)
}

#' @rdname print.gpcseq
#' @export
print.gpcseq_fit <- function(x,...) {
  cat("Look ",x$look,": ",x$method," / ",x$effect,", estimate = ",format(x$estimate),"\n",sep="")
  invisible(x)
}
