#' Describe three restricted win--loss effects using one IPCW fit
#'
#' @description Fit the weighted win and loss probabilities at one calendar
#'   look and report net benefit, win odds and win ratio together. This is a
#'   descriptive calculation: it does not change a monitoring object, spend
#'   error, or test several effects simultaneously.
#' @param object An IPCW trial plan returned by `gpc_design()`. Its horizon,
#'   maturity rule, clinical comparison, biomarker margin and support checks
#'   are used even when its prespecified monitoring effect is only one of
#'   the three descriptive effects.
#' @param look Integer index of the scheduled calendar look, starting at 1.
#' @param data Optional replacement for \code{patients}: a cut from \code{\link{gpc_trial_data}} containing its own histories, or a native patient data frame with separate event or transition arguments. Supply only one of \code{patients} and \code{data}.
#' @param patients One row per enrolled patient as described in `gpc_fit()`:
#'   `id`, `arm`, `entry`, `followup` and `death`. Use only records known
#'   by this calendar cut, including complete histories through follow-up.
#' @param events Event records available at the look, in the same format as
#'   `gpc_fit()`, or `NULL` for none. The clinical hierarchy set by the plan
#'   resolves comparisons after survival, using only histories through the cutoff.
#' @return A list with `look`, calendar `time`, follow-up `horizon`,
#'   `summary`, `influences` and `diagnostics`. The three rows of `summary`
#'   are `net_benefit`, `win_odds` and `win_ratio`. Each contains `estimate`
#'   (on the test scale), `natural` (on the original scale), `se` (on the
#'   test scale), `valid` and `reason`. For either ratio, the tested scale
#'   is log. An estimable constant can have zero standard error.
#'   `influences` contains `id`, `arm` and one test-scale influence column per
#'   effect, in patient order. `diagnostics` includes arm-specific censoring
#'   survival and risk counts at the horizon, and the common pair summaries.
#' @details Clinical comparisons are evaluated once and only one reverse
#'   Kaplan--Meier curve per arm is fit. Net benefit is the weighted win
#'   probability minus the weighted loss probability. Win odds follows by
#'   transformation of net benefit; win ratio uses the two component
#'   probabilities separately. All influences account for the randomness
#'   in both estimated censoring distributions. The variance is the sum
#'   of squared patient influences divided by the square of the respective
#'   arm's enrollment.
#'
#' An undefined transformation is returned with `valid = FALSE`, `NA`
#' estimates, standard error and influences, and an informative reason;
#' it does not invalidate other effects. If the common fit cannot be
#' estimated, a `gpc_deferred` condition is raised instead. Invalid
#' inputs are errors. The results are not multiplicity adjusted.
#' @examples
#' cut <- gpc_example(2)
#' design <- gpc_design(2, .05, "ipcw", horizons = 2)
#' joint <- gpc_ipcw_effects(design, 1, cut$patients, cut$events)
#' joint$summary
#' joint$influences
#' @seealso gpc_design, gpc_fit, gpc_influence, gpc_update
#' @export
gpc_ipcw_effects <- function(object, look, patients = NULL, events = NULL, data = NULL) {
  d <- .gpc_check_look(object, look)
  if (d$method != "ipcw")
    stop("gpc_ipcw_effects requires a design with method = 'ipcw'.", call. = FALSE)
  x <- .gpc_data_snapshot(patients, events, NULL, data, d, look)
  tau <- .gpc_horizon(d, look, x$p)
  .gpc_check_support(x$p, d, tau)
  core <- .pair_core(x$p, x$e, "ipcw", tau, d$biomarker_margin, d$min_risk, d$min_g, hierarchy = d$hierarchy, channels = x$channels)
  effects <- c("net_benefit", "win_odds", "win_ratio")
  sizes <- lengths(core$arms)
  denom <- ifelse(x$p$arm == 1, sizes[1L], sizes[2L])
  ans <- data.frame(effect = effects, estimate = NA_real_, natural = NA_real_,
                    se = NA_real_, valid = FALSE, reason = NA_character_)
  inf <- x$p[, c("id", "arm"), drop = FALSE]
  base_if <- NULL
  for (j in seq_along(effects)) {
    which_effect <- effects[j]
    f <- tryCatch(.pair_effect(core, x$p, which_effect, base_if),
                  gpc_deferred = function(e) e)
    if (inherits(f, "gpc_deferred")) {
      ans$reason[j] <- conditionMessage(f)
      inf[[which_effect]] <- rep(NA_real_, nrow(inf))
    } else {
      if (j == 1L) base_if <- f$influence
      ans$estimate[j] <- f$estimate
      ans$natural[j] <- f$natural
      ans$se[j] <- sqrt(sum((f$influence / denom)^2))
      ans$valid[j] <- TRUE
      inf[[which_effect]] <- f$influence
    }
  }
  list(look = look, time = d$looks[look], horizon = tau,
       summary = ans, influences = inf, diagnostics = core$diagnostics)
}
