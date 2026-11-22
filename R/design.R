#' Plan calendar looks for one pairwise treatment effect
#'
#' @description Prespecify the comparison, scheduled looks, follow-up horizons if needed, and error spending. Calendar time determines who has enrolled. Follow-up time determines which outcomes enter a pairwise comparison.
#' @param looks Positive, finite, strictly increasing numeric calendar times; the last is the final look. Entry times use the same origin and units.
#' @param spending Cumulative error spending: \code{NULL} or \code{"obrien-fleming"} uses the Lan--DeMets O'Brien--Fleming spending function; \code{"pocock"} uses the Pocock function; \code{"power"} uses Kim--DeMets power spending; and \code{"hsd"} uses Hwang--Shih--DeCani spending. \code{"linear"} is power spending with exponent one. Alternatively, supply one cumulative probability per look, strictly increasing between zero and one. Its last value is the overall level. For example, \code{c(0.01, 0.05)} spends 0.01 initially and an additional 0.04 at the final look if the first test was performed.
#' @param method One of \code{"shared"} (default), \code{"ipcw"}, or \code{"aj"}.
#' @param effect Character scalar. For \code{"ipcw"}, choose \code{"net_benefit"} (default), \code{"win_odds"}, or \code{"win_ratio"}. The only choices for \code{"shared"} and \code{"aj"} are \code{"net_benefit"} and \code{"net_win_time"}, respectively, selected automatically when \code{NULL}.
#' @param horizons One positive follow-up time for every look, or a numeric vector of positive times in scheduled order. Times are measured from each patient's entry. Select exactly one of \code{horizons} and \code{maturity} for a restricted effect; omit both for \code{"shared"}.
#' @param maturity A named list with \code{p} strictly between 0 and 1 and positive \code{cap}. Select the longest horizon up to \code{cap} that at least proportion \code{p} of the enrollees in each arm had the opportunity to reach.
#' @param direction \code{"two-sided"} (default), \code{"upper"} for treatment benefit only, or \code{"upper-final-two-sided"} for upper-only interim looks and a two-sided scheduled final look. Set \code{spending} on the intended overall scale.
#' @param min_enrolled Minimum number of enrollees in each arm at a tested look. Positive integer; default 1.
#' @param min_risk For IPCW, minimum number in each arm with observed follow-up extending to the horizon. Positive integer; default 1.
#' @param min_g For IPCW, lower bound on each arm's fitted probability of observation through the horizon. Number in \eqn{[0,1]}; default 0. A zero fitted probability is always unsupported.
#' @param biomarker_margin For win--loss comparisons, an optional nonnegative minimum difference for the last observed biomarker; lower is better. A difference equal to the margin is a tie. With \code{NULL}, no biomarker is added by this argument.
#' @param hierarchy For shared-window and IPCW win--loss comparisons, a rule from \code{\link{gpc_hierarchy}}. The default favors longer survival, fewer admissions, then later last admission. Specify the rule before the first look; it is separate from the observation and weighting method. Omit for Aalen--Johansen.
#' @param state_order For \code{"aj"}, character vector of distinct state labels in clinical order, best first. At least two states are required.
#' @param absorbing For \code{"aj"}, character vector of states known to have no exits. Defaults to none. A state with no observed exits is not automatically absorbing.
#' @param death_state For \code{"aj"}, name of the absorbing state that denotes death. Required if the records include an observed death; otherwise may be \code{NULL}.
#' @param tol Absolute probability tolerance for boundary calibration (default \code{1e-7}). For named rules calibration tightens this bound at each look if the spending increment requires it. For numeric schedules \code{tol} must be smaller than the smallest increment.
#' @param level Overall rejection probability for a named spending rule. Defaults to 0.025 for \code{direction = "upper"} and 0.05 otherwise. Omit for a numeric schedule.
#' @param spending_fractions Planned fractions at which to evaluate a named spending function, one per look, strictly increasing in (0, 1] and ending at 1. Defaults to \code{looks / looks[length(looks)]}, the fractions of planned calendar duration from time zero. Supply protocol-specified fractions to use a different scale, such as anticipated information. Omit for a numeric schedule.
#' @param spending_param Optional shape parameter for \code{"power"} (positive exponent, default 1) or \code{"hsd"} (between -40 and 40, default -4; zero gives linear spending). Omit for all other rules.
#' @details For planned fraction \eqn{t} and total level \eqn{\alpha}, the O'Brien--Fleming spending function is \eqn{2\Phi\{\Phi^{-1}(\alpha/2)/\sqrt{t}\}} and the Pocock function is \eqn{\alpha\log\{1+(e-1)t\}}. Kim--DeMets spending is \eqn{\alpha t^{\rho}} for exponent \eqn{\rho>0}. Hwang--Shih--DeCani spending is \eqn{\alpha(1-e^{-\gamma t})/(1-e^{-\gamma})} for \eqn{\gamma\ne0} and \eqn{\alpha t} for \eqn{\gamma=0}. The same function at a given level applies for all directions. These are cumulative first-crossing probabilities; the stopping boundaries are calibrated separately using the look correlations and the prespecified crossing directions. Calendar fractions do not estimate the statistical information. The first test uses a normal tail quantile even when its continuation probability rounds to one. Named schedules require a positive representable first spend and numerically distinguishable continuation probabilities at later looks.
#'
#' \code{"shared"} compares each treatment--control pair at its shorter observed follow-up time. Its signed mean tests equality of the win and loss probabilities throughout follow-up; its magnitude depends on censoring and is not net benefit at a common horizon.
#'
#' \code{"ipcw"} estimates a common-horizon effect using inverse probabilities of observation through the first death in a pair or the horizon, whichever comes first. \code{"aj"} integrates differences in arm-specific state favorability using Aalen--Johansen state occupation estimates.
#'
#' With \code{maturity}, sort the potential follow-up times (look minus entry) in decreasing order within arm; take the observation at rank \eqn{\lceil p n\rceil} in each and cap the smaller at \code{cap}. Death and dropout do not shorten potential follow-up. Ties are counted without interpolation. New enrollment can shorten a later horizon. For IPCW, inspect the number still observed in each arm and the estimated censoring survival at that horizon. Few patients may remain under observation even when many enrolled early enough to reach the horizon.
#'
#' The procedure monitors one prespecified effect. Stopping for the first of several endpoints to cross requires additional joint calibration. With changing horizons the intended null is that the effect is zero throughout follow-up. The monitoring approximation requires adequate support, independent patients, and within-arm independence of outcome histories from entry and censoring. When few patients or events contribute at an early look, the Gaussian approximation can be inaccurate even if all support checks pass. Before choosing the looks and support thresholds, simulate the proposed trial under the null with its planned accrual, event rates and censoring. Estimate the probability of crossing any boundary over the trial and the frequency of skipped looks.
#' @returns An open \code{gpcseq} object. Its \code{design} contains the cumulative \code{spending}, \code{spending_rule} (\code{"custom"} for numeric input), \code{level}, the \code{spending_fractions} (\code{NULL} for numeric input), and \code{spending_param} when applicable. Pass the object to \code{\link{gpc_update}} to process the first look.
#' @seealso \code{\link{gpc_fit}}, \code{\link{gpc_update}}, \code{\link{gpc_results}}.
#' @examples
#' plan <- gpc_design(c(12, 24, 36))  # Default O'Brien--Fleming spending
#' plan$design$spending
#' plan <- gpc_design(c(1, 2), c(0.01, 0.05), "ipcw", horizons = c(1, 2))
#' early <- gpc_design(c(18, 30, 48), c(0.005, 0.015, 0.05),
#'                     "ipcw", maturity = list(p = 0.25, cap = 24))
#' @importFrom utils head tail
#' @importFrom stats setNames
#' @export
gpc_design <- function(looks, spending = NULL, method = "shared", effect = NULL,
                       horizons = NULL, maturity = NULL, direction = "two-sided",
                       min_enrolled = 1L, min_risk = 1L, min_g = 0,
                       biomarker_margin = NULL, state_order = NULL,
                       absorbing = character(), death_state = NULL, tol = 1e-7,
                       hierarchy = NULL, level = NULL, spending_fractions = NULL,
                       spending_param = NULL) {
  numeric_vector <- function(x) is.numeric(x) && length(x) > 0L && all(is.finite(x))
  choice <- function(x, choices, field) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !(x %in% choices))
      stop(field, " must be one of: ", paste(choices, collapse = ", "), call. = FALSE)
    x
  }
  if (!numeric_vector(looks) || any(looks <= 0) || any(diff(looks) <= 0))
    stop("looks must be finite, positive and strictly increasing.", call. = FALSE)
  method <- choice(method, c("shared", "ipcw", "aj"), "method")
  direction <- choice(direction, c("two-sided", "upper", "upper-final-two-sided"), "direction")
  named_spending <- is.null(spending) || is.character(spending)
  if (named_spending) {
    spending_rule <- if (is.null(spending)) "obrien-fleming" else
      choice(spending, c("obrien-fleming", "pocock", "power", "hsd", "linear"), "spending")
    if (spending_rule %in% c("power", "hsd")) {
      if (is.null(spending_param)) spending_param <- if (spending_rule == "power") 1 else -4
      if (!numeric_vector(spending_param) || !is.null(dim(spending_param)) ||
          length(spending_param) != 1L)
        stop("spending_param must be one finite number.", call. = FALSE)
      if (spending_rule == "power" && spending_param <= 0)
        stop("Kim--DeMets power spending needs a positive spending_param.", call. = FALSE)
      if (spending_rule == "hsd" && abs(spending_param) > 40)
        stop("Hwang--Shih--DeCani spending_param must be between -40 and 40.", call. = FALSE)
    } else if (!is.null(spending_param))
      stop("spending_param is only used with 'power' or 'hsd' spending.", call. = FALSE)
    if (spending_rule == "linear") {
      spending_rule <- "power"
      spending_param <- 1
    }
    if (is.null(level)) level <- if (direction == "upper") 0.025 else 0.05
    if (!numeric_vector(level) || !is.null(dim(level)) || length(level) != 1L ||
        level <= 0 || level >= 1)
      stop("level must be one finite number strictly between zero and one.", call. = FALSE)
    if (is.null(spending_fractions)) spending_fractions <- looks / looks[length(looks)]
    if (!numeric_vector(spending_fractions) || !is.null(dim(spending_fractions)) ||
        length(spending_fractions) != length(looks) ||
        any(spending_fractions <= 0) || any(diff(spending_fractions) <= 0) ||
        spending_fractions[length(looks)] != 1)
      stop("spending_fractions must be positive and strictly increasing, one per look, ending at 1.",
           call. = FALSE)
    # expm1(x)/x tends to 1 at zero, including when gamma*t underflows.
    exprel <- function(x) {
      out <- rep(1, length(x))
      nonzero <- x != 0
      out[nonzero] <- expm1(x[nonzero]) / x[nonzero]
      out
    }
    spending <- switch(spending_rule,
      "obrien-fleming" = 2 * stats::pnorm(stats::qnorm(level / 2) / sqrt(spending_fractions)),
      "pocock" = level * log1p(expm1(1) * spending_fractions),
      "power" = level * spending_fractions ^ spending_param,
      "hsd" = if (spending_param == 0) level * spending_fractions else
                level * spending_fractions * (exprel(-spending_param * spending_fractions) /
                                                exprel(-spending_param)))
    spending[length(looks)] <- level
    # The first marginal quantile uses alpha directly. Later joint boundaries
    # require numerically distinct successive continuation probabilities.
    unresolved <- which(c(spending[1L] <= 0, diff(1 - spending) >= 0))
    if (length(unresolved))
      stop("The planned spending at look ", unresolved[1L],
           " cannot be resolved at numerical precision. Revise the planned looks or spending rule.",
           call. = FALSE)
  } else {
    if (!is.null(level) || !is.null(spending_fractions) || !is.null(spending_param))
      stop("For numeric spending, omit level, spending_fractions and spending_param; the vector specifies the entire schedule.",
           call. = FALSE)
    if (!numeric_vector(spending) || length(spending) != length(looks) ||
        any(spending <= 0 | spending >= 1) || any(diff(spending) <= 0))
      stop("spending must increase strictly between zero and one, once per look.", call. = FALSE)
    spending_rule <- "custom"
    level <- spending[length(spending)]
  }
  allowed <- switch(method, shared = "net_benefit", ipcw = c("net_benefit","win_odds","win_ratio"),
                    aj = "net_win_time")
  if (is.null(effect)) effect <- allowed[1L]
  effect <- choice(effect, allowed, "effect")
  if (method == "shared" && (!is.null(horizons) || !is.null(maturity)))
    stop("Shared windows are pair dependent; omit horizons and maturity.", call. = FALSE)
  if (method != "shared" && is.null(horizons) == is.null(maturity))
    stop("A restricted analysis needs exactly one of horizons and maturity.", call. = FALSE)
  if (!is.null(horizons)) {
    if (!numeric_vector(horizons) || any(horizons <= 0) ||
        !(length(horizons) %in% c(1L, length(looks))))
      stop("horizons must be positive finite follow-up times: one common value, or one per look in scheduled order.",
           call. = FALSE)
    horizons <- rep_len(horizons, length(looks))
  }
  if (!is.null(maturity)) {
    if (!is.list(maturity) || !setequal(names(maturity), c("p","cap")) ||
        length(maturity) != 2L || !numeric_vector(maturity$p) ||
        length(maturity$p) != 1L || maturity$p <= 0 || maturity$p >= 1 ||
        !numeric_vector(maturity$cap) || length(maturity$cap) != 1L || maturity$cap <= 0)
      stop("maturity needs p strictly between 0 and 1 and a positive cap.", call. = FALSE)
  }
  posint <- function(x) numeric_vector(x) && length(x)==1L && x >= 1 && x == floor(x)
  if (!posint(min_enrolled) || !posint(min_risk))
    stop("min_enrolled and min_risk must be positive whole numbers.", call. = FALSE)
  if (!numeric_vector(min_g) || length(min_g)!=1L || min_g<0 || min_g>1)
    stop("min_g must be between 0 and 1.", call. = FALSE)
  if (!numeric_vector(tol) || length(tol)!=1L || tol<=0 || tol>=1)
    stop("tol must be a positive number smaller than one.", call. = FALSE)
  if (!named_spending && tol>=min(diff(c(0,spending))))
    stop("tol must be positive and smaller than the smallest spending increment.", call. = FALSE)
  if (!is.null(biomarker_margin) && (!numeric_vector(biomarker_margin) ||
      length(biomarker_margin)!=1L || biomarker_margin<0 || method=="aj"))
    stop("biomarker_margin must be nonnegative and used with a win-loss method.", call. = FALSE)
  if (!is.null(hierarchy)) {
    if (method == "aj") stop("hierarchy is for win-loss methods; use state_order for AJ.", call. = FALSE)
    if (!inherits(hierarchy, "gpcseq_hierarchy")) stop("Use gpc_hierarchy() to specify the comparison.", call. = FALSE)
    if (!is.null(biomarker_margin)) stop("Specify biomarker_margin inside hierarchy or in gpc_design(), not both.", call. = FALSE)
  }
  if (is.null(hierarchy) && method != "aj")
    hierarchy <- gpc_hierarchy(biomarker_margin = biomarker_margin)
  if (method == "aj") {
    if (!is.character(state_order) || length(state_order)<2L || anyNA(state_order) ||
        any(!nzchar(state_order)) || anyDuplicated(state_order))
      stop("state_order must name at least two unique states, best first.", call. = FALSE)
    if (!is.character(absorbing) || anyNA(absorbing) || anyDuplicated(absorbing) ||
        any(!absorbing %in% state_order))
      stop("absorbing must be distinct state names, each listed in state_order.", call. = FALSE)
    if (!is.null(death_state) && (!is.character(death_state) || length(death_state)!=1L ||
        !death_state %in% absorbing)) stop("death_state must be one declared absorbing state.", call. = FALSE)
  } else if (!is.null(state_order) || length(absorbing) || !is.null(death_state))
    stop("State definitions apply to AJ only.", call. = FALSE)
  spec <- list(looks=looks,spending=spending,method=method,effect=effect,horizons=horizons,
               maturity=maturity,direction=direction,min_enrolled=min_enrolled,
               min_risk=min_risk,min_g=min_g,biomarker_margin=biomarker_margin,
               state_order=state_order,absorbing=absorbing,death_state=death_state,tol=tol,
               hierarchy=hierarchy,spending_rule=spending_rule,level=level,
               spending_fractions=spending_fractions,spending_param=spending_param)
  structure(list(design=spec, fits=list(), log=list(), covariance=matrix(numeric(),0,0),
                 registry=NULL, closed=FALSE), class="gpcseq")
}

.gpc_defer <- function(reason) {
  stop(structure(list(message=as.character(reason),call=NULL),
                 class=c("gpc_deferred","error","condition")))
}

.gpc_horizon <- function(d,k,p) {
  if (d$method=="shared") return(NA_real_)
  if (!is.null(d$horizons)) return(d$horizons[k])
  candidate <- vapply(0:1,function(a) {
    v <- d$looks[k]-p$entry[p$arm==a]
    if (!length(v)) return(NA_real_)
    # Compare the enrolled fraction directly: e.g. 14/100 == .14, even
    # though .14 * 100 rounds just above 14.
    n <- length(v)
    m <- max(1, floor(d$maturity$p * n))
    if (m / n < d$maturity$p) m <- m + 1
    i <- n - m + 1
    sort.int(v, partial=i)[i]
  },numeric(1))
  if (anyNA(candidate)) .gpc_defer("An arm has no enrolled patients for horizon selection.")
  min(candidate,d$maturity$cap)
}
