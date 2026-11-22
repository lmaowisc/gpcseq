#' Extract arm-specific Aalen--Johansen state occupation probabilities
#'
#' @param object A standalone Aalen--Johansen fit or a monitoring object with
#'   a performed Aalen--Johansen test. Set `keep_occupation = TRUE` when
#'   calling [gpc_fit()] or [gpc_update()] for the look of interest.
#' @param look Optional scheduled look index. For a monitoring object the
#'   default is the last performed test. For a standalone fit the index, if
#'   supplied, must match that fit.
#' @param times Optional numeric vector of elapsed follow-up times from
#'   entry, between zero and the fitted horizon, inclusive. Supply a common
#'   grid to compare the arms at the same times. The order is preserved
#'   within each arm. By default each arm is reported at entry, at its
#'   observed transition times through the horizon and at the horizon.
#' @returns A data frame with columns `look`, `arm` (0 control, 1 treatment),
#'   `time` (elapsed time), `state` (in the prespecified clinical order) and
#'   `probability`. Each arm's probabilities sum to one at each time.
#' @details The curves are right-continuous. Clinical transitions exactly at
#'   a requested time are included before tied censoring. A transition at the
#'   horizon affects the reported probability there and contributes zero
#'   duration to restricted net win time. The function reads stored curves;
#'   an earlier saved fit without curves requires the original data cut to
#'   obtain them with [gpc_fit()].
#' @examples
#' cut <- gpc_example(2)
#' plan <- gpc_design(2, .05, "aj", horizons = 2,
#'   state_order = c("well", "hospitalized", "dead"),
#'   absorbing = "dead", death_state = "dead")
#' fit <- gpc_fit(plan, 1, cut$patients, transitions = cut$transitions,
#'                keep_occupation = TRUE)
#' gpc_occupation(fit, times = c(0, 1, 2))
#' @export
gpc_occupation <- function(object, look = NULL, times = NULL) {
  if (!is.null(look) && (!is.numeric(look) || length(look) != 1L ||
      !is.finite(look) || look != floor(look) || look < 1L))
    stop("look must be one positive scheduled look index.", call. = FALSE)
  if (inherits(object, "gpcseq")) {
    if (is.null(look)) look <- utils::tail(names(object$fits), 1L)
    if (length(look) != 1L || !as.character(look) %in% names(object$fits))
      stop("No completed test at this look.", call. = FALSE)
    object <- object$fits[[as.character(look)]]
  }
  if (!inherits(object, "gpcseq_fit") || object$method != "aj")
    stop("Supply an Aalen--Johansen fit or monitoring object.", call. = FALSE)
  if (!is.null(look) && as.character(look) != as.character(object$look))
    stop("The requested look does not match this fit.", call. = FALSE)
  if (is.null(object$occupation))
    stop("No curves stored at this look. Use keep_occupation = TRUE with its original data cut.",
         call. = FALSE)
  if (!is.null(times) && (!is.numeric(times) || !is.null(dim(times)) ||
      any(!is.finite(times)) || any(times < 0 | times > object$horizon)))
    stop("times must be finite follow-up times between zero and the fitted horizon.",
         call. = FALSE)
  # A common elapsed-time grid has the same state labels and times in both arms.
  # Assemble its public table once. Preserve base row naming for named grids.
  if (!is.null(times) && is.null(names(times)) && !is.object(times)) {
    states <- colnames(object$occupation[[1L]]$probability)
    n <- length(times) * length(states)
    probabilities <- lapply(object$occupation, function(curve) {
      index <- findInterval(times, curve$time)
      as.vector(t(curve$probability[index, , drop = FALSE]))
    })
    return(data.frame(look = rep.int(object$look, 2 * n),
                      arm = rep(0:1, each = n),
                      time = rep(rep(times, each = length(states)), 2L),
                      state = rep(states, 2 * length(times)),
                      probability = unlist(probabilities, use.names = FALSE),
                      stringsAsFactors = FALSE))
  }
  rows <- lapply(0:1, function(arm) {
    curve <- object$occupation[[arm + 1L]]
    at <- if (is.null(times)) curve$time else times
    P <- curve$probability[findInterval(at, curve$time), , drop = FALSE]
    nr <- length(P)
    data.frame(look = rep.int(object$look, nr), arm = rep.int(arm, nr),
               time = rep(at, each = ncol(P)), state = rep(colnames(P), nrow(P)),
               probability = as.vector(t(P)), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
