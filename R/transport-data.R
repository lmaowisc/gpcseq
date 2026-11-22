#' Export fully ascertained death and first hospitalization at one horizon
#'
#' @description Convert prepared longitudinal trial data to the five outcome
#' columns consumed by pairtransport's hierarchical_survival() comparison.
#' Patient identifiers are retained so the outcome can be joined to treatment
#' and baseline covariates from the trial. No target sample is constructed.
#' @param data A gpcseq_trial_data object returned by gpc_trial_data(). Its
#' selected event channel must contain the complete history of qualifying
#' hospitalizations through each patient's supplied follow-up.
#' @param horizon One positive finite follow-up time from entry, in the same
#' units as the patients' follow-up and admission times.
#' @details A participant must be observed through the earlier of death and
#' the horizon. Early censoring cannot be represented by this complete-outcome
#' format and is an error. A death or first hospitalization exactly at the
#' horizon is an event. Among survivors without an admission, the event flag
#' is zero and the time is the horizon. For a decedent, an observed admission
#' is retained; if none occurred, the two unused admission fields are missing.
#' With multiple recurrent-event channels, explicitly select the hospital
#' admission channel in gpc_trial_data() for this export.
#' No absent report is interpreted as an event-free interval. The caller must
#' prepare a data cut with complete ascertainment; this converter does not
#' select a calendar look or retrieve later reports.
#'
#' pairtransport compares survival, then the time of first hospitalization
#' among survivors. By default gpcseq compares admission counts and, on a
#' count tie, compares the times of the last admission. To match the transport comparison through
#' a common fully observed horizon, use gpc_hierarchy(admission_count = FALSE,
#' admission_time = "first") with no biomarker tier. A last-admission or count
#' tier need not agree. This converter does not fit a censoring model, calculate a
#' sequential score, or infer target covariates.
#' @return A data frame in the input patient order with character id and the
#' numeric columns death_time, death_event, hospitalization_time,
#' hospitalization_event, and followup. Event flags are integers (0 or 1).
#' Follow-up is truncated at the horizon.
#' @examples
#' p <- data.frame(id = c("T", "C"), arm = c(1L, 0L), entry = 0,
#'                 followup = c(2, 1), death = c(FALSE, TRUE))
#' admissions <- data.frame(id = "T", time = 0.5)
#' d <- gpc_trial_data(p, events = list(admission = admissions))
#' gpc_transport_outcome(d, horizon = 2)
#' @export
gpc_transport_outcome <- function(data, horizon) {
  if (!inherits(data, "gpcseq_trial_data") || !is.list(data) ||
      !all(c("patients", "events") %in% names(data)))
    stop("data must be a prepared gpcseq_trial_data object from gpc_trial_data().",
         call. = FALSE)
  if (!is.numeric(horizon) || !is.null(dim(horizon)) || length(horizon) != 1L ||
      !is.finite(horizon) || horizon <= 0)
    stop("horizon must be one positive finite follow-up time.", call. = FALSE)
  if (isTRUE(data$selection$event_ambiguous))
    stop("Select the hospitalization event channel in gpc_trial_data() before transport export.", call. = FALSE)
  if (!is.null(data$channels) && length(data$channels$events) && is.null(data$selection$event))
    stop("The hospitalization channel was omitted; select its name before transport export.", call. = FALSE)
  .gpc_frame(data$events, c("id", "time", "type"), "prepared event records")
  # Revalidate the fields used by this export, including event/ID alignment.
  # There is no new calendar cut: Inf disables only the calendar upper bound.
  adm <- data$events[data$events$type == "hospitalization", , drop = FALSE]
  if (isTRUE(data$selection$mixed_available) && "available" %in% names(adm) &&
      all(is.na(adm$available))) adm$available <- NULL
  normalized <- .gpc_snapshot(data$patients, adm, NULL,
                              list(looks = Inf, method = "shared"), 1L)
  p <- normalized$p
  e <- normalized$e
  .gpc_bad(!p$death & p$followup < horizon, p$id,
           "Incomplete follow-up: observation must extend to death or the horizon")
  died <- p$death & p$followup <= horizon
  last <- pmin(p$followup, horizon)
  h <- e[e$type == "hospitalization" & e$time <= horizon,
         c("id", "time"), drop = FALSE]
  h <- h[order(h$time, method = "radix"), , drop = FALSE]
  first <- match(p$id, h$id)
  observed <- !is.na(first)
  ht <- rep.int(as.numeric(horizon), nrow(p))
  he <- integer(nrow(p))
  ht[observed] <- h$time[first[observed]]
  he[observed] <- 1L
  ht[died & !observed] <- NA_real_
  he[died & !observed] <- NA_integer_
  data.frame(id = p$id, death_time = last, death_event = as.integer(died),
             hospitalization_time = ht, hospitalization_event = he,
             followup = last, stringsAsFactors = FALSE)
}
