# Resolve the public input once, then use exactly the native validation path.
.gpc_data_snapshot <- function(patients, events, transitions, data, d, k) {
  if (!is.null(data) && !is.null(patients))
    stop("Supply data or patients, not both.", call. = FALSE)
  if (is.null(data)) data <- patients
  prepared <- inherits(data, "gpcseq_trial_data")
  channels <- NULL
  binding <- list()
  rule <- d$hierarchy
  named_rule <- !is.null(rule$endpoints) || !is.null(rule$compare)
  alias_callback <- !is.null(rule$compare) || any(vapply(rule$endpoints, function(f)
    is.null(.pair_endpoint_spec(f)), logical(1L)))
  needs_event <- is.null(rule) || isTRUE(rule$admission_count) || rule$admission_time != "none"
  needs_marker <- !is.null(d$biomarker_margin) || !is.null(rule$biomarker_margin) ||
    !is.null(rule$biomarker_threshold)
  if (prepared) {
    if (!is.null(events) || !is.null(transitions))
      stop("A gpc_trial_data cut contains its own histories; omit separate events and transitions.", call. = FALSE)
    patients <- data$patients
    events <- data$events
    transitions <- data$transitions
    channels <- data$channels
    selection <- data$selection
    if (d$method != "aj" && !named_rule) {
      if (needs_event && isTRUE(selection$event_ambiguous))
        stop("The preset compares a recurrent event. Select event = 'channel name' in gpc_trial_data(), or specify named endpoints in the design.", call. = FALSE)
      if (needs_marker && isTRUE(selection$measurement_ambiguous))
        stop("The preset compares a biomarker. Select measurement = 'channel name' in gpc_trial_data(), or specify named endpoints in the design.", call. = FALSE)
      if (!is.null(channels)) {
        if (needs_event && is.null(selection$event) && length(channels$events))
          stop("The active preset needs a recurrent event channel; choose it explicitly.", call. = FALSE)
        if (needs_marker && is.null(selection$measurement) && length(channels$measurements))
          stop("The active preset needs a biomarker channel; choose it explicitly.", call. = FALSE)
        # A declared zero-row channel retains its name. Do not infer a
        # clinical name from a cut that has no named channel of this kind.
        if (needs_event) binding["event"] <- list(selection$event)
        if (needs_marker) binding["measurement"] <- list(selection$measurement)
      }
    }
    if (d$method != "aj" && alias_callback && !is.null(channels)) {
      # Undeclared callbacks may use the legacy shortcuts as well as
      # explicit names. Keep the meaning of those shortcuts stable.
      binding["event"] <- list(selection$event)
      binding["measurement"] <- list(selection$measurement)
    }
    if (!is.null(channels)) {
      .gtd_validate_look(channels, patients, d$looks[k])
      attr(channels, "selection") <- selection
      if (isTRUE(selection$mixed_available) && "available" %in% names(events)) {
        # Some channels carry calendar timestamps and others do not. Check
        # every reported date; no imputation for the remaining records.
        ok <- !is.na(events$available)
        key <- match(events$id[ok], patients$id)
        .gpc_bad(!is.finite(events$available[ok]) |
                   events$available[ok] < patients$entry[key] + events$time[ok] - 1e-12 |
                   events$available[ok] > d$looks[k], events$id[ok],
                 "Invalid event availability at this look")
        events$available <- NULL
      }
    }
  } else {
    if (!is.data.frame(data))
      stop("Supply a patient data frame or a cut prepared by gpc_trial_data().", call. = FALSE)
    patients <- data
  }
  if (is.null(channels) && d$method != "aj") {
    if ((!named_rule && needs_event) || alias_callback) binding$event <- "hospitalization"
    if ((!named_rule && needs_marker) || alias_callback) binding$measurement <- "biomarker"
  }
  x <- .gpc_snapshot(patients, events, transitions, d, k)
  if (d$method != "aj") .pair_check_endpoint_channels(d$hierarchy, channels)
  x$channels <- channels
  x$binding <- binding
  x
}

.gpc_frame <- function(x, cols, name) {
  if (!is.data.frame(x) || anyDuplicated(names(x)))
    stop(name," must be a data frame with unique column names.",call.=FALSE)
  absent <- setdiff(cols,names(x))
  if (length(absent)) stop(name," lacks: ",paste(absent,collapse=", "),call.=FALSE)
}
.gpc_bad <- function(condition,id,why) {
  w <- which(condition | is.na(condition))
  if (length(w)) stop(why," [ID: ", paste(head(unique(id[w]),4L),collapse=", "),"]",call.=FALSE)
}
.gpc_finite <- function(v,field) {
  if (!is.numeric(v) || any(!is.finite(v))) stop(field," must contain finite numbers.",call.=FALSE)
}

# Mark repeated (patient, time) keys after ID and time validation. A stable
# radix sort puts equal numeric times together within person. Scatter the
# marks back to input order so the diagnostic retains the first repeated ID.
.gpc_repeated_measurement <- function(id, time) {
  if (is.object(id) || is.object(time) || !is.null(dim(time)) || !is.null(dim(id)))
    return(duplicated(data.frame(id = id, time = time)))
  n <- length(time)
  if (n < 2L) return(logical(n))
  person <- match(id, id)
  order_key <- order(person, time, method = "radix")
  if (is.integer(order_key) && is.integer(person))
    return(.Call(gpc_mark_measurement_repeats, person, as.double(time), order_key))
  # Keep the R indexing path for long vectors.
  repeated <- logical(n)
  previous <- order_key[-n]
  current <- order_key[-1L]
  repeated[current] <- person[current] == person[previous] & time[current] == time[previous]
  repeated
}

.gpc_snapshot <- function(p, events, transitions, d, k) {
  .gpc_frame(p,c("id","arm","entry","followup","death"),"patients")
  if (!is.character(p$id) || anyNA(p$id) || any(!nzchar(p$id)) || anyDuplicated(p$id))
    stop("patients$id must contain unique nonempty character identifiers.",call.=FALSE)
  if (!is.numeric(p$arm) || anyNA(p$arm) || any(!p$arm %in% 0:1))
    stop("patients$arm must be numeric 0 (control) or 1 (treatment).",call.=FALSE)
  .gpc_finite(p$entry,"patients$entry"); .gpc_finite(p$followup,"patients$followup")
  if (!is.logical(p$death) || anyNA(p$death)) stop("patients$death must be TRUE or FALSE.",call.=FALSE)
  .gpc_bad(p$entry<0 | p$entry>d$looks[k],p$id,"Patient is outside this calendar cut")
  .gpc_bad(p$followup<0 | p$followup>d$looks[k]-p$entry+1e-12,p$id,
           "followup must be nonnegative and cannot extend beyond the look")
  if (is.null(events)) events <- data.frame(id=character(),time=numeric(),type=character(),value=numeric())
  if (is.null(transitions)) transitions <- data.frame(id=character(),time=numeric(),from=character(),to=character())
  validate_timing <- function(x,extra,name) {
    .gpc_frame(x,c("id","time",extra),name)
    if (!is.character(x$id) || anyNA(x$id) || any(!x$id %in% p$id))
      stop(name," contains a missing or unenrolled ID.",call.=FALSE)
    .gpc_finite(x$time,paste0(name,"$time"))
    j <- match(x$id,p$id)
    .gpc_bad(x$time<0 | x$time>p$followup[j]+1e-12,x$id,
             paste(name,"time outside observed follow-up"))
    if ("available" %in% names(x)) {
      .gpc_finite(x$available,paste0(name,"$available"))
      .gpc_bad(x$available < p$entry[j]+x$time-1e-12 | x$available>d$looks[k],x$id,
               paste(name,"availability must follow occurrence and precede or equal the look"))
    }
  }
  prepared <- NULL
  if (d$method=="aj") {
    .gpc_frame(p,"state0","patients for AJ")
    if (!is.character(p$state0) || anyNA(p$state0) || any(!p$state0 %in% d$state_order))
      stop("patients$state0 must contain declared state labels.",call.=FALSE)
    validate_timing(transitions,c("from","to"),"transitions")
    prepared <- .aj_events(p,transitions,d$state_order,d$absorbing,d$death_state)
  } else {
    validate_timing(events,"type","events")
    if (!is.character(events$type) || anyNA(events$type) ||
        any(!events$type %in% c("hospitalization","biomarker")))
      stop("events$type must be hospitalization or biomarker.",call.=FALSE)
    if (!"value" %in% names(events)) events$value <- rep(NA_real_,nrow(events))
    if (is.logical(events$value) && all(is.na(events$value))) events$value <- as.numeric(events$value)
    if (!is.numeric(events$value)) stop("events$value must be numeric.",call.=FALSE)
    bio <- events$type=="biomarker"
    if (!is.null(d$biomarker_margin) || !is.null(d$hierarchy$biomarker_margin) ||
        !is.null(d$hierarchy$biomarker_threshold) || !is.null(d$hierarchy$compare) ||
        !is.null(d$hierarchy$endpoints)) {
      .gpc_bad(bio & is.infinite(events$value),events$id,"Infinite biomarker value")
      ids <- events$id[bio]
      dup <- .gpc_repeated_measurement(ids, events$time[bio])
      if (any(dup)) stop("Resolve biomarker readings at the same time for ID ",
                                  ids[which(dup)[1L]],".",call.=FALSE)
    }
  }
  list(p=p,e=events,tr=transitions,aj_events=prepared)
}

.gpc_registry <- function(p,old) {
  fields <- intersect(c("id","entry","arm","followup","death","state0"),names(p))
  if (!is.null(old)) {
    idx <- match(old$id,p$id)
    if (anyNA(idx)) stop("Each cut must retain previously enrolled patients: ",
                         paste(head(old$id[is.na(idx)],4L),collapse=", "),call.=FALSE)
    now <- p[idx,,drop=FALSE]
    .gpc_bad(now$entry!=old$entry | now$arm!=old$arm,old$id,"Entry or arm changed")
    if ("state0" %in% names(old) && "state0" %in% names(now))
      .gpc_bad(now$state0!=old$state0,old$id,"Baseline state changed")
    .gpc_bad(now$followup<old$followup-1e-12 |
               (old$death & (!now$death | abs(now$followup-old$followup)>1e-12)),old$id,
             "Previously observed follow-up or death was reversed")
  }
  p[,fields,drop=FALSE]
}
