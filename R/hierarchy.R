#' Specify a survival-prioritized clinical comparison
#'
#' @description Prespecify how to compare two patients who are alive at the
#' comparison time. Death has priority under every rule: an earlier death
#' loses, including at the cutoff, and simultaneous deaths tie.
#' @param admission_count Logical. Compare the number of hospitalizations
#' through the comparison time, favoring fewer. The default is TRUE.
#' @param admission_time One of "last" (default), "first" or "none". Compare
#' the elapsed times of the indicated admissions through the cutoff, favoring
#' later. An admission-free history beats one with an admission, including an
#' admission exactly at the cutoff. If neither patient was admitted, pass to
#' the next tier.
#' @param biomarker_margin Optional minimum difference for the biomarker. Give
#' a finite nonnegative margin in its clinical units; zero compares unequal
#' last observed values. A difference equal to the margin ties. A missing
#' measurement when this tier must decide the pair defers the look.
#' @param biomarker_direction "lower" (default) or "higher", the favorable
#' direction for the biomarker. Specify a margin or threshold to include this tier.
#' @param compare Optional function of three positional arguments: the
#' treatment history, the control history and the follow-up cutoff. It
#' replaces all built-in tiers below survival. Each history contains events
#' (named time vectors) and measurements (named time/value data frames),
#' sorted and clipped at the cutoff. The admissions and biomarker aliases
#' refer to the selected convenience channels, when available. Use exact
#' channel names; an undeclared endpoint is not known to be event-free.
#' Return one number: 1 for a treatment win, -1 for a loss or 0 for a tie.
#' Use the same deterministic clinical rule for both patients; swapping them
#' must reverse the sign. The function may use only these histories and the
#' cutoff. Do not combine it with nondefault built-in tier settings.
#' @param biomarker_threshold Optional finite clinical cutoff in biomarker
#' units. Instead of comparing a pairwise difference, compare whether each
#' latest value is strictly above this cutoff. A value equal to the cutoff
#' is in the at-or-below group; equal categories tie. For "lower", the
#' at-or-below group wins; "higher" reverses this preference. Specify one
#' of biomarker_margin and biomarker_threshold, or leave both NULL to omit
#' this tier.
#' @param endpoints Optional named list of lower-tier functions in clinical
#' priority order. Use \code{\link{gpc_endpoint}} for named event counts,
#' admission or other event times, and serial measurements, or provide a
#' function of the two histories and time for a derived endpoint. Each
#' function returns 1, -1, or 0. Only a tie advances to the next tier.
#' Endpoints replace all built-in lower tiers, and cannot be combined with
#' compare. Names must be distinct. An empty list uses survival alone.
#' @details The built-in order is survival, admission count, time to first
#' or last admission, then biomarker. The first unequal tier settles the
#' comparison. Only specified tiers are used. Both treatment effects and
#' their null hypotheses depend on the rule; use the same prespecified rule
#' throughout a trial. A custom function is called only for an ascertained
#' pair in which both patients are alive. IPCW and the shared-window test
#' apply the same rule; they differ in observation window and weighting.
#' @returns A gpcseq_hierarchy specification for the hierarchy argument
#' of \code{\link{gpc_design}}. It contains no patient data.
#' @examples
#' gpc_hierarchy()  # Survival, fewer admissions, then later last admission
#' gpc_hierarchy(biomarker_threshold = 10)  # Then favor a value at or below 10
#' rule <- gpc_hierarchy(admission_time = "first", biomarker_margin = 0)
#' plan <- gpc_design(2, 0.05, "shared", hierarchy = rule)
#' # Survival, then later last admission; do not compare counts.
#' gpc_hierarchy(admission_count = FALSE, admission_time = "last")
#' # Custom lower tier: fewer admissions occurring after time 1.
#' gpc_hierarchy(compare = function(treatment, control, time)
#'   sign(sum(control$admissions > 1) - sum(treatment$admissions > 1)))
#' @seealso \code{\link{gpc_design}}, \code{\link{gpc_fit}}.
#' @export
gpc_hierarchy <- function(admission_count = TRUE,
                          admission_time = c("last", "first", "none"),
                          biomarker_margin = NULL,
                          biomarker_direction = c("lower", "higher"),
                          compare = NULL, biomarker_threshold = NULL, endpoints = NULL) {
  # Use exact choices; a misspelled protocol option should not partially match.
  if (missing(admission_time)) admission_time <- "last"
  if (missing(biomarker_direction)) biomarker_direction <- "lower"
  if (!is.logical(admission_count) || length(admission_count) != 1L ||
      is.na(admission_count))
    stop("admission_count must be TRUE or FALSE.", call. = FALSE)
  if (!is.character(admission_time) || length(admission_time) != 1L ||
      !(admission_time %in% c("none", "first", "last")))
    stop("admission_time must be none, first or last.", call. = FALSE)
  if (!is.character(biomarker_direction) || length(biomarker_direction) != 1L ||
      !(biomarker_direction %in% c("lower", "higher")))
    stop("biomarker_direction must be lower or higher.", call. = FALSE)
  if (!is.null(biomarker_margin) && (!is.numeric(biomarker_margin) ||
      length(biomarker_margin) != 1L || !is.finite(biomarker_margin) || biomarker_margin < 0))
    stop("biomarker_margin must be NULL or a finite nonnegative number.", call. = FALSE)
  if (!is.null(biomarker_threshold) && (!is.numeric(biomarker_threshold) ||
      length(biomarker_threshold) != 1L || !is.finite(biomarker_threshold)))
    stop("biomarker_threshold must be NULL or a finite clinical cutoff.", call. = FALSE)
  if (!is.null(biomarker_threshold) && !is.null(biomarker_margin))
    stop("Choose either biomarker_threshold or biomarker_margin, not both.", call. = FALSE)
  if (!is.null(compare) && (!is.function(compare) || !is.null(endpoints)))
    stop("Specify a single compare function or ordered endpoints, not both.", call. = FALSE)
  if (!is.null(endpoints) && (!is.list(endpoints) || is.data.frame(endpoints) ||
      (length(endpoints) && (is.null(names(endpoints)) || anyNA(names(endpoints)) ||
       any(!nzchar(names(endpoints))) || anyDuplicated(names(endpoints)) ||
       !all(vapply(endpoints, is.function, logical(1L)))))))
    stop("endpoints must be a list of functions with distinct, nonempty names.", call. = FALSE)
  if ((!is.null(compare) || !is.null(endpoints)) && (!admission_count ||
      admission_time != "last" || !is.null(biomarker_margin) || !is.null(biomarker_threshold) ||
      biomarker_direction != "lower"))
    stop("Do not combine compare or endpoints with nondefault built-in lower tiers.", call. = FALSE)
  if (is.null(biomarker_margin) && is.null(biomarker_threshold) && biomarker_direction != "lower")
    stop("Supply a biomarker margin or threshold to use a biomarker tier.", call. = FALSE)
  structure(list(admission_count = admission_count, admission_time = admission_time,
                 biomarker_margin = biomarker_margin, biomarker_direction = biomarker_direction,
                 compare = compare, biomarker_threshold = biomarker_threshold, endpoints = endpoints), class = "gpcseq_hierarchy")
}

#' Compile one clinical score provider for a validated cut
#'
#' @param patients,events Validated look-specific patient and event tables.
#' @param ti,ci Treatment and control row indices in patients.
#' @param ordered_ends Treatment observation ends in nondecreasing order;
#' IPCW caps them at the horizon. Used only to index the control histories.
#' @param rule Prespecified gpcseq_hierarchy with a clinical comparison.
#' @param restricted Whether all living, ascertained pairs are scored at the same horizon.
#' @param channels Optional named event and measurement table lists for independent endpoints.
#' @return A function of the treatment rank r in ordered_ends, its arm-local
#' row index i, a vector j of living control indices, and their corresponding
#' pair cutoffs. Returns one signed lower-tier score in each j position.
#' Ranks must be requested in ascending order (some may be skipped).
#' @details Only the built-in provider advances event counts and times once
#' across the control histories. The custom provider supplies per-pair clipped
#' lists and directly validates the return. Neither reads weights or fits a
#' censoring model; the driver handles ascertainment and first death.
#' @keywords internal
#' @noRd
.pair_score_provider <- function(patients, events, ti, ci, ordered_ends, rule, channels = NULL, restricted = FALSE) {
  if (!inherits(rule, "gpcseq_hierarchy"))
    stop("hierarchy must be obtained from gpc_hierarchy().", call. = FALSE)
  if (!is.null(rule$compare) || !is.null(rule$endpoints))
    return(.pair_tier_provider(patients, events, ti, ci, rule, channels, restricted))
  marker <- !is.null(rule$biomarker_margin) || !is.null(rule$biomarker_threshold)
  counts <- rule$admission_count
  timing <- rule$admission_time
  if (!marker && !counts && timing == "none")
    return(function(r, i, j, at) numeric(length(j)))
  ev <- .pair_events(patients, events, ti, ci, marker)
  nc <- length(ci)
  # The control sweep cannot pass a control's own last observed event.
  # Both arms use events only through min(treatment end, control end).
  need_admission <- counts || timing != "none"
  hc <- integer(nc); when <- rep(NA_real_, nc); bc <- rep(NA_real_, nc)
  ih <- ib <- 0L
  if (need_admission) h_end <- findInterval(ordered_ends, ev$ch)
  if (marker) b_end <- findInterval(ordered_ends, ev$cb)
  function(r, i, j, at) {
    if (need_admission && h_end[r] > ih) {
      z <- seq.int(ih + 1L, h_end[r]); ids <- ev$chi[z]
      if (timing != "none") {
        use <- if (timing == "first") !duplicated(ids) & is.na(when[ids]) else !duplicated(ids, fromLast = TRUE)
        when[ids[use]] <<- ev$ch[z[use]]
      }
      hc <<- hc + tabulate(ids, nc)
      ih <<- h_end[r]
    }
    if (marker && b_end[r] > ib) {
      z <- seq.int(ib + 1L, b_end[r]); ids <- ev$cbi[z]
      last <- !duplicated(ids, fromLast = TRUE)
      bc[ids[last]] <<- ev$cv[z[last]]
      ib <<- b_end[r]
    }
    ans <- numeric(length(j))
    if (need_admission) nt <- findInterval(at, ev$th[[i]])
    if (counts) ans <- sign(hc[j] - nt)
    if (timing != "none") {
      a <- which(ans == 0)
      ct <- hc[j[a]]; tt <- nt[a]
      ans[a] <- as.integer(tt == 0L) - as.integer(ct == 0L)
      both <- a[tt > 0L & ct > 0L]
      if (length(both)) {
        time_t <- if (timing == "first") ev$th[[i]][1L] else ev$th[[i]][nt[both]]
        ans[both] <- sign(time_t - when[j[both]])
      }
    }
    if (marker) {
      a <- which(ans == 0)
      if (length(a)) {
        pos <- findInterval(at[a], ev$tb[[i]])
        bt <- c(NA_real_, ev$tv[[i]])[1L + pos]
        bcontrol <- bc[j[a]]
        absent <- is.na(bt) | is.na(bcontrol)
        if (any(absent)) {
          s <- a[which(absent)[1L]]
          .gpc_defer(sprintf("Missing biomarker when needed for patients %s and %s at follow-up %g.",
                             patients$id[ti[i]], patients$id[ci[j[s]]], at[s]))
        }
        if (!is.null(rule$biomarker_threshold)) {
          delta <- as.integer(bcontrol > rule$biomarker_threshold) -
            as.integer(bt > rule$biomarker_threshold)
          margin <- 0
        } else {
          delta <- bcontrol - bt
          margin <- rule$biomarker_margin
        }
        if (rule$biomarker_direction == "higher") delta <- -delta
        ans[a] <- (delta > margin) - (delta < -margin)
      }
    }
    ans
  }
}

#' Define a comparison tier from a named endpoint
#'
#' @param channel One exact name in the named events (for count, first and
#' last) or measurements (for value) supplied to gpc_trial_data. On native
#' patient/event input use hospitalization and biomarker, respectively.
#' @param type "count" (default) compares the number of events; "first" and
#' "last" compare the elapsed time of the indicated event; "value" compares
#' the last measured value. Include records at the pair cutoff.
#' @param direction "lower" or "higher". The defaults favor fewer events,
#' a later event, and a lower measurement. An event-free history is later
#' than any observed event, including an event exactly at the cutoff. Two
#' event-free histories tie for timing. A time margin applies only when
#' both have an event.
#' @param threshold Optional finite clinical cutoff for a value tier. A value
#' strictly above the cutoff belongs to the upper category; at-or-below is
#' the lower category. Equal categories tie. Use margin zero with threshold.
#' @param margin Finite nonnegative pairwise margin, in the units of the
#' summary being compared. Differences at or below the margin tie.
#' @description Construct a reusable function comparing one endpoint. List
#' these functions, in priority order, in the endpoints argument of
#' \code{\link{gpc_hierarchy}}. An unknown channel is an input error.
#' A known event channel can have zero records. A deciding value tier with
#' no observed measurement for either patient defers the look.
#' @returns A signed comparison function of treatment history, control
#' history and follow-up time. The package indexes the channel specified
#' in its metadata for repeated use at a look.
#' @examples
#' rule <- gpc_hierarchy(endpoints = list(
#'   relapse = gpc_endpoint("relapse", "count"),
#'   function_score = gpc_endpoint("function_score", "value", "higher")
#' ))
#' @seealso \code{\link{gpc_hierarchy}}, \code{\link{gpc_trial_data}}.
#' @export
gpc_endpoint <- function(channel, type = c("count", "first", "last", "value"),
                         direction = NULL, threshold = NULL, margin = 0) {
  if (missing(type)) type <- "count"
  if (!is.character(channel) || length(channel) != 1L || is.na(channel) || !nzchar(channel))
    stop("channel must be one exact nonempty endpoint name.", call. = FALSE)
  if (!is.character(type) || length(type) != 1L || !(type %in% c("count", "first", "last", "value")))
    stop("type must be count, first, last or value.", call. = FALSE)
  if (is.null(direction)) direction <- if (type %in% c("first", "last")) "higher" else "lower"
  if (!is.character(direction) || length(direction) != 1L || !(direction %in% c("lower", "higher")))
    stop("direction must be lower or higher.", call. = FALSE)
  scalar <- function(x) is.numeric(x) && is.null(dim(x)) && length(x) == 1L && is.finite(x)
  if (!scalar(margin) || margin < 0) stop("margin must be a finite nonnegative number.", call. = FALSE)
  if (!is.null(threshold) && (!scalar(threshold) || type != "value" || margin != 0))
    stop("A clinical threshold requires a value tier with zero pairwise margin.", call. = FALSE)
  spec <- list(channel = channel, type = type, direction = direction, threshold = threshold, margin = margin)
  ans <- function(treatment, control, time) {
    group <- if (type == "value") "measurements" else "events"
    pick <- function(h) {
      z <- h[[group]]
      if (!channel %in% names(z)) stop("Unknown ", group, " channel: ", channel, ".", call. = FALSE)
      .pair_tier_value(z[[channel]], time, type)
    }
    .pair_tier_contrast(pick(treatment), pick(control), spec)
  }
  attr(ans, "gpcseq_endpoint") <- spec
  ans
}

#' Exact lookup of a declared history channel in a custom comparison
#' @param x Checked named list of measurement histories, event histories or history fields.
#' @param i One exact character name or an in-range integer index.
#' @param ... Not used.
#' @param exact Must be TRUE; channel names do not support partial matching.
#' @return The selected history or component, including a genuinely empty history.
#' @export
#' @noRd
`[[.gpcseq_channels` <- function(x, i, ..., exact = TRUE) {
  z <- unclass(x)
  if (!isTRUE(exact) || length(i) != 1L || is.na(i) ||
      !(if (is.character(i)) i %in% names(z) else
        is.numeric(i) && is.finite(i) && i == floor(i) && i >= 1 && i <= length(z)))
    stop("Unknown or unavailable history channel/field: ", paste(i, collapse = ", "),
         ". Use an exact declared name.", call. = FALSE)
  z[[i]]
}

#' @param name Exact history channel or field name.
#' @rdname internal-gpcseq-channels
#' @export
#' @noRd
`$.gpcseq_channels` <- function(x, name) x[[name, exact = TRUE]]

#' Retrieve the declared tier for a function created by gpc_endpoint
#' @param f A clinical comparison function.
#' @return The captured specification when the function is still the original
#'   endpoint closure, or NULL when its body must be evaluated as a callback.
#' @noRd
.pair_endpoint_spec <- function(f) {
  spec <- attr(f, "gpcseq_endpoint", exact = TRUE)
  if (is.null(spec) || !is.function(f)) return(NULL)
  frame <- environment(f)
  if (is.null(frame) || identical(frame, emptyenv()) ||
      !identical(parent.env(frame), environment(gpc_endpoint))) return(NULL)
  # The constructor retains the original closure in its own frame. Replacing
  # a function body or copying its attributes to a callback must not preserve
  # the original scoring shortcut. Saved genuine endpoints retain this frame.
  if (!identical(frame$ans, f) || !identical(frame$spec, spec) ||
      !identical(frame$channel, spec$channel) || !identical(frame$type, spec$type))
    return(NULL)
  spec
}

#' Validate named references before fitting a look
#' @param rule A clinical hierarchy; NULL represents an older preset design.
#' @param channels Named event and measurement lists, or NULL for native tables.
#' @return Invisibly NULL; an absent declared endpoint is a data error.
#' @noRd
.pair_check_endpoint_channels <- function(rule, channels) {
  fs <- if (is.null(rule$compare)) rule$endpoints else list(rule$compare)
  for (f in fs) {
    sp <- .pair_endpoint_spec(f)
    if (is.null(sp)) next
    group <- if (sp$type == "value") "measurements" else "events"
    known <- if (is.null(channels)) {
      if (group == "events") "hospitalization" else "biomarker"
    } else names(channels[[group]])
    if (!(sp$channel %in% known))
      stop("Unknown ", group, " channel: ", sp$channel, ". Declared: ",
           paste(known, collapse = ", "), ".", call. = FALSE)
  }
  invisible(NULL)
}

#' Summarize a single observed endpoint through a follow-up cutoff
#' @param z Sorted event times or a sorted time/value measurement table.
#' @param at One or more pair cutoffs (elapsed times).
#' @param type One declared endpoint type.
#' @return Numeric summaries; Inf means event-free for first/last, and NA
#' means no available deciding measurement. Cutoff inclusion is exact.
#' @noRd
.pair_tier_value <- function(z, at, type) {
  if (type == "value") {
    n <- findInterval(at, z$time)
    ans <- rep(NA_real_, length(n)); hit <- n > 0L
    ans[hit] <- z$value[n[hit]]
    return(ans)
  }
  n <- findInterval(at, z)
  if (type == "count") return(n)
  if (type == "first") return(ifelse(n > 0L, z[1L], Inf))
  ans <- rep(Inf, length(n)); hit <- n > 0L
  ans[hit] <- z[n[hit]]
  ans
}

#' Apply a tier's direction, margin or value cutoff to two summaries
#' @param x,y Numeric summaries aligned to one or more treatment/control pairs.
#' @param sp Metadata of gpc_endpoint.
#' @return Signed scores; a required missing measurement defers the look.
#' @noRd
.pair_tier_contrast <- function(x, y, sp) {
  if (anyNA(x) || anyNA(y))
    .gpc_defer(paste0("Missing deciding measurement for channel ", sp$channel, "."))
  if (!is.null(sp$threshold)) {
    x <- as.integer(x > sp$threshold)
    y <- as.integer(y > sp$threshold)
  }
  d <- x - y
  # Both patients event-free: no observed timing distinction.
  d[x == y] <- 0
  if (sp$direction == "lower") d <- -d
  (d > sp$margin) - (d < -sp$margin)
}

#' Build a provider for an explicitly ordered or user-supplied comparison
#' @param patients,events Validated native look records.
#' @param ti,ci Treatment/control positions in the patient table.
#' @param rule Clinical hierarchy with ordered endpoints or a single callback.
#' @param restricted Whether living pairs share the look-specific restriction horizon.
#' @param channels All validated named events and measurements for this cut, or
#' NULL on the native input path. The selection attribute identifies aliases.
#' @return A scorer with the same row and cutoff interface as the preset
#' provider. It indexes records by endpoint and subject once. For declarative
#' rules, only surviving ties need later summaries; for arbitrary functions,
#' both clinical histories are truncated once per visited pair and shared
#' across its tiers. It does not retain a matrix or a list of all pairs.
#' @noRd
.pair_tier_provider <- function(patients, events, ti, ci, rule, channels, restricted = FALSE) {
  fs <- if (!is.null(rule$compare)) list(compare = rule$compare) else rule$endpoints
  if (!length(fs)) return(function(r, i, j, at) numeric(length(j)))
  .pair_check_endpoint_channels(rule, channels)
  specs <- lapply(fs, .pair_endpoint_spec)
  pure <- !any(vapply(specs, is.null, logical(1L)))
  if (is.null(channels)) {
    if (is.null(events)) events <- data.frame(id = character(), time = numeric(),
                                             type = character(), value = numeric())
    if (!("value" %in% names(events))) events$value <- rep(NA_real_, nrow(events))
    channels <- list(events = list(hospitalization = events[events$type == "hospitalization",
                                                            c("id", "time"), drop = FALSE]),
      measurements = list(biomarker = events[events$type == "biomarker",
                                              c("id", "time", "value"), drop = FALSE]))
    selection <- list(event = "hospitalization", measurement = "biomarker")
  } else {
    selection <- attr(channels, "selection", exact = TRUE)
  }
  if (pure) {
    needed_e <- unique(unlist(lapply(specs, function(x) if (x$type != "value") x$channel)))
    needed_m <- unique(unlist(lapply(specs, function(x) if (x$type == "value") x$channel)))
    channels <- list(events = channels$events[needed_e], measurements = channels$measurements[needed_m])
  }

  n <- nrow(patients)
  index <- function(tab, measurement) {
    k <- match(tab$id, patients$id)
    if (anyNA(k)) stop("Unenrolled or unknown patient in named endpoint table.", call. = FALSE)
    # Use exact clinical cutoffs even if the validation allowance for a
    # floating-point representation admits a time rounded just above an end.
    ix <- which(tab$time <= patients$followup[k])
    ix <- ix[order(tab$time[ix])]
    rows <- unname(split(ix, factor(k[ix], levels = seq_len(n))))
    if (measurement) {
      # Validated plain columns already have the required types and lengths.
      # Retain data.frame conversion for attributed or named numeric inputs.
      if (is.null(attributes(tab$time)) && is.null(attributes(tab$value)))
        return(lapply(rows, function(v)
          list2DF(list(time = tab$time[v], value = tab$value[v]))))
      return(lapply(rows, function(v)
        data.frame(time = tab$time[v], value = tab$value[v])))
    }
    lapply(rows, function(v) tab$time[v])
  }
  eh <- lapply(channels$events, index, measurement = FALSE)
  mh <- lapply(channels$measurements, index, measurement = TRUE)
  if (pure) {
    sources <- lapply(specs, function(sp) if (sp$type == "value") mh[[sp$channel]] else eh[[sp$channel]])
    return(function(r, i, j, at) {
      ans <- numeric(length(j))
      for (l in seq_along(fs)) {
        todo <- which(ans == 0)
        if (!length(todo)) break
        sp <- specs[[l]]; z <- sources[[l]]
        a <- at[todo]
        x <- .pair_tier_value(z[[ti[i]]], a, sp$type)
        y <- vapply(seq_along(todo), function(s)
          .pair_tier_value(z[[ci[j[todo[s]]]]], a[s], sp$type), numeric(1L))
        ans[todo] <- .pair_tier_contrast(x, y, sp)
      }
      ans
    })
  }

  checked <- function(x) structure(x, class = "gpcseq_channels")
  clip <- function(person, at) {
    e <- checked(lapply(eh, function(z) {
      v <- z[[person]]
      v[seq_len(findInterval(at, v))]
    }))
    m <- checked(lapply(mh, function(z) {
      v <- z[[person]]
      # Plain indexed histories have automatic row names. Slice their columns
      # directly; preserve generic data-frame behavior for attributed inputs.
      ix <- seq_len(findInterval(at, v$time))
      if (is.null(attributes(v$time)) && is.null(attributes(v$value)) &&
          .row_names_info(v, 1L) <= 0L)
        list2DF(list(time = v$time[ix], value = v$value[ix]))
      else v[ix, , drop = FALSE]
    }))
    list(events = e, measurements = m)
  }
  history <- function(person, at) {
    out <- get_channels(person, at)
    if (!is.null(selection$event)) out$admissions <- out$events[[selection$event]]
    if (!is.null(selection$measurement)) out$biomarker <- out$measurements[[selection$measurement]]
    checked(out)
  }
  # A living restricted pair is ascertained only when both patients reach
  # the horizon. All callback histories for a given patient thus have one
  # cutoff within this fit. Keep invoking every callback in its original order.
  get_channels <- clip
  if (restricted) {
    clipped <- vector("list", n)
    get_channels <- function(person, at) {
      if (is.null(clipped[[person]])) clipped[[person]] <<- clip(person, at)
      # A callback may modify a measurement table by reference. Preserve a fresh
      # pair-local history, including shared aliases within that callback.
      .Call(gpc_copy_channels, clipped[[person]])
    }
  }
  function(r, i, j, at) {
    vapply(seq_along(j), function(s) {
      t <- at[s]; k <- j[s]
      left <- history(ti[i], t); right <- history(ci[k], t)
      for (l in seq_along(fs)) {
        z <- fs[[l]](left, right, t)
        if (!is.numeric(z) || length(z) != 1L || !is.null(dim(z)) || !(z %in% c(-1, 0, 1)))
          stop("Custom compare tier ", names(fs)[l],
               " must return a single number -1, 0 or 1, for patients ",
               patients$id[ti[i]], " and ", patients$id[ci[k]], " at ", t, ".", call. = FALSE)
        if (z != 0) return(as.numeric(z))
      }
      0
    }, numeric(1L))
  }
}
