#' Prepare patient histories for a single sequential trial data cut
#'
#' @description Convert a patient table with named clinical history channels,
#'   or a `gpclearn` `gpc_outcome` object built from such tables, to the
#'   native patient, event and transition inputs and retain all named channels
#'   for a hierarchy with multiple prespecified endpoints. This conversion
#'   changes representation only unless administrative follow-up is explicitly
#'   declared below; it does not select a look, censor follow-up,
#'   choose an endpoint hierarchy, or infer an event-free interval.
#' @param subjects Data frame with one row per enrolled randomized patient,
#'   or an existing `gpc_outcome`. For the latter, use its stored IDs,
#'   follow-up, death and named channels; `treatment`, `entry` and `state0`
#'   refer to columns in its retained subject table. No dependency on
#'   gpclearn is needed to use this representation.
#' @param id,followup,death Column names in a raw subject table. IDs must be
#'   stable unique character labels or whole numbers; factor ID **labels**
#'   are accepted. `followup` is elapsed time since entry to the observed
#'   death or end of complete observation at this cut. `death` is logical or
#'   numeric 0/1 (1 = observed death at that time). Do not respecify these
#'   selectors for a `gpc_outcome`.
#' @param treatment,entry Column names for randomized assignment (1/TRUE =
#'   treatment, 0/FALSE = control) and nonnegative calendar time of entry
#'   on the design clock. Dates must already be converted to the declared
#'   units. Treatment and death labels are never decoded from factor levels.
#' @param events,measurements Named lists of clinical data frames for raw
#'   subjects. Event channels have `id,time`; each row counts as
#'   one occurrence of the prespecified recurrent event (e.g. admissions or
#'   exacerbations). Measurement channels add numeric `value`.
#'   A named channel can have zero rows. `time` is
#'   elapsed follow-up. Optional `available` is the calendar reporting
#'   time. When an outcome object is supplied, use its own channels.
#' @param event,measurement Optional convenience selections for the built-in
#'   single-event/marker rule. `NULL` selects a sole channel; with multiple
#'   channels it leaves the alias unresolved (an error only when the rule
#'   needs that alias). A name selects one; `FALSE` explicitly omits the
#'   alias. Neither option removes channels from the full named lists.
#'   Use `gpc_hierarchy(endpoints = ...)` to specify several channels directly.
#' @param state0 Optional subject-column name for each patient's known
#'   initial state (character labels). Required to supply transition rows.
#' @param transitions Optional data frame with `id,time,from,to` and optional
#'   calendar `available`. States are named explicitly; arbitrary numeric
#'   state measurements in a `gpc_outcome` are not reconstructed as paths.
#' @param administrative Optional nonmissing logical vector in subject-row order.
#'   TRUE declares that the patient is alive and completely observed through
#'   `look`, with administrative censoring at that time. Supply `look` together
#'   with this vector. Named vectors must have names identical to the canonical
#'   subject IDs in order. FALSE leaves the patient's follow-up unchanged.
#' @param look Calendar time of this already prepared cut, in the same units
#'   and origin as `entry`. Used only with `administrative`.
#' @details For explicitly declared administrative rows, follow-up is computed
#'   as `look - entry`, using the same subtraction as the maturity rule. This
#'   prevents independently scaled durations from breaking mathematical ties.
#'   The supplied duration must agree within eight machine epsilons times the
#'   largest magnitude of look, entry and duration. A nonzero correction must
#'   also be no larger than the square root of machine epsilon times the larger
#'   of the supplied and recomputed durations. Otherwise, reconcile the source
#'   records and use a clock with an origin close to the trial before conversion.
#'   A recorded death cannot be declared administratively censored. All stored
#'   history times for declared patients must be at or before the recomputed
#'   duration, and reporting times must be at or before the look. Histories are
#'   retained exactly; a conflicting record requires source-data adjudication.
#'
#'   The tables must represent information known as of the intended
#'   look. The constructor does not turn a completed trial into a historical
#'   cut. Never retain a later death by clipping its follow-up to an earlier
#'   look, or remove a delayed report and declare the resulting gap event-free.
#'   Within a snapshot, known histories extend through the reported follow-up;
#'   actual calendar-limit and design-dependent clinical checks occur at fit
#'   or update. Each channel preserves its own optional reporting times;
#'   all present dates are checked when the full object is passed as `data`.
#'   Undated channels do not inherit a reporting date. When selected nonempty
#'   channels differ in whether they are dated, the pooled convenience view
#'   keeps the known dates and uses NA for undated rows. Use `data = cut`;
#'   a native fit rejects an incomplete reporting-time column.
#' @return A `gpcseq_trial_data` list with `patients` (standard id, arm,
#'   entry, followup, death and optional state0), `events` (id, time, type,
#'   value and optional available), and `transitions` (id, time, from, to
#'   and optional available). `channels$events` and `channels$measurements`
#'   retain one canonical table per original name. `selection` records the
#'   convenience channel names, their ambiguity flags, and whether the flat
#'   reporting-time column contains undated rows because coverage differed. Input
#'   patient order and stable IDs are retained. With administrative declarations,
#'   `administrative` records the look and an audit table containing each
#'   declared ID, original and canonical durations, correction and error bound.
#' @seealso `gpc_fit()`, `gpc_update()`, `gpc_ipcw_effects()`
#' @examples
#' # A frozen cut: follow-up is complete through the time in each subject row.
#' subjects <- data.frame(pid = 1:4, assigned = c(1, 1, 0, 0), entry = 0,
#'                        end = c(2, 1.5, 1, 2), died = c(0, 1, 1, 0))
#' admissions <- data.frame(id = 4L, time = 0.5)
#' dat <- gpc_trial_data(subjects, id = "pid", followup = "end", death = "died",
#'                       treatment = "assigned", events = list(admit = admissions))
#' plan <- gpc_design(2, .05, "ipcw", horizons = 2)
#' gpc_fit(plan, 1, data = dat)$estimate
#' # Both recurrent endpoints can be used by name, with no default alias.
#' many <- gpc_trial_data(subjects, id = "pid", followup = "end", death = "died",
#'                        treatment = "assigned",
#'                        events = list(admit = admissions, flare = admissions[0, ]))
#' names(many$channels$events)
#' # Explicit administrative censoring for subjects 1 and 4 at calendar time 2.
#' cut <- gpc_trial_data(subjects, id = "pid", followup = "end", death = "died",
#'                       treatment = "assigned", look = 2,
#'                       administrative = c(TRUE, FALSE, FALSE, TRUE))
#' cut$administrative$audit
#' @export
gpc_trial_data <- function(subjects, id = "id", followup = "followup", death = "death",
                           treatment = "arm", entry = "entry", events = list(),
                           measurements = list(), event = NULL, measurement = NULL,
                           state0 = NULL, transitions = NULL,
                           administrative = NULL, look = NULL) {
  if (xor(is.null(administrative), is.null(look)))
    stop("Supply administrative and look together, or omit both.", call. = FALSE)
  outcome <- inherits(subjects, "gpc_outcome")
  channel_ids <- NULL
  if (outcome) {
    if (!missing(id) || !missing(followup) || !missing(death) ||
        !missing(events) || !missing(measurements))
      stop("A gpc_outcome supplies its own IDs, follow-up, death and history channels; ",
           "omit id, followup, death, events and measurements.", call. = FALSE)
    o <- subjects
    need <- c("ids","end","death","subjects","id_column","events","measurements","states")
    if (!is.list(o) || !all(need %in% names(o)) || anyDuplicated(names(o)))
      stop("Unsupported gpc_outcome: the validated history fields are missing or duplicated.", call. = FALSE)
    subjects <- o$subjects; .gpc_frame(subjects, character(), "gpc_outcome$subjects")
    k <- .gtd_field(subjects, o$id_column, "gpc_outcome id_column")
    stored <- .gtd_ids(o$ids, "gpc_outcome$ids", unique = TRUE)
    key <- .gtd_ids(k, "gpc_outcome subject ID column", unique = TRUE)
    # gpclearn uses as.character(original_id), which can record large numeric
    # IDs in scientific notation depending on the display settings at creation.
    # Verify numeric identity, retain original side-table keys for matching,
    # and expose the same decimal IDs as the raw-data constructor.
    if (is.numeric(k) && !is.object(k)) {
      z <- suppressWarnings(as.numeric(stored))
      aligned <- length(z) == length(k) && all(is.finite(z)) &&
        identical(unname(as.numeric(k)), unname(z))
    } else aligned <- identical(key, stored)
    if (!aligned)
      stop("gpc_outcome$ids do not match the retained subject IDs in order.", call. = FALSE)
    channel_ids <- stored
    fup <- o$end; dead <- o$death; events <- o$events; measurements <- o$measurements
    if (length(fup) != length(key) || length(dead) != length(key))
      stop("gpc_outcome has inconsistent subject-level vector lengths.", call. = FALSE)
  } else {
    .gpc_frame(subjects, character(), "subjects")
    key <- .gtd_ids(.gtd_field(subjects, id, "id"), "subjects id", unique = TRUE)
    fup <- .gtd_field(subjects, followup, "followup")
    dead <- .gtd_field(subjects, death, "death")
  }
  # Scientific endpoint/arm codes must be explicit. Do not convert factor
  # levels by their integer positions or use final-trial status to replace NA.
  arm <- .gtd_binary(.gtd_field(subjects, treatment, "treatment"), "treatment")
  dead <- .gtd_binary(dead, "death")
  en <- .gtd_field(subjects, entry, "entry")
  .gtd_time(en, "entry"); .gtd_time(fup, "followup")
  p <- data.frame(id = key, arm = as.integer(arm), entry = as.numeric(en),
                  followup = as.numeric(fup), death = as.logical(dead),
                  stringsAsFactors = FALSE)
  if (!is.null(state0)) {
    st <- .gtd_field(subjects, state0, "state0")
    if (!is.character(st) || anyNA(st) || any(!nzchar(st)))
      stop("Initial states must be explicit, nonmissing character labels.", call. = FALSE)
    p$state0 <- st
  } else if (!is.null(transitions) && is.data.frame(transitions) && nrow(transitions)) {
    stop("Supply state0 (the name of the subject entry-state column) with transitions.", call. = FALSE)
  }
  audit <- NULL
  if (!is.null(administrative)) {
    adjusted <- .gtd_administrative(p, look, administrative)
    p <- adjusted$patients
    audit <- adjusted$audit
  }
  # The two legacy convenience channels are views of the full named data,
  # not a filter on what the arbitrary clinical rule may subsequently use.
  channels <- list(events = .gtd_channels(events, "events", p, channel_ids),
                   measurements = .gtd_channels(measurements, "measurements", p, channel_ids))
  a <- .gtd_select(channels$events, event, "events")
  b <- .gtd_select(channels$measurements, measurement, "measurements")
  eh <- .gtd_view(channels$events, a$name, "hospitalization")
  bh <- .gtd_view(channels$measurements, b$name, "biomarker")
  ed <- "available" %in% names(eh); bd <- "available" %in% names(bh)
  mixed <- nrow(eh) > 0L && nrow(bh) > 0L && xor(ed, bd)
  if (mixed) {
    # Preserve the dates even if this table is manually unpacked. NA means
    # the original channel supplied no timestamp; it is not an assumed date.
    # The prepared data= path validates each original table before fitting;
    # unpacking this partial column into a native fit must instead fail.
    if (!ed) eh$available <- rep(NA_real_, nrow(eh))
    if (!bd) bh$available <- rep(NA_real_, nrow(bh))
  } else {
    keep <- (ed && nrow(eh)>0L) || (bd && nrow(bh)>0L) ||
            (!nrow(eh) && !nrow(bh) && (ed || bd))
    if (keep) {
      if (!ed) eh$available <- numeric()
      if (!bd) bh$available <- numeric()
    } else {
      eh$available <- NULL; bh$available <- NULL
    }
  }
  e <- rbind(eh, bh); rownames(e) <- NULL
  tr <- .gtd_long(transitions, "transition", p, channel_ids)
  selection <- list(event = a$name, measurement = b$name, event_ambiguous = a$ambiguous,
                    measurement_ambiguous = b$ambiguous, mixed_available = mixed)
  out <- structure(list(patients = p, events = e, transitions = tr,
                        channels = channels, selection = selection), class = "gpcseq_trial_data")
  if (!is.null(administrative)) {
    .gtd_administrative_histories(out, look, administrative)
    out$administrative <- list(look = look, audit = audit)
  }
  out
}

#' Select an explicitly named column from a validated data frame
#' @param x Subject data frame; names have already been checked for uniqueness.
#' @param sel One nonmissing, nonempty column name. No partial matching.
#' @param role Field description used in an error identifying the missing name.
#' @return The selected unconverted vector, with its original type.
#' @noRd
.gtd_field <- function(x, sel, role) {
  if (!is.character(sel) || length(sel) != 1L || is.na(sel) ||
      !nzchar(sel) || !sel %in% names(x))
    stop(role, " must name an existing column in the subject data.", call. = FALSE)
  x[[sel]]
}

#' Normalize identifiers without deriving identity from row positions
#' @param x Character labels, factor ID labels, or ordinary finite whole numbers
#'   exactly representable as doubles. Blanks, text NA sentinels, edge whitespace
#'   and control characters are rejected; numeric 001 and character "001" are
#'   not treated as interchangeable identifiers.
#' @param name Diagnostic field name.
#' @param unique Require unique subject keys; long-table records may legitimately repeat IDs.
#' @return Unnamed character IDs in original order. Character IDs are not
#'   trimmed or recoded. Numeric IDs are formatted in full decimal notation.
#' @noRd
.gtd_ids <- function(x, name, unique = FALSE) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) && !is.object(x) && is.null(dim(x))) {
    if (any(!is.finite(x)) || any(x != trunc(x)) || any(abs(x) > 2^53 - 1))
      stop(name, " numeric IDs must be finite exact whole numbers.", call. = FALSE)
    x[x == 0] <- 0
    x <- sprintf("%.0f", x)
  }
  if (!is.character(x) || !is.null(dim(x)) || anyNA(x) ||
      any(!nzchar(x)) || any(trimws(x) != x) || any(x %in% c("NA","<NA>")) ||
      any(grepl("[[:cntrl:]]", x)))
    stop(name, " needs nonmissing IDs without blank/sentinel labels, edge whitespace or control characters.",
         call. = FALSE)
  x <- unname(x)
  if (unique && anyDuplicated(x))
    stop(name, " has a duplicate normalized ID: ", x[which(duplicated(x))[1L]], ".", call. = FALSE)
  x
}

#' Require explicit logical or numeric binary coding of an endpoint or arm
#' @param x Input subject vector; classed vectors such as factors are rejected.
#' @param role Diagnostic field name.
#' @return Unnamed 0/1 integer vector, aligned to subject order.
#' @noRd
.gtd_binary <- function(x, role) {
  if (is.object(x) || !is.null(dim(x)) || !(is.logical(x) || is.numeric(x)) ||
      anyNA(x) || any(!x %in% c(0,1)))
    stop(role, " must be logical or explicitly coded numeric 0/1; recode labels before use.", call. = FALSE)
  as.integer(x)
}

#' Require finite, nonnegative time in one already established unit
#' @param x Raw elapsed or calendar vector; dates, factors, and matrices must
#'   be mapped externally, using the same time origin and unit for the design.
#' @param name Field used in any diagnostic.
#' @return Invisibly NULL. This validation neither truncates nor imputes times.
#' @noRd
.gtd_time <- function(x, name) {
  if (!is.numeric(x) || is.object(x) || !is.null(dim(x)) || any(!is.finite(x)) || any(x < 0))
    stop(name, " must contain finite nonnegative numeric times in the design units.", call. = FALSE)
  invisible(NULL)
}

#' Validate the names and storage of a collection of clinical channels
#' @param channels An empty or uniquely named list of data frames.
#' @param what Collection type (events or measurements) for diagnostics.
#' @return Invisibly NULL; malformed, duplicated or blank names are rejected.
#' @noRd
.gtd_channel_names <- function(channels, what) {
  if (!is.list(channels) || is.data.frame(channels))
    stop(what, " must be a named list of data frames.", call. = FALSE)
  nm <- names(channels)
  if (length(channels) && (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) ||
                          anyDuplicated(nm)))
    stop(what, " must have distinct, nonempty channel names.", call. = FALSE)
  invisible(NULL)
}

#' Normalize every supplied history channel without combining endpoint types
#' @param channels Named source list; a configured zero-row data frame is retained.
#' @param what One of events or measurements.
#' @param p Canonical enrolled patient table.
#' @param channel_ids Optional source ID aliases supplied by gpclearn.
#' @return Named list in input order with id,time for event channels or
#'   id,time,value for measurement channels; available is retained where supplied.
#'   Patient records are kept in source row order without pooling channels.
#' @noRd
.gtd_channels <- function(channels, what, p, channel_ids) {
  .gtd_channel_names(channels, what)
  biomarker <- identical(what,"measurements")
  out <- lapply(names(channels), function(nm) {
    tab <- channels[[nm]]
    if (!is.data.frame(tab))
      stop(what," channel '",nm,"' must be a data frame, including when empty.",call.=FALSE)
    what_type <- if (biomarker) "biomarker" else "hospitalization"
    z <- tryCatch(.gtd_long(tab,what_type,p,channel_ids),
                  error=function(err) stop(what," channel '",nm,"': ",conditionMessage(err),call.=FALSE))
    if (biomarker && any(.gpc_repeated_measurement(z$id, z$time)))
      stop("measurement channel '", nm, "' has duplicate values for an ID at the same time.",
           call.=FALSE)
    z
  })
  names(out) <- names(channels)
  out
}

#' Resolve one convenience alias without discarding the original source channels
#' @param channels Normalized named tables of one clinical type.
#' @param selector Single channel name, NULL for the sole channel, FALSE to omit.
#' @param what Label used for invalid selection diagnostics.
#' @return A list of name (NULL if omitted or ambiguous) and ambiguous (TRUE
#'   only for multiple supplied channels with NULL selector).
#' @noRd
.gtd_select <- function(channels, selector, what) {
  if (identical(selector,FALSE)) return(list(name=NULL, ambiguous=FALSE))
  if (is.null(selector)) {
    if (length(channels)==1L) return(list(name=names(channels)[1L], ambiguous=FALSE))
    return(list(name=NULL, ambiguous=length(channels)>1L))
  }
  if (!is.character(selector) || length(selector)!=1L || is.na(selector) ||
      !selector %in% names(channels)) stop("Unknown ",what," channel.",call.=FALSE)
  list(name=unname(selector), ambiguous=FALSE)
}

#' Form the legacy one-type view over the selected clinical history
#' @param channels Canonical named clinical tables.
#' @param key Selected channel name or NULL when no convenience channel is known.
#' @param type Historical native tag (hospitalization or biomarker).
#' @return A data frame id,time,type,value and optional available; an ambiguous
#'   or explicitly absent selection yields zero rows, never a merged channel.
#' @noRd
.gtd_view <- function(channels,key,type) {
  if (is.null(key))
    return(data.frame(id=character(),time=numeric(),type=character(),value=numeric()))
  z <- channels[[key]]
  z$type <- rep(type,nrow(z))
  if (type=="hospitalization") z$value <- rep(NA_real_,nrow(z))
  z[,c("id","time","type","value",intersect("available",names(z))),drop=FALSE]
}

#' Check every stored clinical channel against the actual calendar data cut
#' @param channels The full canonical list of event and measurement channel
#'   lists. Null is allowed only for compatibility with old/native inputs.
#' @param p Validated look-specific patient table with stable id, entry and
#'   followup. An unseen patient may not appear in any channel.
#' @param look Scheduled numeric calendar time on the same clock as entry and
#'   available, not its ordinal index; Inf skips only the calendar limit.
#' @return Invisibly NULL. Checks all channels, including those not selected as
#'   convenience aliases. Dated and undated channels remain distinct. A
#'   measurement channel may have at most one record per patient and time.
#'   This checks chronology; it never infers that an undated report was timely.
#' @noRd
.gtd_validate_look <- function(channels,p,look) {
  if (is.null(channels)) return(invisible(NULL))
  if (!is.numeric(look) || length(look)!=1L || is.na(look) || look<0)
    stop("Invalid numeric calendar time in channel validation.",call.=FALSE)
  if (!is.list(channels) || is.null(names(channels)) || anyDuplicated(names(channels)) ||
      !all(c("events","measurements") %in% names(channels)))
    stop("Expected the original named event and measurement channel lists.",call.=FALSE)
  for (kind in c("events","measurements")) {
    tabs <- channels[[kind]]; .gtd_channel_names(tabs,kind)
    for (nm in names(tabs)) {
      z <- tabs[[nm]]; what <- paste0(kind," channel '",nm,"'")
      measure <- kind=="measurements"
      .gpc_frame(z,c("id","time",if (measure) "value"),what)
      if (!is.character(z$id)) stop(what," must retain canonical character IDs.",call.=FALSE)
      idx <- match(z$id,p$id)
      .gpc_bad(is.na(idx),z$id,paste(what,"has an unenrolled ID"))
      .gtd_time(z$time,paste(what,"time"))
      occurrence <- p$entry[idx] + z$time
      .gpc_bad(z$time>p$followup[idx]+1e-12 | occurrence>look+1e-12,
               z$id,paste(what,"occurrence outside observed follow-up at this look"))
      if ("available" %in% names(z)) {
        .gtd_time(z$available,paste(what,"available"))
        .gpc_bad(z$available < occurrence-1e-12 | z$available>look,
                 z$id,paste(what,"availability must follow occurrence and precede or equal the look"))
      }
      if (measure) {
        if (!is.numeric(z$value) || is.object(z$value) || !is.null(dim(z$value)) ||
            any(is.infinite(z$value))) stop(what," values must be numeric, finite or NA.",call.=FALSE)
        dup <- .gpc_repeated_measurement(z$id, z$time)
        .gpc_bad(dup,z$id,paste(what,"has more than one measurement at the same time"))
      }
    }
  }
  invisible(NULL)
}

#' Standardize a selected clinical event, measurement or transition channel
#' @param tab Selected long table (or NULL). id is the standard side-table
#'   key even when the main patient table has a user-specified ID column name.
#' @param type One of hospitalization, biomarker or transition.
#' @param p Standardized patient table used to validate IDs and observed times.
#' @param channel_ids Optional original stored gpclearn labels, in patient
#'   order. These allow an ID recorded in scientific notation by gpclearn
#'   to match its exact decimal numeric subject ID. Text IDs stay literal.
#' @return Standard rows in source order; optional calendar available column
#'   preserved. Clinical event/measurement times cannot exceed observed end.
#'   Calendar look and design state checks remain the responsibility of fitting.
#' @noRd
.gtd_long <- function(tab, type, p, channel_ids = NULL) {
  transition <- type == "transition"; biomarker <- type == "biomarker"
  need <- c("id","time", if (transition) c("from","to") else if (biomarker) "value")
  if (is.null(tab)) {
    tab <- data.frame(id=character(),time=numeric(),stringsAsFactors=FALSE)
    if (transition) {tab$from <- character(); tab$to <- character()}
    if (biomarker) tab$value <- numeric()
  }
  .gpc_frame(tab, need, type)
  key <- .gtd_ids(tab$id, paste(type, "id")); j <- match(key,p$id)
  if (!is.null(channel_ids)) {
    old <- match(key, channel_ids)
    if (any(!is.na(old) & !is.na(j) & old != j))
      stop(type, " contains an ambiguous stored ID.", call. = FALSE)
    missing <- is.na(j)
    j[missing] <- old[missing]
  }
  if (anyNA(j)) stop(type," contains an unknown patient ID: ",key[which(is.na(j))[1L]],".",call.=FALSE)
  key <- p$id[j]
  .gtd_time(tab$time, paste(type,"time"))
  .gpc_bad(tab$time > p$followup[j] + 1e-12, key, paste(type,"outside observed follow-up"))
  z <- data.frame(id=key,time=as.numeric(tab$time),stringsAsFactors=FALSE)
  if (transition) {
    for (nm in c("from","to")) {
      y <- tab[[nm]]
      if (!is.character(y) || anyNA(y) || any(!nzchar(y)))
        stop("transition ",nm," must contain nonmissing state labels.",call.=FALSE)
      z[[nm]] <- y
    }
  } else if (biomarker) {
    val <- tab$value
    if (is.logical(val) && all(is.na(val))) val <- as.numeric(val)
    if (!is.numeric(val) || is.object(val) || !is.null(dim(val)) || any(is.infinite(val)))
      stop("Biomarker values must be numeric, finite or NA.",call.=FALSE)
    z$value <- as.numeric(val)
  }
  if ("available" %in% names(tab)) {
    .gtd_time(tab$available, paste(type,"available"))
    .gpc_bad(tab$available < p$entry[j] + z$time - 1e-12,key,
             paste(type,"reporting precedes occurrence"))
    z$available <- as.numeric(tab$available)
  }
  z
}
