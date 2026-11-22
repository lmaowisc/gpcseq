# Explicit preparation of administrative durations; inference stays unchanged.
.gtd_administrative <- function(p, look, administrative) {
  .gtd_time(look, "look")
  if (length(look) != 1L)
    stop("look must be one finite nonnegative calendar time.", call. = FALSE)
  if (!is.logical(administrative) || is.object(administrative) ||
      !is.null(dim(administrative)) || length(administrative) != nrow(p) ||
      anyNA(administrative))
    stop("administrative must be a nonmissing logical vector, one value per subject.", call. = FALSE)
  if (!is.null(names(administrative)) && !identical(names(administrative), p$id))
    stop("administrative names must match canonical subject IDs in row order.", call. = FALSE)
  .gpc_bad(p$entry > look, p$id, "entry is after the declared look")
  ix <- which(administrative)
  .gpc_bad(p$death[ix], p$id[ix], "a recorded death cannot be declared administratively censored")
  target <- look - p$entry[ix]
  original <- p$followup[ix]
  delta <- target - original
  bound <- 8 * .Machine$double.eps * pmax(abs(look), abs(p$entry[ix]),
                                         abs(original), .Machine$double.xmin)
  relative_bound <- sqrt(.Machine$double.eps) * pmax(target, original)
  .gpc_bad(abs(delta) > bound | abs(delta) > relative_bound, p$id[ix],
           "administrative followup is inconsistent with look - entry; reconcile source records and clock origin")
  audit <- data.frame(id = p$id[ix], original = original, canonical = target,
                      delta = delta, bound = pmin(bound, relative_bound))
  p$followup[ix] <- target
  list(patients = p, audit = audit)
}

.gtd_administrative_histories <- function(cut, look, administrative) {
  p <- cut$patients
  check <- function(tab, label) {
    if (!nrow(tab)) return(invisible(NULL))
    j <- match(tab$id, p$id) # Constructor already validates every ID and time.
    .gpc_bad(administrative[j] & tab$time > p$followup[j], tab$id,
             paste(label, "exceeds canonical administrative followup; adjudicate source records"))
    if ("available" %in% names(tab))
      .gpc_bad(tab$available > look, tab$id,
               paste(label, "reporting availability is after the declared look"))
  }
  for (kind in c("events", "measurements"))
    for (nm in names(cut$channels[[kind]]))
      check(cut$channels[[kind]][[nm]], paste(kind, nm))
  check(cut$transitions, "transitions")
  invisible(NULL)
}
