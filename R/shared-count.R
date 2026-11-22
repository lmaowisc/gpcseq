#' @useDynLib gpcseq, .registration = TRUE
NULL

# The plain count-only preset compares survival then the number of admissions.
# Custom callbacks and declarative endpoint rules keep their own scorer.
.shared_count_eligible <- function(restricted, margin, rule) {
  !restricted && is.null(margin) && inherits(rule, "gpcseq_hierarchy") &&
    isTRUE(rule$admission_count) && identical(rule$admission_time, "none") &&
    is.null(rule$biomarker_margin) && is.null(rule$biomarker_threshold) &&
    is.null(rule$compare) && is.null(rule$endpoints)
}

.shared_count_summaries <- function(patients, events, ti, ci) {
  use <- events$type == "hospitalization"
  records <- .Call(gpc_shared_count, as.double(patients$followup),
                   as.logical(patients$death), as.integer(patients$arm),
                   match(events$id[use], patients$id), as.double(events$time[use]))
  pairs <- as.double(length(ti)) * length(ci)
  mean <- colSums(records[ti,,drop=FALSE]) / pairs
  records[ti,] <- records[ti,,drop=FALSE] / length(ci)
  records[ci,] <- records[ci,,drop=FALSE] / length(ti)
  list(mean=mean, average=records, n_pairs=pairs, n_ascertained=pairs)
}
