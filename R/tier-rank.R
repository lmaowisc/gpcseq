# Describe the ordered comparisons after survival in a restricted analysis.
# NULL selects the existing pair evaluator. An empty list gives survival alone.
.rank_tiers <- function(rule) {
  if (is.null(rule)) return(NULL) # Older saved plans retain reference behavior.
  if (!is.null(rule$compare) || !is.null(rule$endpoints)) {
    tiers <- if (is.null(rule$compare)) rule$endpoints else list(rule$compare)
    specs <- lapply(tiers, .pair_endpoint_spec)
    if (any(vapply(specs, is.null, logical(1)))) return(NULL)
    # A nonterminal margin can produce overlapping ties, so it cannot refine
    # the subjects into disjoint blocks. The last margin needs only two counts.
    if (length(specs) > 1L && any(vapply(utils::head(specs, -1L),
        function(s) s$margin > 0, logical(1)))) return(NULL)
    return(specs)
  }
  tier <- function(channel, type, direction, threshold = NULL, margin = 0) {
    list(channel = channel, type = type, direction = direction,
         threshold = threshold, margin = margin)
  }
  specs <- list()
  if (rule$admission_count) specs <- list(tier("hospitalization", "count", "lower"))
  if (rule$admission_time != "none")
    specs <- c(specs, list(tier("hospitalization", rule$admission_time, "higher")))
  if (!is.null(rule$biomarker_margin) || !is.null(rule$biomarker_threshold)) {
    margin <- if (is.null(rule$biomarker_margin)) 0 else rule$biomarker_margin
    specs <- c(specs, list(tier("biomarker", "value", rule$biomarker_direction,
                                rule$biomarker_threshold, margin)))
  }
  specs
}

# Use the validated named channels for declarative tiers. Preset tiers always
# use the selected native view, even when the input contains additional names.
.rank_channel_key <- function(spec, ids, tau, events, channels) {
  measured <- spec$type == "value"
  if (is.null(channels)) {
    tab <- events
    if (!is.null(tab)) tab <- tab[tab$type == spec$channel, , drop = FALSE]
  } else {
    nameset <- if (measured) channels$measurements else channels$events
    tab <- nameset[[spec$channel]]
  }
  known <- tab$time <= tau
  k <- match(tab$id[known], ids, nomatch = 0L)
  use <- k != 0L
  k <- k[use]; t <- tab$time[known][use]
  if (spec$type == "count") {
    x <- tabulate(k, length(ids))
  } else {
    x <- rep(if (measured) NA_real_ else Inf, length(ids))
    if (length(k)) {
      ord <- order(k, t, method = "radix")
      k <- k[ord]
      use_row <- !duplicated(k, fromLast = spec$type != "first")
      vals <- if (measured) tab$value[known][use] else t
      x[k[use_row]] <- vals[ord[use_row]]
    }
    if (!is.null(spec$threshold)) x <- as.numeric(x > spec$threshold)
  }
  if (spec$direction == "lower") -as.double(x) else as.double(x)
}

# For every focal patient, count opposite-arm scores lower/higher by strictly
# more than the terminal margin. Earlier prefix groups sort by absolute
# rank; binary searches operate only within the tied prefix. Compare the
# original subtraction with the margin as in .pair_tier_contrast. Adding a
# margin to a patient's score before searching can change a floating-point tie.
.rank_margin_counts <- function(group, key, focal, other, margin) {
  ord <- order(group[other], key[other], method = "radix")
  other_key <- key[other[ord]]
  offsets <- c(0L, cumsum(tabulate(group[other], max(group))))
  g <- group[focal]
  bound <- function(upper) {
    lo <- offsets[g]; hi <- offsets[g + 1L]
    repeat {
      todo <- which(lo < hi)
      if (!length(todo)) break
      mid <- lo[todo] + ceiling((hi[todo] - lo[todo]) / 2)
      left <- key[focal[todo]]; right <- other_key[mid]
      delta <- left - right
      delta[left == right] <- 0 # both event-free infinities are tied.
      before <- if (upper) delta >= -margin else delta > margin
      lo[todo[before]] <- mid[before]
      hi[todo[!before]] <- mid[!before] - 1
    }
    lo
  }
  list(better = bound(FALSE), worse = length(other) - bound(TRUE))
}

# Return one ordered group for each survivor, or use a positive terminal
# margin to return the opposite-arm counts directly. Compute a tier only for
# groups still containing both arms; unmixed groups cannot raise a missing
# deciding measurement. Costs depend on records, tiers and subjects rather
# than on the cross product of the two arms.
.rank_tier_counts <- function(specs, ids, nt, tau, events, channels) {
  n <- length(ids)
  tr <- seq_len(nt); co <- seq.int(nt + 1L, n)
  group <- rep.int(1L, n)
  for (sp in specs) {
    mixed <- tabulate(group[tr], max(group)) > 0L & tabulate(group[co], max(group)) > 0L
    active <- which(mixed[group])
    if (!length(active)) break
    key <- numeric(n)
    key[active] <- .rank_channel_key(sp, ids[active], tau, events, channels)
    if (anyNA(key)) return(NULL)
    if (sp$margin > 0) {
      t <- .rank_margin_counts(group, key, tr, co, sp$margin)
      c <- .rank_margin_counts(group, key, co, tr, sp$margin)
      return(list(win1 = t$better, loss1 = t$worse, loss0 = c$better, win0 = c$worse))
    }
    ord <- order(group, key, method = "radix")
    # Equality is safe for adjacent event-free histories (+Inf == +Inf).
    next_group <- c(TRUE, group[ord[-1L]] != group[ord[-n]] |
                            key[ord[-1L]] != key[ord[-n]])
    group[ord] <- cumsum(next_group)
  }
  p1 <- c(0, cumsum(as.double(tabulate(group[tr], max(group)))))
  p0 <- c(0, cumsum(as.double(tabulate(group[co], max(group)))))
  list(win1 = p0[group[tr]], loss1 = length(co) - p0[group[tr] + 1L],
       loss0 = p1[group[co]], win0 = length(tr) - p1[group[co] + 1L])
}
