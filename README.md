# gpcseq

Monitor a randomized trial in which a treatment–control comparison can change as follow-up accrues. At each scheduled analysis, `gpcseq` estimates the prespecified effect, accounts for its correlation with earlier estimates, and reports whether the trial should stop.

## Installation

Install the package from GitHub with R 4.1.0 or later. Installation from source
requires a C/C++ toolchain. To build the tutorials, also install Pandoc or use
the copy bundled with RStudio.

```r
install.packages("remotes")
remotes::install_github("lmaowisc/gpcseq", dependencies = TRUE,
                        build_vignettes = TRUE)
```

Start with `vignette("gpcseq", package = "gpcseq")`. For a worked trial with a
deferred interim look and a two-sided final analysis, open
`vignette("monitoring", package = "gpcseq")`.

## Plan the analyses

The analyst supplies the calendar times in `looks`, including the final analysis. For example, `looks = c(12, 24, 36)` specifies two interim looks and a final analysis at month 36. By default the trial uses O’Brien–Fleming-like spending at a two-sided level of 0.05. The other named choices are Pocock-like (`"pocock"`), Kim–DeMets (`"power"`), and Hwang–Shih–DeCani (`"hsd"`). Alternatively, `spending = c(.005, .02, .05)` specifies the cumulative alpha at those three looks directly.

## Monitor a mock trial

The package ships `topcat_mock`, a simulated dataset patterned after the
TOPCAT trial in heart failure: 500 patients randomized 1:1 (`arm` 1 =
active treatment, 0 = placebo), enrolled over 66 months, with analyses
scheduled at calendar months 18, 36, 54, 72 and 84. `death_time` and
`hosp_time` are months from the patient's own entry; `NA` means the event
was not observed by the end of the study. All values are simulated.

```r
library(gpcseq)
data(topcat_mock)
head(topcat_mock)
#>      id arm entry death_time hosp_time
#> 1 P0001   1  35.5       41.7        NA
#> 2 P0002   0  22.7         NA        NA
#> 3 P0003   1  24.4       32.6      24.2
#> 4 P0004   0  24.7         NA       0.1
#> 5 P0005   1  65.2         NA       6.6
#> 6 P0006   0  41.8       21.4       6.1
```

A monitoring look at calendar month `u` sees only what has happened by
then: the patients enrolled so far, each followed from entry to death or
to month `u`, and the hospitalizations observed in that window. The
helper below reconstructs this data cut from the final database.

```r
snapshot <- function(u) {
  p <- topcat_mock[topcat_mock$entry <= u, ]
  reach <- u - p$entry
  p$followup <- pmin(p$death_time, reach, na.rm = TRUE)
  p$death <- !is.na(p$death_time) & p$death_time <= reach
  seen <- !is.na(p$hosp_time) & p$hosp_time <= p$followup
  list(patients = p[, c("id", "arm", "entry", "followup", "death")],
       events = data.frame(id = p$id[seen], time = p$hosp_time[seen],
                           type = "hospitalization"))
}

peek <- function(trial) {   # Latest decision, rounded for display
  out <- gpc_results(trial)[, c("look", "time", "n0", "n1", "estimate",
                                "se", "z", "critical", "spent", "action")]
  out[, 5:9] <- round(out[, 5:9], 3)
  tail(out, 1)
}
```

Schedule the five analyses, then feed each data cut to `gpc_update()` as
its calendar time arrives. The shared-window win–loss test compares each
treatment–control pair over its shared follow-up: survival first, then
fewer hospitalizations, then a later last hospitalization if counts tie.

```r
trial <- gpc_design(looks = c(18, 36, 54, 72, 84), method = "shared")

cut1 <- snapshot(18)
trial <- gpc_update(trial, cut1$patients, cut1$events)
peek(trial)
#>   look time n0 n1 estimate    se     z critical spent   action
#> 1    1   18 66 66    0.073 0.050 1.453    4.234     0 continue
```

At month 18, 132 patients are enrolled and the net benefit — the win
probability minus the loss probability for a random treatment–control
pair — is 0.073. The O'Brien–Fleming-like boundary at this early look is
4.234, nearly uncrossable by design: almost no error is spent while most
information is still to come. The trial continues.

```r
cut2 <- snapshot(36)
trial <- gpc_update(trial, cut2$patients, cut2$events)
peek(trial)
#>   look time  n0  n1 estimate    se     z critical spent   action
#> 2    2   36 129 124    0.092 0.047 1.959    2.996 0.003 continue

cut3 <- snapshot(54)
trial <- gpc_update(trial, cut3$patients, cut3$events)
peek(trial)
#>   look time  n0  n1 estimate    se     z critical spent   action
#> 3    3   54 198 196    0.057 0.040 1.428    2.481 0.015 continue
```

The boundary descends as planned while the evidence drifts favorable but
inconclusive — the estimate even dips at month 54, when many recent
enrollees contribute only short follow-up. Each boundary is calibrated
against the estimated correlation with all earlier test statistics, so no
multiplicity is lost to the repeated analyses.

```r
cut4 <- snapshot(72)
trial <- gpc_update(trial, cut4$patients, cut4$events)
peek(trial)
#>   look time  n0  n1 estimate    se     z critical spent       action
#> 4    4   72 250 250    0.136 0.043 3.187    2.234 0.034 reject_upper
```

At month 72, with enrollment complete and follow-up matured, z = 3.187
crosses the boundary 2.234: the trial stops for benefit a year before its
scheduled end. The design is now closed — calling `gpc_update()` with the
month-84 cut would signal an error ("no update is possible after
stopping"), so the final analysis is never performed. The full monitoring
history remains available:

```r
gpc_results(trial)[, c("look", "time", "estimate", "se", "z",
                       "critical", "spent", "action")]
gpc_covariance(trial)   # Correlation structure behind the boundaries
gpc_influence(trial)    # Each patient's contribution at each look
```

For a four-patient data cut with hand-calculable answers under all three
methods, see `?gpc_example`. The [group sequential
tutorial](vignettes/gpcseq.Rmd) follows a six-patient trial through two
completed tests and shows how the correlation between looks determines
the final boundary.

## Choose the treatment effect

- **Shared-window win–loss** (`method = "shared"`): compare patients through their shared follow-up and test whether the win and loss probabilities are equal over time.
- **Restricted win–loss** (`method = "ipcw"`): estimate net benefit, win odds or win ratio through a common follow-up horizon. The standard error accounts for estimating the censoring weights.
- **Restricted net win time** (`method = "aj"`): compare time spent in more favorable clinical states, using each arm's estimated state occupation probabilities.

For either restricted method, `horizons` specifies follow-up time since entry, distinct from the calendar times in `looks`. Alternatively, `maturity = list(p = .25, cap = 24)` selects at each look a horizon up to 24 time units that at least a quarter of enrollees in each arm have had the opportunity to reach. Prespecify one primary test for stopping.

At an early IPCW look, few patients may still be observed at the horizon. The normal stopping boundary can then be crossed too often under the null, even when both arms meet the minimum risk count. The planning example in `vignette("monitoring", package = "gpcseq")` shows how to estimate the probability of any rejection under the null with the planned accrual, censoring and endpoint distributions.

To specify your own endpoint order, see `?gpc_hierarchy`. To map patient, event and measurement tables to a data cut, see `?gpc_trial_data`. Completed tests can be inspected with `gpc_results()`, `gpc_influence()` and `gpc_covariance()`.
