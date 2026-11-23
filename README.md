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

In this four-patient example we schedule two looks, at times 1 and 2.
The shared-window test compares each treated--control pair over its
shared follow-up: survival first, then fewer admissions, then a later
last admission if counts tie. Each look passes the patient records
known by that time, plus the events observed so far.

```r
library(gpcseq)                              # Load the package

events <- data.frame(id = "C2", time = 0.5,  # C2 is hospitalized at time 0.5
                     type = "hospitalization")

look1 <- data.frame(                         # Records known at calendar time 1
  id       = c("T1", "T2", "C1", "C2"),      # Two treated, two control patients
  arm      = c(1, 1, 0, 0),                  # 1 = treatment, 0 = control
  entry    = 0,                              # All enter at trial opening
  followup = c(1, 1, 1, 1),                  # Each observed through time 1
  death    = c(FALSE, FALSE, TRUE, FALSE))   # C1 has died, at time 1

look2 <- data.frame(                         # Records known at calendar time 2
  id       = c("T1", "T2", "C1", "C2"),      # Same four patients
  arm      = c(1, 1, 0, 0),
  entry    = 0,
  followup = c(2, 1.5, 1, 2),                # Follow-up ends at death or at 2
  death    = c(FALSE, TRUE, TRUE, FALSE))    # T2 has now died, at time 1.5

trial <- gpc_design(looks = c(1, 2),         # Analyses scheduled at times 1 and 2
                    method = "shared")       # Shared-window win-loss test
trial <- gpc_update(trial, look1, events)    # Look 1: treatment wins every pair --
                                             #   no score variation, test deferred
trial <- gpc_update(trial, look2, events)    # Look 2: T2's death reverses one pair
gpc_results(trial)[, c("look", "estimate", "se", "z",   # The monitoring log:
                       "critical", "spent", "action")]  #   one row per look
```

The same data cut, including the state transitions used by the
Aalen--Johansen method, ships as `gpc_example()`; its help page
(`?gpc_example`) gives hand-calculable answers for all three methods.

At look 1, treatment wins all four comparisons. With no variation in the scores, the test is deferred and spends no error. At look 2, one treatment patient has died, reversing one comparison. Three wins and one loss give a mean score of 0.5 and a standard error of 0.5. The test statistic is 1, below the final boundary of 1.96; the trial ends without rejection.

A larger simulated trial ships as the dataset `topcat_mock` (`?topcat_mock`): 500 patients monitored at five scheduled looks, stopping early at the fourth. The package tests reproduce its full monitoring path.

The [group sequential tutorial](vignettes/gpcseq.Rmd) follows a six-patient trial through **two completed tests**. It shows how to read the interim decision and use the correlation between looks to determine the final boundary.

## Choose the treatment effect

- **Shared-window win–loss** (`method = "shared"`): compare patients through their shared follow-up and test whether the win and loss probabilities are equal over time.
- **Restricted win–loss** (`method = "ipcw"`): estimate net benefit, win odds or win ratio through a common follow-up horizon. The standard error accounts for estimating the censoring weights.
- **Restricted net win time** (`method = "aj"`): compare time spent in more favorable clinical states, using each arm's estimated state occupation probabilities.

For either restricted method, `horizons` specifies follow-up time since entry, distinct from the calendar times in `looks`. Alternatively, `maturity = list(p = .25, cap = 24)` selects at each look a horizon up to 24 time units that at least a quarter of enrollees in each arm have had the opportunity to reach. Prespecify one primary test for stopping.

At an early IPCW look, few patients may still be observed at the horizon. The normal stopping boundary can then be crossed too often under the null, even when both arms meet the minimum risk count. The planning example in `vignette("monitoring", package = "gpcseq")` shows how to estimate the probability of any rejection under the null with the planned accrual, censoring and endpoint distributions.

To specify your own endpoint order, see `?gpc_hierarchy`. To map patient, event and measurement tables to a data cut, see `?gpc_trial_data`. Completed tests can be inspected with `gpc_results()`, `gpc_influence()` and `gpc_covariance()`.
