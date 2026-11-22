# gpcseq

Monitor a randomized trial in which a treatment–control comparison can change as follow-up accrues. At each scheduled analysis, `gpcseq` estimates the prespecified effect, accounts for its correlation with earlier estimates, and reports whether the trial should stop.

**Companion paper.** Lu Mao and KyungMann Kim (2026). Group sequential analysis
of pairwise treatment effects under continuing follow-up. Submitted to *Biometrics*.

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

In this four-patient example we use two looks, at times 1 and 2. All four
patients enter at calendar time 0: treated patient T1 stays alive with no
admissions, treated patient T2 dies at time 1.5, control C1 dies at time 1,
and control C2 is admitted at time 0.5 and remains alive. The shared-window
test compares each pair over its available follow-up: survival first, then
fewer admissions, then a later last admission if counts tie.

```r
library(gpcseq)

death  <- c(T1 = Inf, T2 = 1.5, C1 = 1, C2 = Inf)    # Times of death
events <- data.frame(id = "C2", time = 0.5, type = "hospitalization")
snapshot <- function(u)                   # Records known by calendar time u
  data.frame(id = names(death), arm = c(1, 1, 0, 0), entry = 0,
             followup = unname(pmin(u, death)), death = unname(death <= u))

trial <- gpc_design(looks = c(1, 2), method = "shared")
for (k in 1:2) {
  trial <- gpc_update(trial, snapshot(k), events)
  if (!tail(gpc_results(trial)$action, 1) %in% c("continue", "deferred")) break
}
gpc_results(trial)[, c("look", "estimate", "se", "z", "critical", "spent", "action")]
```

The same data cut, including the state transitions used by the
Aalen--Johansen method, ships as `gpc_example()`; its help page
(`?gpc_example`) gives hand-calculable answers for all three methods.

At look 1, treatment wins all four comparisons. With no variation in the scores, the test is deferred and spends no error. At look 2, one treatment patient has died, reversing one comparison. Three wins and one loss give a mean score of 0.5 and a standard error of 0.5. The test statistic is 1, below the final boundary of 1.96; the trial ends without rejection.

The [group sequential tutorial](vignettes/gpcseq.Rmd) follows a six-patient trial through **two completed tests**. It shows how to read the interim decision and use the correlation between looks to determine the final boundary.

## Choose the treatment effect

- **Shared-window win–loss** (`method = "shared"`): compare patients through their shared follow-up and test whether the win and loss probabilities are equal over time.
- **Restricted win–loss** (`method = "ipcw"`): estimate net benefit, win odds or win ratio through a common follow-up horizon. The standard error accounts for estimating the censoring weights.
- **Restricted net win time** (`method = "aj"`): compare time spent in more favorable clinical states, using each arm's estimated state occupation probabilities.

For either restricted method, `horizons` specifies follow-up time since entry, distinct from the calendar times in `looks`. Alternatively, `maturity = list(p = .25, cap = 24)` selects at each look a horizon up to 24 time units that at least a quarter of enrollees in each arm have had the opportunity to reach. Prespecify one primary test for stopping.

At an early IPCW look, few patients may still be observed at the horizon. The normal stopping boundary can then be crossed too often under the null, even when both arms meet the minimum risk count. The planning example in `vignette("monitoring", package = "gpcseq")` shows how to estimate the probability of any rejection under the null with the planned accrual, censoring and endpoint distributions.

To specify your own endpoint order, see `?gpc_hierarchy`. To map patient, event and measurement tables to a data cut, see `?gpc_trial_data`. Completed tests can be inspected with `gpc_results()`, `gpc_influence()` and `gpc_covariance()`.
