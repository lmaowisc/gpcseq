# Aalen–Johansen net win time: patient influences

This is the implementation of Example 2 and Supporting information C.3. One call to the estimator returns the restricted net win time through the requested horizon and one influence for every enrolled patient, in the input order and in the same time units. The monitoring calculation matches these values by patient and arm at different calendar looks.

## What is observed

Supply an entry state and all observed transitions, in time since entry, for each patient. Specify the state ordering from best to worst, the absorbing states, and the label for death if applicable. An absorbing state never has departures; it is not identified from a lack of observed transitions. If death is recorded, the transition into its absorbing state must agree with that patient's terminal time and death indicator. State at time zero belongs in the entry-state variable. Changes for different patients at one time use a common pretransition risk set. Within a patient, the source state must match the previous state and distinct changes need distinct adjudicated times. Clinical transitions through an observation endpoint are counted before that patient is removed. Records may include readmissions, discharges and other reversals.

Before a calendar look can be declared untestable, the public snapshot routine validates all submitted paths and terminal records, including events after the chosen follow-up horizon. A malformed record raises a data error; the open trial can be retried after correction. The fitted AJ engine accepts the indexed event table from this check, so it does not validate and sort the same paths a second time. These indices belong only to the same snapshot and patient order.

The grid starts at zero, ends at the horizon, and includes intervening transition times from both arms. We also split at observation exits to check whether an occupied transient state has become unsupported. At a grid time, the probabilities used for the next interval are **after** that time's transitions. A change exactly at the horizon has no duration in the integral.

## Forward estimate

Let B(y,z) be +1 if state y is better than z, -1 if worse, and 0 for a tie. Within each arm, start the probability column p at the entry-state proportions. At a time with transitions, the off-diagonal entry of A(y,z) is the number moving from y to z divided by the number observed in y just before that time; each diagonal is the negative sum of its off-diagonal entries. Update p to (I+A)'p. The fitted effect is the sum over the common grid of the interval width times p1' B p0. Medical follow-up after death is unnecessary, while positive fitted probability for a transient state without a patient under observation on an interval prevents testing at that horizon.

## Same derivative, calculated backwards

C.3 writes a forward occupation influence for each subject. Computing every subject's entire occupation curve is unnecessary for this scalar effect. Instead, for each arm compute how an increase in probability in each state at a grid time would change the final integrated effect. At a given interval, the instantaneous derivative v is B p0 for treatment, or B' p1 for control. Starting at the horizon with zero future value, walk backwards:

1. At the start of the current interval, the derivative **after the current jump** is the interval width times v plus the future value from the next grid time. Denote this vector by h.
2. A transition from y to z changes future value by h[z] - h[y]. Add (h[z] - h[y]) / r[y] to the pre-jump value for state y for each observed transition, where r[y] is its pretransition risk count. Together with h, this is the matrix multiplication (I+A)h.
3. For a patient who moved from y to z, the transition numerator contributes m p[y] (h[z] - h[y]) / r[y], where m is the enrolled arm size and p[y] is the **pre-jump** probability. The transition denominator also changes when a subject is upweighted. For everyone at risk in y, subtract the sum of those numerator contributions over all transitions out of y, divided once more by r[y]. A patient who is censored still contributes to pre-censoring transitions and risk sets.
4. At zero, add the derivative of the empirical baseline distribution: h[patient's entry state] minus its arm average. This term vanishes with a common initial state.

For speed, store the risk corrections as cumulative sums within state. Each observed state interval (start, end] then requires two array lookups; no loop over all patients at every transition is needed. Algebraically, the transition residual and empirical baseline term equal m times the derivative obtained by changing one patient’s weight from 1 to 1+epsilon. The unit tests check this identity against separate weighted refits.

For S states, L intervals, M changes and n enrollees, forward and reverse sums require O(SL+M+n) work once indexed. Constructing the grid/ordering and locating interval endpoints adds sorting and lookup costs. Memory is O(SL+M+n); no n-by-L-by-S array or pair matrix is formed. Each arm’s hazards are fitted once. The reverse pass reuses the fitted hazard counts, risk denominators and pre-jump probabilities.

## Internal function contracts

The seven helpers are private. Analysts call `gpc_fit()` or `gpc_update()` and extract patient contributions with `gpc_influence()`.

| Helper | Inputs | Returns and convention |
| --- | --- | --- |
| `.aj_events` | Look-specific patient and transition tables; ordered, absorbing and death-state labels | Validated records ordered by patient and time, with integer patient/state indices. Validates the full observed path, including death agreement and duplicate times, before any support-based deferral. |
| `.aj_group` | A grid index per record and the number of grid times | A list, in grid order, of the corresponding input row positions. Times without records have an empty element. |
| `.aj_add` | A numeric vector and positions with signed increments | Adds by position, summing repeated indices so tied transitions are pooled. Preserves the original vector order. |
| `.aj_forward` | One arm’s entry states, observed follow-up ends, indexed transitions, integration grid and absorbing-state declarations | The interval-by-state occupation matrix and the stored risk counts and pre-jump factors needed for derivatives. Loss of support in an occupied transient state defers the look. |
| `.aj_payoff` | The opposite arm’s interval-by-state occupation matrix, states best first | For each state and interval, the probability the opposite arm is worse minus the probability it is better. The caller negates this for control contributions. |
| `.aj_influence` | One forward fit, its signed payoff matrix, grid starts and widths | One arm-size-scaled influence per patient in arm input order. Sums to zero up to rounding. Uses cached hazards and risk sums, without case-by-case refits. |
| `.aj_fit` | Both arms at one look, transitions, horizon, state ordering and terminal-state definitions; optional snapshot-specific indexed events | Net win time in follow-up units on both reported scales, one influence per input patient and support diagnostics. Reuses prevalidated events when supplied; otherwise validates paths before testing support. |

## Four patients

All start at home (H); the state order is H, hospital (U), dead (D). T1 remains home through time 2; T2 dies at 1.5; C1 dies at 1; C2 enters hospital at 0.5 and remains there through 2. With the horizon fixed at 2 and observation to death or horizon for all four patients, the exact signed win times are:

| | C1 | C2 | Patient pair mean |
|---|---:|---:|---:|
| T1 | 1.0 | 1.5 | 1.25 |
| T2 | 0.5 | 0.5 | 0.50 |
| Patient pair mean | 0.75 | 1.00 | 0.875 |

The net win time is 0.875. In this complete-data case the influence is also each patient's pair mean minus 0.875, yielding (0.375, -0.375, -0.125, 0.125) for (T1,T2,C1,C2). This is a calculation example; it need not pass the sample-size rule of a planned trial. A second test independently refits both arm distributions after individual weight perturbations at three horizons, with censoring, tied times across different patients, varied entry states and reversible transitions.
