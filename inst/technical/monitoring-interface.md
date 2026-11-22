# Data cuts, fits and the monitoring log

All times use one unit. Calendar time is measured from trial opening; follow-up time is measured from the patient's entry. The scheduled look index is an integer and is neither of these times. The package accepts already prepared snapshots: it never asks the analyst to supply unobserved potential censoring or survival times. The analyst must establish that each clinical history is known through the supplied follow-up; a missing case report does not establish that no event occurred.

## Snapshot helpers

| Internal call | Input, rule and result |
|---|---|
| `.gpc_frame(x, cols, name)` | Checks that x is a data frame with unique column names and includes the named fields. An error identifies the offending table or absent fields. Extra columns are permitted. |
| `.gpc_finite(v, field)` | Requires finite numeric values. An empty vector is valid for an empty event table. A malformed column raises an input error rather than a statistical deferral. |
| `.gpc_bad(condition, id, why)` | Takes a logical predicate per patient or event row. TRUE and NA are invalid. An error gives the rule and up to four affected IDs; otherwise returns without altering data. |
| `.gpc_data_snapshot(patients, events, transitions, data, d, k)` | Accepts one patient table or one prepared trial-data object. Rejects competing sources of patients or event histories, then delegates to the same snapshot and calendar checks for all methods; no estimation or truncation occurs in the adapter. |
| `.gpc_snapshot(p, events, transitions, d, k)` | Validates the enrolled patients and the event table needed by the active method at calendar look k. Checks ID, arm, enrollment and observed follow-up, clinical dates, optional reporting dates, and baseline state where relevant. Returns the three tables as p, e and tr, supplying empty defaults; for AJ also returns aj_events, validated transition paths encoded by patient and state. Invalid paths are rejected before support checks, and the fit reuses these encoded events. It neither fits a model nor reconstructs histories using future data. |
| `.gpc_registry(p, old)` | Joins the previous cut's IDs to the new cut. All earlier participants must remain, arm and entry must agree, follow-up cannot decrease, and a known death cannot disappear or move. Returns the reduced registry for the next update. A change in row order is immaterial. An invalid input can be corrected before repeating the same scheduled look. |

A biomarker may be recorded at entry. A later visit supersedes it. An unobserved required value is not a clinical tie, but is unnecessary when survival or admissions decide a pair. For AJ, specify state ordering and absorbing states: an absence of observed exits does not establish absorption. Use a separate design object for each method. Stopping on whichever of several methods first rejects is a different trial design.

## From fit to decision

| Internal call | Input, computation and return |
|---|---|
| `.gpc_horizon(d,k,p)` | Returns NA for a shared window, the scheduled horizon for fixed restriction, or the lower of the clinical cap and the two arm-specific opportunity order statistics. For n enrollees, the descending rank is the smallest m with m/n at least p. Compare the fraction directly at an attainable threshold. There is no interpolation or monotonicity constraint across looks. An empty arm raises a deferral for maturity selection. |
| `.gpc_check_look(object,look)` | Requires a gpcseq object and an integer index between 1 and the number of scheduled looks; returns the validated design specification. An invalid request is an input error. |
| `.gpc_check_support(p,d,horizon)` | Counts enrolled participants in arm 0 then arm 1, checks the declared enrollment floor, and requires a positive finite follow-up horizon for a restricted analysis. Returns the two counts invisibly; insufficient support raises a statistical deferral. |
| `.gpc_fit_cut(x,d,k,horizon)` | Takes a validated snapshot and horizon. Checks enrollment and horizon support, dispatches once to the selected estimator, and verifies a finite estimate and one finite influence per patient. The returned fit's subjects table preserves input order. No boundary is calculated. |
| `.gpc_logrow(k,time,n,alpha,spent,direction)` | Creates a typed scheduling record; n contains control then treatment enrollment. Unavailable numbers start as NA. The scheduled spending target and last used spending are separate fields. A zero-row slice defines the result schema before any look. |
| `.gpc_defer(reason)` | Signals a gpc_deferred condition. A standalone fit exposes it; a monitoring update logs its reason, retains prior tested looks and spending, and advances to the next scheduled look, if any. Ordinary input errors do not trigger this rule. |

Standalone fits report the test-scale estimate, natural effect and keyed patient influences. The monitoring update divides the influences by the arm- and look-specific enrollment counts to form covariance. Calibration uses the boundaries of earlier tested looks. A skipped interim is not retroactively tested. If a final look cannot be tested, the trial ends without rejection. A raw estimate can exist even if it has zero variance: a standalone fit returns it, whereas monitoring records that no test was possible. The action and reason are authoritative.

Validating a snapshot requires scanning its rows; individual estimators may also index or sort follow-up times. The saved monitor retains the registry, completed fits with patient contributions, covariance and scheduling log. It does not hold full clinical histories. See the other technical guides for method-specific formulas, computational cost and worked examples. For syntax, defaults and return types of the public routines, see the installed help pages.


For named histories, the data preparation keeps each observed event and measurement channel separate. Reporting-time and follow-up checks apply before a rule can inspect them. For a preset the monitor remembers the selected active event and biomarker channel names across named cuts: an observed-event channel may be empty at a later look, but it cannot be replaced by a different outcome. For a declarative ordered rule, its endpoint specifications fix the names for every look. When an ordered rule includes a bare callback, the selected convenience names also remain fixed because the callback may use them. Use only explicitly named channels and omit the aliases with event=FALSE and measurement=FALSE if the callback does not need those shortcuts. The raw table interface uses the names hospitalization and biomarker; keep an endpoint name consistent if you mix input forms.
