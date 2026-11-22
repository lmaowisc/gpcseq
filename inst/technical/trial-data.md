# Clinical history exchange at a calendar look

The method answers a randomized treatment comparison at a scheduled calendar look. It uses two time coordinates: each patient's entry and the look are on one calendar clock; the observation end, admissions, markers and transitions are measured **from that patient's entry**. An entry date must be supplied even when the patient table used for a fixed analysis did not need one. A real trial date may be converted to numeric days from trial start, but the constructor does not invent entry dates.

## Supported inputs

The native interface accepts patient, clinical-event and transition data frames. `gpc_trial_data()` builds these from a raw patient table and named clinical channels, or directly from an existing `gpclearn::gpc_outcome()`. The constructor accepts either form as its first argument and returns the same class: the three native convenience tables plus the complete named histories and convenience-selection metadata. It fits no model and does not modify the original.

| Source | Target | Interpretation |
|---|---|---|
| Subject ID or stored outcome IDs | patients$id | Permanent patient identity, normalized to character; repeated event rows join by ID. |
| Randomized assignment | patients$arm | 0/1 or logical; 1 represents treatment. Never use factor order or a trial/target population flag. |
| Entry | patients$entry | Common calendar time of enrollment. Not time since entry. |
| End and death | patients$followup, patients$death | Last complete observation at this cut; death is true only for a documented death at that time. |
| Selected event channel | events of type hospitalization | Each `id,time` row counts one occurrence of the selected recurrent outcome, for example admission or exacerbation. The native code `hospitalization` is historical; labeling a channel as exacerbations preserves its scientific meaning. |
| Selected measurement channel | events of type biomarker | `id,time,value` supplies the last recorded scalar biomarker value. Under the revised default it is used, if configured, only after survival, fewer occurrences of the selected recurrent event and a later last occurrence on an equal positive count. The plan specifies the biomarker margin. |
| Explicit starting state and dated transitions | patients$state0 and transitions | Labels, not numeric factor codes; supply `id,time,from,to` for each observed change. The design specifies state favorability, absorbing states, and the death state. |

For a `gpc_outcome`, normalized `ids`, `end` and `death` are authoritative; their order is checked against the original ID column of `subjects` before joining randomized assignment and entry. If the original subject key was numeric, stored scientific-notation labels are mapped to the same exact decimal IDs produced for a raw subject table, regardless of current R display settings. The object carries the named clinical tables. Do not repeat those selectors or pass competing tables. Using an outcome object does **not** import its associated preference rules or prediction model. Visit-level numeric state values are not exact transition records; in particular, hospital **admissions** alone do not identify time currently hospitalized without **discharges**.

A named channel remains available even when it has zero records; an undeclared channel is different. A list with a single channel provides the built-in convenience alias. With several channels, a NULL selector leaves that alias unresolved; a fit using the built-in alias must select one by name, whereas an ordered rule names its own endpoint channels. `event = FALSE` or `measurement = FALSE` omits only the convenience alias and preserves the original named tables; it is the analyst's responsibility to align this with the clinical rule in the design. Extra baseline predictors are retained in the original source object and do not enter gpcseq's within-trial estimates. All long tables use a canonical `id` column even when the subject ID column has a different name. Keep digit strings such as "001" as characters throughout.

## Data available at the look

The conversion is deliberately not a historical slicing algorithm. Supply a cut whose patients, death status and follow-up were actually known at that time. Future entries, events, reporting times or observation ends are rejected when fitting; data from a closed trial must not be made plausible merely by capping follow-up at the look. For example, if a participant's only documented death occurred at study time 8, clipping their follow-up to 4 and retaining `death = TRUE` fabricates an early death even though the resulting single row might pass time-range validation.

An optional long-table `available` reporting date is measured on the calendar clock and is retained unchanged. A report about month 2 recorded at calendar month 5 cannot be used at calendar month 4. Absence of a report is not proof that there was no event: ascertainment of the complete history through the chosen observation end must be defensible. Where reporting dates do not exist, a retrospective immediate-reporting assumption must be declared by the analyst; gpcseq cannot establish it. Each stored named channel preserves its original reporting-time column, if supplied. All supplied dates are checked at the actual calendar cut, even for channels not selected as convenience aliases. When selected nonempty channels differ in reporting metadata, the flattened table retains known dates and marks rows from the undated channel with `NA_real_`; it never fills them with an assumed date. Pass the whole cut as `data = cut`. Unpacking this mixed flattened table into the native interface will reject its missing reporting times rather than silently using a delayed record.

Compatibility concerns data representation, not statistical targets. Only use a hospitalization-specific transport export if the selected recurrent event actually denotes hospital admission; the internal slot name is not evidence of that fact. An external fixed-analysis object with only endpoint totals or durations cannot support honest replay at calendar looks without actual entry and ascertainment history. A transport dataset's target-only subjects and an ML model's covariates, risk scores and trained comparison rules are not automatically inputs to the sequential estimators.

## Using an existing outcome object

If `outcome` is a gpclearn outcome **already prepared from this look's cut**, and `arm`, `enrolled` are columns of its subject table:

```r
cut <- gpc_trial_data(outcome, treatment = "arm", entry = "enrolled",
                      event = "hospitalisation", measurement = "biomarker")
fit <- gpc_fit(plan, look = 1, data = cut)
```

For Aalen–Johansen, provide an explicit baseline-state column and the separate observed transition log. A data format cannot compensate for an unidentified state-occupancy process.

## Using several endpoints without pooling them

`channels$events` is the full named list of standardized `id,time[,available]` tables; `channels$measurements` is the full named list of `id,time,value[,available]` tables. Measurement ties within the same patient, time and channel must be adjudicated; blood pressure and symptom score measured at the same time are separate observations. A typed empty channel identifies a known registered outcome with no recorded events in this cut. The constructor cannot certify that an omitted report is an observed zero.

For example, put admissions and exacerbations into separate named event tables and a symptom score into its own measurement channel. A prespecified rule can compare survival first, then admissions, then exacerbations, then symptoms by naming the relevant channel in each ordered tier. No channel need be selected as the one default event or biomarker for this rule. Patient histories supplied to a pair comparison are also truncated at its comparison time; availability by the calendar look alone does not entitle the score to use a later follow-up measurement.
