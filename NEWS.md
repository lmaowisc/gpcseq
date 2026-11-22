# gpcseq 0.4.0

* Added `topcat_mock`, a simulated monitoring dataset patterned after
  the TOPCAT heart-failure trial. The README now walks through its
  five-look analysis one data cut at a time, showing each boundary and
  decision, with the path pinned by a regression test.
* Boundary calibration now budgets the trivariate backend's effective
  tolerance floor (1e-14 per inclusion-exclusion corner), so a request
  tighter than the backend can honor defers the look instead of
  appearing to succeed. Default tolerances are unaffected.
* `gpc_design()`, `gpc_fit()` and `gpc_update()` report more specific
  messages for invalid `horizons`, `absorbing` and `look` arguments.
* First public release on GitHub, with `URL`, `BugReports`, ORCID and
  citation metadata.

# gpcseq 0.3.0

* Opt-in arm-specific Aalen-Johansen occupation curves
  (`keep_occupation` in `gpc_fit()`/`gpc_update()`, `gpc_occupation()`).
* Exact shared-window scorer for the count-only preset.
* Faster singleton Aalen-Johansen updates and revised handling of
  extreme sequential probability integrals.

# gpcseq 0.2.0

* First reviewed installable release: shared-window, IPCW and
  Aalen-Johansen monitoring of one prespecified pairwise effect with
  calendar-time error spending, cross-look covariance from patient
  influences, and four tutorials.
