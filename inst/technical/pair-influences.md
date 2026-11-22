# Reproducing shared-window and IPCW patient influences

At one look, the input is one row per enrolled patient, with a permanent ID, randomized arm (1 treatment, 0 control), observed follow-up since entry, and a death indicator. The clinical event table supplies hospitalizations and, if required, dated biomarker measurements. An event can be used only if it was known by this look and occurred by the patient's observed cutoff. All calculations below use treatment-minus-control signs, including for control patients. The returned influence is on the same scale as the tested estimate; it has **not** been divided by arm size or by the standard error.

## Shared-window win–loss

For treatment patient i and control patient j write X_i and X_j for their observed end times. Compare their clinical histories at min(X_i,X_j), yielding J_ij = 1, −1, or 0. A first death by the comparison time settles the score; simultaneous deaths tie. If both are alive, the default favors fewer admissions, then the later last admission. The protocol can substitute the first admission or omit timing. If specified, a biomarker is compared only after remaining ties, using the latest values by that time and the prespecified clinical cutoff or pairwise difference margin. A missing deciding biomarker defers the analysis; absence of a biomarker is irrelevant if survival or hospitalization has settled the pair.

Let U be the average of J over **all** treatment–control pairs. The treatment influence for i is the row mean J_i. minus U; the control influence for j is the column mean J_.j minus U. There is no censoring-weight correction. This statistic tests equality of the win and loss curves; its mean generally depends on the censoring pattern and is not a net benefit at a fixed horizon.

## Restricted win–loss: pair scores and censoring weights

Take one horizon tau for the look. A pair is ascertained if both patients were observed through tau, or an observed first death by tau settled the comparison while its partner remained under observation. For such a pair T_ij = min(X_i,X_j,tau) is observed. This is the **first death or horizon**, not the later death. Set

    Q_W,ij = I(treatment wins)  / [G_1(T_ij) G_0(T_ij)]
    Q_L,ij = I(treatment loses) / [G_1(T_ij) G_0(T_ij)].

For a pair not ascertained, assign both components zero without evaluating an unavailable history. The sample averages W and L divide by n_1 n_0. For an observed clinical tie both components are also zero, but the pair counts as ascertained.

Fit G_a separately in the enrolled samples at this look. At each censoring time s before tau, let dN_i(s) indicate censoring at s and Y_i(s) indicate membership in the censoring risk set. With a documented death-before-censoring convention,

    Y_i(s) = I(X_i >= s) - I(X_i = s and death observed),
    Ybar(s) = sum_i Y_i(s)/n_a,
    h_s = sum_i dN_i(s) / sum_i Y_i(s),
    G_a(t) = product of (1 - h_s) over s < t.

The strict upper limit gives the estimate of P(C >= t): a patient censored exactly at t can supply a comparison at t. Before inverting the curve, require the specified minimum observed risk count and positive censoring survival at tau in both arms. If a death and censoring have indistinguishable dates and their order is unknown, the data must first be adjudicated; the stated risk convention is meaningful only for records known through the recorded time.

## Reconstructing the estimated-weight term

Choose a contrast vector v for the desired test scale:

| Effect | Estimate on test scale | v = (win coefficient, loss coefficient) |
|---|---|---|
| Net benefit | b = W − L | (1, −1) |
| Win odds | log(1+b) − log(1−b) | [2/(1−b²)] (1, −1) |
| Win ratio | log(W) − log(L) | (1/W, −1/L) |

Require −1 < b < 1 for win odds, or W,L > 0 for win ratio. The natural-scale result for either log contrast is its exponential. No clipping or pseudocount is applied. Write Q_v,ij = v_W Q_W,ij + v_L Q_L,ij. Holding the fitted weights fixed, the subject's influence is the centered row or column mean of Q_v. For a log contrast **center at the average of Q_v**, not at the reported log estimate.

For a patient i in arm a, the exact fitted relative reverse-KM influence at t is

    gamma_a(t;i)/G_a(t) =
        -sum_{s<t} [dN_i(s) - Y_i(s) h_s] / [Ybar(s) (1 - h_s)].

This follows by differentiating the empirical hazard ratio, then differentiating each product factor. The factor 1−h_s is retained in finite samples and when several censorings share a time. The influence due to estimated weights is the all-pair mean of **minus** Q_v times this relative influence, evaluated at the pair's settlement time. Thus it is equivalently

    H_v(s) = sum_{pairs with T_ij > s} Q_v,ij / (n_1 n_0)
    weight contribution of i =
       sum_s H_v(s) [dN_i(s) - Y_i(s) h_s] / [Ybar(s) (1 - h_s)].

Only ascertained pairs with a nonzero score contribute to H_v. Add this term to the fixed-weight row or column influence. The sum of the resulting influences in each arm is zero. The win-ratio contrast uses the win and loss components from the same pairs and both fitted curves; applying its coefficients before calculating the influence retains their covariance. No separate censoring fit is needed for win odds or ratio.

For a direct small-sample check, raise and lower one enrolled patient's empirical weight from 1 to 1 ± epsilon within their arm. Refit **both the arm-normalized pair means and the censoring curve** and transform the result to the desired test scale. The slope multiplied by that arm's enrollment converges to the stored subject influence. A derivative with weights held fixed checks only the first term.

## Calculation and reuse

| Internal calculation | What it receives and returns |
|---|---|
| `.pair_km()`, `.pair_G()` | Observed end times, death indicators and horizon from one arm yield sorted censoring jumps, risk counts and survival products; evaluating these at the distinct pair cutoffs gives the weights. |
| `.pair_events()` | Validated patient and event tables become treatment-patient histories and a single sorted control event stream. No event after the current look enters the index. |
| `.pair_summaries()`, `.rank_eligible()` | From the validated method and clinical hierarchy, select a calculation for the two component means, each patient's opposite-arm means, ascertainment count and settlement-time masses. |
| `.rank_ipcw()` | For supported restricted ordered tiers, compute those summaries from arm risk counts, observed deaths, censoring weights and ranks of the surviving patients. |
| `.pair_endpoint_spec()`, `.rank_tiers()` | Inspect whether each clinical function remains the endpoint originally constructed and extract its ordered comparison. A free-form or modified function executes its own body. |
| `.rank_channel_key()`, `.rank_tier_counts()` | Build one score per eligible patient and refine tied ranks only where both arms remain. Ascertainment and censoring weights are handled separately. |
| `.rank_margin_counts()` | Within each group tied on earlier tiers, count opposite-arm scores strictly outside a final difference margin by two binary searches. |
| `.pair_summaries_reference()`, `.pair_endpoint_totals()` | Score each pair for shared-window comparisons and the other restricted hierarchies. Assign its time-specific contribution once, using the treatment endpoint when both times coincide. |
| `.if_shared_win_loss()` | Shared-window summaries give the unweighted signed row/column projection. |
| `.if_ipcw_benefit()` | The restricted summaries and arm KM fits give the signed influence; requesting log odds applies the scalar derivative to the same result. |
| `.if_ipcw_win_ratio()` | The same restricted summaries and KM fits give the bivariate delta-method influence on the log scale. |
| `.if_ipcw_linear()`, `.pair_weight_if()` | The contrast coefficients turn common win/loss aggregates into one scalar influence. A reverse cumulative sum computes H_v for every relevant censoring jump; a forward cumulative sum of hazard residual terms evaluates every patient without materializing a patient-by-pair or patient-by-jump array. |
| `.pair_fit()` | Selects the authorized method/effect, checks follow-up support, returns estimate, natural-scale estimate, influences in input row order, and diagnostics (ascertainment, means, risk and survival). |

For the built-in restricted hierarchy, let R_a(s) count patients observed at least through s and D_a(s) count documented deaths at s in arm a. At each s through the horizon the death comparison gives D_0(s)[R_1(s) − D_1(s)] wins, D_1(s)[R_0(s) − D_0(s)] losses, and D_0(s)D_1(s) simultaneous-death ties. A patient censored exactly at s is included in the opposite arm's risk count. Multiplying by the inverse censoring survival at s gives the weighted component totals. Prefix sums of opposite-arm deaths give the contribution for each patient. Patients still observed alive through the horizon are ordered by admission count and, within a count, last admission time. The numbers of opposite-arm patients with better and worse ranks supply the remaining losses and wins. These pairs all settle at the horizon. A declared event count, first or last event time, measurement value, or threshold category can refine ranks within a group tied on earlier tiers. A positive terminal difference margin uses two opposite-arm range counts within each group. A nonterminal positive margin or arbitrary callback uses direct scoring. A missing deciding value also uses direct scoring to identify the original deferral reason. For the default hierarchy with n patients and e event rows, the sort and accumulation use O(n log n + e log e) time and O(n + e) working storage.

The general scorer compares all n_1 n_0 pairs while holding one row of pairs in memory. It handles shared windows, clinical rules requiring direct comparison, and the count-only default in an older saved design. Both calculations return two patient means and both grid masses to the same censoring influence and effect routines.

To calculate covariance between two performed looks, match by arm and patient ID and add the products of patient influences divided by the product of the two arm enrollment counts. Use the influence and test scale for the chosen method at both looks. The calculations here do not assume independent increments and do not jointly monitor multiple effects.

## Reporting three descriptive effects from one fit

`gpc_ipcw_effects(plan, look, patients, events)` validates the cut, obtains its horizon, then calls `.pair_core()` once. It fits each arm's reverse-KM curve once and computes the two sets of win and loss summaries together. `.pair_effect()` derives a requested contrast from that common core. For all three effects together, net benefit supplies the signed influence once; log win odds scales that same vector. Log win ratio applies a second scalar combination to stored win and loss components. Normal monitoring still computes only the selected effect.

`summary$se` and the three patient influence columns are on the test scale: ordinary net benefit, log win odds, and log win ratio. A constant but estimable quantity may have zero variance descriptively; the monitoring layer cannot standardize it. A log win-odds estimate is deferred when the estimated benefit is indistinguishable from ±1 at working precision (within 32 double-precision machine epsilons), rather than converting a unanimous result into a spurious large finite log odds by roundoff. A failed transform gives one invalid row and missing subjects' contributions for that transform; a common support failure affects the whole look. These are descriptive alternatives to report, not three tests that can each trigger stopping.
