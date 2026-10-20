# Post-result sensitivity plan: unmatched W3 outcomes within supported schools/grades

Written 2026-10-20 before running the new sensitivity script. This is a **post-result, exploratory bias analysis**, not a preregistered or independent confirmatory test. The original participant data and current manuscript/package remain unchanged.

## Target and missingness boundary

- Source: the frozen private `questionnaire_analysis.csv` (SHA-256 `5342bd167496989d0ebf083894d24547c50b210356715cdeceb49b76b5ad352b`). Read in place; export aggregates only.
- Start with complete W2 CTQ and DASS baseline records in grades 4, 5, 7, 8, 10, 11. A W3 outcome is observed only when the same participant ID has a complete W3 DASS record.
- Restrict the sensitivity target to W2 school-by-grade cells with at least one observed matched W3 outcome. This is the **4,519-person supported baseline domain**, with **4,219 matched** and **300 unmatched**. It is not the full 8,902-person W2 baseline or all 5,941 continuing-grade participants. Entirely unsupported school/grade cells are excluded; no source-code equivalence between W2-only `site_09` and W3-only `site_03` will be assumed.

Pre-fit source-audit amendment (also 2026-10-20): among individually linked records, 228 have W2 school code `site_01` and W3 code `site_03`; the other 3,991 keep the same code. The participant IDs provide individual linkage, but the available file does not establish the reason for the school-code change or that the schools are globally interchangeable. The analysis will retain the **W2 baseline school** for its covariate and clustering, as in the existing paper, and flag the code discrepancy for author verification. No school-level relabeling or extension to `site_09` will be made.

## Anchor regressions and reference means

1. Reconstruct the current complete-case W2→3 general-distress OLS contrast: standardized W3 general distress as outcome; W2 neglect and non-sexual abuse standardized from the **full 8,902 W2 baseline**, W2 general distress, grade × sex, and school indicators. The neglect coefficient per baseline SD is the anchor.
2. Reconstruct the current complete-case raw-item-mean content contrast: W3 anhedonic-content mean minus autonomic-content mean (each on 0–3), regressed on W2 neglect and abuse standardized from the same full baseline, all **three W2 DASS domain scores**, grade × sex and school. The neglect coefficient multiplied by `7 / SD(W2 neglect)` gives the 17-vs-10 contrast. This difference-outcome OLS fit is algebraically equivalent to the current fully interacted stacked common-prior contrast, provided both content outcomes use the same participants and covariates. Verify anchors numerically against frozen published source data before proceeding.
3. Fit **three reference conditional-mean models to observed W3 outcomes only**, each fractional-logit GLM using the same covariates as its anchor: W3 general distress divided by 42; W3 anhedonic and autonomic item means each divided by 3. The link keeps reference predictions in the outcome range. This reference is a model-based MAR-style assumption, **not evidence that missingness is MAR**. Reference predictions are used for bias scenarios, never passed off as observed values.

## One pattern-mixture sensitivity and tipping points

Let `z` be W2 neglect standardized using all 8,902 W2 participants and `z̄_mis` its mean among the 300 unmatched in supported cells. For an unmatched person, set the *hypothetical conditional mean* of W3 distress to `clip(reference_mean + δ_D × (z − z̄_mis), 0, 42)`. `δ_D` is an unobserved, exposure-dependent departure in **raw DASS distress points per 1 SD of baseline neglect**. For content, set the two hypothetical conditional means to `clip(reference_A + δ_C/2 × (z − z̄_mis), 0, 3)` and `clip(reference_U − δ_C/2 × (z − z̄_mis), 0, 3)`; `δ_C` is a departure in the **anhedonic-minus-autonomic item mean per 1 SD of neglect**. The centering avoids deliberately changing the missing group's mean at the mean of its own neglect distribution; boundary clipping can still change its aggregate mean. Observed W3 outcomes are never changed.

For each scenario, regress the observed outcomes plus these **counterfactual conditional means**, not invented individual observations, on the unchanged original design matrix over all 4,519 supported W2 participants. Report the resulting *scenario point estimate only*. This is deterministic bias analysis; ordinary OLS/cluster confidence intervals or p values on these pseudo-completed outcomes would understate uncertainty and will **not** be reported.

Fixed diagnostic grids, chosen before inspecting sensitivity results:

- `δ_D = −42, −41.5, …, +42` raw distress points per baseline-neglect SD; this wide grid is an algebraic stress test, not a claim of plausible departures. Report the local scenarios −8, −4, 0, +4, +8 separately for interpretation.
- `δ_C = −3, −2.975, …, +3` content item-mean units per baseline-neglect SD; again a stress test. Report the local scenarios −0.5, −0.25, 0, +0.25, +0.5 separately.
- Find a negative-departure value at which the **scenario point estimate** reaches zero only if the bounded scenario curve actually crosses zero; otherwise state that no point-estimate tipping point was reached within the grid. Do not infer a p-value tipping point.
- Export aggregate scenario curves, selected-scenario summary, aggregate domain counts, model/anchor diagnostics, and a readme. Do not export individual rows, IDs, conditional predictions, or fitted weights.

## Interpretation limits

These scenarios assess one specified exposure-dependent unobserved-outcome mechanism among the 300 unmatched participants **only**. They cannot validate DASS-21 for younger children, recover outcomes for completely unsupported schools/grades, identify why follow-up is missing, prove MAR/MNAR, or establish generalizability to the full W2 cohort. Results depend on the reference conditional-mean model and its outcome bounds. A scenario point estimate is not a corrected estimate or a new confirmatory model; it describes how much assumed differential missingness would change the existing association.
