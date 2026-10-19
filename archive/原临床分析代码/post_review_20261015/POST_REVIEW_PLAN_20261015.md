# Fixed post-review analysis plan, 15 October 2026

This plan was written before running the new analyses below. It responds to the 14 October 2026 mock reviews, after the investigators/analysts have seen the existing full-sample results. It is neither historical preregistration nor independent confirmation. Existing primary models, units, data and earlier outputs will remain unchanged. Outputs contain aggregates only.

## Inputs and invariant checks

Use the corrected de-identified `questionnaire_analysis.csv` from the 12 October working package, whose DASS items already use 0–3 coding. Do not subtract one again. Reuse the earlier `content_analysis.prepare`, `school_fit` and restricted studentized `wild_matrix` functions without editing their source. Record input hashes, script/dependency versions and source-data boundaries. Require baseline n=8,902; matched W2–W3 n=4,219; >=7-grade n=2,406 in seven schools. All full-sample wave-2 scales remain fixed (ddof=0).

## A. Secondary-grade common-prior content comparisons

Restrict the paired sample to baseline grade >=7. Fit anhedonic and autonomic item means with the same three baseline DASS domain scores, both CTQ composites, grade-by-sex terms and school fixed effects. Fully interact these terms with content in stacked models, retaining the content-by-exposure contrast. Fit standardized outcomes and raw 0–3 outcomes for each of neglect and abuse. Use school CR1 covariance with t(G−1) intervals, and enumerate all 2^7=128 Rademacher school-sign patterns in the existing null-restricted studentized wild bootstrap. These are computational enumerations, not treatment-randomization tests. Report four new p values; define the principal new BH family as all four contrasts, with the two-exposures-within-units BH values additionally labelled for comparison with earlier tables. Do not infer age differences from separate significance decisions.

## B. Scaling register and algebraic identities

Recompute each school reference population from the canonical corrected table: for W1–W2 use eligible complete wave-1 baseline observations; for W2–W3 use eligible complete wave-2 baseline observations. Report interval, reference population, baseline N, variable nonmissing N, mean, SD and ddof=0. Retain those references even in restricted samples. Historical adult reference means/SDs are copied from existing aggregate records, with reference cohort sizes and variable-specific N explicitly marked unavailable where not recorded; do not infer variable-specific N from cohort size. Check raw-to-standardized component and direct-contrast identities against stored original and common-prior results, including joint contrast covariance. This is a computational scaling audit, not new raw-source adjudication.

## C. Function form and influence, fixed limited scope

Evaluate the full W2–W3 matched sample only. Use two model systems: general distress with prior distress, and common-prior anhedonic/autonomic content. Both adjust for both CTQ composites, grade-by-sex and school fixed effects. Compare the linear model with an additive natural cubic regression spline model having three degrees of freedom per continuous covariate: both CTQ composites and either baseline distress or the three baseline DASS domains. To prevent coincident quantile knots caused by floor concentrations, fix two interior knots at the one-third and two-thirds quantiles of each variable's **sorted unique full-wave-2 baseline values**; fix lower/upper boundaries to that same full-baseline support. Patsy `cr(..., knots=..., lower_bound=..., upper_bound=..., constraints='center')` then gives three centered degrees of freedom. Store every knot and boundary. No alternative knot rules, interactions beyond those already specified, alternative degrees of freedom, or search for smaller p values.

For each exposure, predict the conditional difference from its full-W2 Q25 to Q75, holding the other covariates at each matched participant's observed values and averaging; compare identical IQR estimands across linear and spline fits. For additive models the averaged contrast equals the common contrast in the exposure design basis. For content, report the difference between the two outcome-specific IQR contrasts in both standardized and raw units, with fully interacted stacked school covariance and t(G−1) intervals. These are conditional model-based contrasts, not causal intervention effects, and they differ from constant per-SD slopes. No new p-value or BH claims will be based on these diagnostic contrasts. Report predictions outside the outcome scale as diagnostics, without treating floor concentration alone as an OLS violation.

Compute ordinary least-squares leverage and conventional residual/studentized DFBETAs only as diagnostics. For stacked models calculate exact single-participant-deletion changes through the two-row Woodbury formula, because one participant contributes both content rows. Report only aggregate maxima/quantiles/counts and changes relative to the original school SE. Flag conventional h>2p/n and |DFBETA|>2/sqrt(n) descriptively. A severe influence deletion sensitivity will be run only if h>0.5 or |DFBETA|>1 is observed (participant union across content rows); report any deletion counts. Such a fit is diagnostic and cannot replace the primary estimate. No IDs or individual observations will be exported.

## D. Optional subscale contrast

If implementation and time permit, a formal EN−PN contrast may be computed from the mutually adjusted four-subscale general-distress model, using the already fixed full-W2 per-subscale SD units and school inference. This contrasts specified standardized coefficients only; it is not a contribution decomposition. The manuscript will remove causal-driver/contribution-share wording whether or not this optional analysis is run.

## Interpretation

Report the findings regardless of direction or statistical boundary. Unavailable source documents, actual questionnaire versions/administration, ethics, consent, author contribution and approval remain factual author-confirmation matters. Computation cannot close them.
