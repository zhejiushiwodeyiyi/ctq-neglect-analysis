# Post-review observation-weighting analysis plan

Date: 2026-10-12. Written before fitting the additional models. These are exploratory analyses added in response to review, not preregistered analyses or a reconstruction of a historical plan.

## Population and alignment

Use the harmonized, deidentified questionnaire table in `working_package/data/questionnaire_analysis.csv`. Define the baseline as wave 2, and observation as availability of a wave-3 questionnaire for the same random participant code. Check uniqueness of participant code × wave, harmonized DASS reconstruction and source-completeness flags. Use grades 4, 5, 7, 8, 10 and 11 and explicitly mapped sex only; grades 6, 9 and 12 have no observed follow-up and are outside the weighting target. Standardize baseline covariates using all wave-2 records and population SD (ddof = 0), matching the supplied outcome models. Preserve original source data and frozen results.

## Observation model

First refit the exact visible historical specification: `LogisticRegression(C=1, max_iter=2000)` on a Patsy matrix containing `C(grade)*sex + C(site) + abuse_z + neglect_z + dep_score_z + anxiety_score_z + stress_score_z`. Record resolved parameters, installed versions, warning messages, iteration count, rank/columns, retained vs lost probabilities, and agreement with frozen reported weight summaries. Retain the Patsy intercept plus the sklearn intercept in this reproduction to avoid an undisclosed change to the historical penalized parameterization.

For the revised sensitivity, keep this same design matrix, L2 penalty and C = 1; explicitly use lbfgs, fit_intercept=True, max_iter=10000 and tol=1e-10. If numerical convergence is still not achieved, report the limitation and do not silently substitute a selected model. Check the first-order gradient and compare predictions from a second start/solver-independent penalized likelihood optimization before using weights. Do not tune penalty or select specifications based on outcome coefficients.

## Weights and diagnostics

For observed participants use unstabilized weights 1/p and stabilized weights r/p, where r is the observation fraction in the eligible baseline population. Preserve raw weights. Define the primary revised sensitivity using an upper cap at the empirical 99th percentile of observed-participant stabilized weights, as in the supplied historical code. Report raw and truncated distributions, truncation count, sum, effective sample size, baseline/observed/weighted grade × school composition, observed vs mean predicted follow-up probability in each grade × school cell, and retained-versus-lost probability summaries. Empty or zero-observation cells remain visible and are not claimed recoverable by weighting.

Report baseline-target standardized mean differences before and after weighting for sex; continuous grade, neglect, abuse and three DASS domains; all grade and school indicators; and nonempty grade × school cells. Use the eligible baseline population SD as the fixed denominator for continuous variables and sqrt(p_target(1-p_target)) for binary indicators. Signed SMD = (observed or weighted mean − target mean)/denominator. Provide maxima and the share of listed indicators above |SMD| = 0.1; this threshold is a diagnostic description, not proof that selection is removed.

## Outcome sensitivity

Within the same paired eligible sample, compare unweighted, raw-stabilized weighted and 99%-capped weighted regressions of each of the seven prespecified-in-the-supplied-code follow-up outcomes (distress, depression, anxiety, stress, anhedonic content, autonomic content, negative-self content) on both combined adversity scores, baseline depression, anxiety and stress, grade × sex and school. Fit the identical covariate matrix for each weight specification. Report neglect and abuse coefficients, school-cluster CR1 standard errors and t_(G−1) 95% intervals, school-cluster p values, N and rank. Also refit the supplied own-prior-symptom unweighted distress model to show that differences from the principal estimate include a covariate-specification change. No new outcome hypothesis is selected from these sensitivities. No bootstrap/BH claim will be attached to these new weighted estimates; they are descriptive robustness analyses.

## Verification and inference boundary

Independently check sample counts, reindexing of weights to the participant-indexed outcome frame, weight sums and ESS, WLS invariance to multiplying weights by r, and school CR1 covariance against a direct weighted sandwich calculation. Intervals treat fitted weights as fixed and do not incorporate uncertainty in the observation model; only ten observed schools support the broad estimates. A model-based sensitivity can address measured follow-up selection among continuing grades only. It cannot identify missing outcomes in absent grades or cells, establish missing-at-random, restore unmeasured covariates or establish a causal effect.

## Outputs

Create only aggregate CSV/JSON and review text under `analysis/results/ipw*` and `analysis/ipw_revision_text.md`. Do not export participant identifiers, individual predictions or individual weights.

## Post-diagnostic amendment (2026-10-12)

The first diagnostic run found two continuing-grade school × grade cells with zero observed follow-up, both in coded school `site_09` (grades 10 and 11, 1,422 baseline records). The original 5,941-person target is therefore unsupported in these cells. Add a descriptive comparison to the supported-cell baseline target (4,519 records) by removing these two zero-observation cells from the target used for balance and grade composition. Apply the same already fitted probabilities and weights; do not refit or select a new observation model, do not change the 4,219-person outcome sample or the principal unweighted analysis. Report both targets and explain that acceptable measured balance in the supported target cannot establish generalization to the absent school. This amendment is explicitly after seeing the support diagnostics, and is not preregistration.
