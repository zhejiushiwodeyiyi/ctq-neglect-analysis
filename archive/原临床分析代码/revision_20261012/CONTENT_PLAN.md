# CTQ revision: fixed analysis plan

Date: 2026-10-12. This plan was written after peer review and before fitting the new models below. These are post-review exploratory sensitivity analyses, not preregistered analyses or independent replications.

## Inputs and fixed population

Read only `working_package/data/questionnaire_analysis.csv`. Use harmonized DASS responses (0–3); do not subtract 1 again. Baseline: wave 2 with complete source CTQ and DASS flags. Follow-up: wave 3 with complete source DASS flag. Intersect participants, require observed grade, harmonized sex, school and all model values. Expect 4,219 participants in 10 schools. Derive school codes from the supplied coded school field. Keep grade categorical and grade-by-sex interactions. Do not impute responses, alter source values or release participant-level output.

Standardization uses the eligible full wave-2 baseline population, with population SD (`ddof=0`). Apply the same baseline mean/SD to follow-up responses. Broad distress is the mean of DASS depression, anxiety and stress (each 0–42). Neglect and abuse are the sums of their two corresponding CTQ subscales.

## 1. Common-baseline direct content contrasts

Anhedonic content is the mean of DASS items 3, 10 and 16; autonomic content is the mean of items 2, 4, 7 and 19. Fit both follow-up outcomes with exactly the same baseline covariate set: all three baseline DASS domains, neglect, abuse, categorical grade-by-sex, and school fixed effects. Fit a fully interacted stacked model so each nuisance slope may differ between contents; the exposure-by-content coefficient is the direct difference (anhedonic minus autonomic). Report each content slope and the direct contrast. Analyse standardized response units and raw 0–3 response units. Both use an exposure SD from the full baseline population. These unit representations are the same underlying data and are not independent corroboration.

For each units family, apply BH to the two restricted studentized school-wild-bootstrap P values (neglect and abuse). Use all 1,024 Rademacher sign patterns for 10 schools, restricted null fits, CR1 studentization and two-sided exceedance. Report school-CR1 SE and 95% t(9) confidence intervals separately from bootstrap P values. Retain original own-prior-content contrasts as explicitly labelled comparisons; do not overwrite frozen original tables. Require stacked estimates to equal the corresponding two separate OLS slope differences to within 1e-10. No selective alternative baseline models will be searched.

## 2. Mutually adjusted CTQ subscale estimates

Fit follow-up standardized broad distress on all four baseline CTQ subscales together, baseline standardized broad distress, categorical grade-by-sex, and school fixed effects. Express each subscale in its own full-baseline SD. Report all four coefficients, school-CR1/t(9) intervals, two-sided restricted studentized school-wild-bootstrap P values (1,024 patterns) and BH across the four subscale tests. Also rescale all coefficients and intervals to raw distress points per five CTQ-subscale points. Do not infer that subscales differ merely because their P values differ. Do not treat an emotional/physical subscale association as a validated deprivation mechanism or claim a formal difference without a separate direct test.

## 3. Grade-stratified measurement descriptions

Describe paired students according to baseline grade 4–5 versus grade ≥7, at both baseline and follow-up: participant count, school count, raw content mean/SD, zero-score count and percentage, raw Cronbach alpha, average within-content Pearson inter-item correlation and between-content Pearson correlation. Also report baseline CTQ-subscale score distributions (n, mean, sample SD, median, IQR, min/max and score=5 floor) within these paired strata and within the full eligible wave-2 grade strata. No between-grade significance tests, reliability rankings or age-validity claims will be made. Grade is only an available proxy for educational stage; chronological age and actual questionnaire administration remain author-confirmation items. Internal consistency alone cannot validate constructs or establish suitability in younger children.

## Reproducibility and interpretation

Write only aggregate CSV/JSON files and code under `analysis/`. Record source SHA-256, input sample sizes, baseline scaling, versions, full-rank checks, bootstrap pattern counts, and direct-contrast identity checks. Record all planned outcomes and analyses regardless of result. Distinguish these new sensitivity analyses from the original discovery search and limited within-family BH adjustments; do not represent BH as correction of the manuscript's entire exploration.
