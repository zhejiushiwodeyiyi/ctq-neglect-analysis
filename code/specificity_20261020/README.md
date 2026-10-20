# Post-result symptom-content checks (20 October 2026)

These scripts implement the dated, **post-result exploratory** plan in
`ANALYSIS_PLAN.md`. The public plan is a path-scrubbed copy of the original
locked document; its original SHA-256 is recorded there. This was not a
prospective registration or an independent replication.

The four analysis and audit Python files are portable CLI adaptations of the scripts actually
executed in the local analysis package. The adaptations replace machine-specific
input and output paths with explicit arguments, refuse to write inside this
public repository, and preserve the locked input hash and frozen-result checks.
They do not change the model formulas, samples, contrasts, item definitions, or
statistical calculations. Using the authorized original inputs, all ten
generated aggregate CSVs were compared field by field with the original run;
all shared numeric fields matched exactly. The output files themselves are
**not** distributed here.

`build_figure_s4.py` is the same figure rendering code with explicit aggregate
input and output paths and the repository's licensed panel-alignment helper.
Its figure-source CSVs and 600-dpi PNG matched the local submission figure
byte for byte in the portability check. No image is stored here.

## Required local inputs

- The authorized participant-level `questionnaire_analysis.csv`, with SHA-256
  `5342bd167496989d0ebf083894d24547c50b210356715cdeceb49b76b5ad352b`.
- Frozen aggregate files: `spline_knots_frozen.json`,
  `manuscript_scaling_frozen.json`, `function_form_iqr_contrasts.csv`,
  `secondary_common_prior_contrasts.csv`, and `ipw_summary.json`.
- Optional cross-check: `measurement_grade_content_descriptives.csv`.
- Optional byte-identical original plan and previous Word files for audit.

These inputs are intentionally absent from this code-only repository. Access to
the participant data and redistribution of any resulting tables are separate
decisions by the data custodian. In particular, `measurement_audit.py` writes
item-category frequencies that can contain small cells; its local output must
not be uploaded or released without disclosure review. The leave-one-school-out
and balance outputs also deserve review before release.

## Run with authorized inputs

The following is a path template. All output directories must be new and
outside this repository. Use the project environment described in the root
README, and run from the repository root:

```sh
python code/specificity_20261020/measurement_audit.py \
  --data /path/to/authorized/questionnaire_analysis.csv \
  --reference-descriptives /path/to/authorized/measurement_grade_content_descriptives.csv \
  --output /path/to/new-local-run/measurement

python code/specificity_20261020/bounded_models.py \
  --data /path/to/authorized/questionnaire_analysis.csv \
  --knots /path/to/authorized/spline_knots_frozen.json \
  --scaling /path/to/authorized/manuscript_scaling_frozen.json \
  --frozen /path/to/authorized/function_form_iqr_contrasts.csv \
  --output /path/to/new-local-run/fractional

python code/specificity_20261020/ancillary_checks.py \
  --data /path/to/authorized/questionnaire_analysis.csv \
  --ipw-summary /path/to/authorized/ipw_summary.json \
  --secondary /path/to/authorized/secondary_common_prior_contrasts.csv \
  --output /path/to/new-local-run/ancillary

python code/specificity_20261020/verify_run.py \
  --run /path/to/new-local-run \
  --data /path/to/authorized/questionnaire_analysis.csv \
  --frozen /path/to/authorized/function_form_iqr_contrasts.csv

python code/specificity_20261020/build_figure_s4.py \
  --model-grid /path/to/new-local-run/fractional/four_model_grid.csv \
  --ancillary /path/to/new-local-run/ancillary/ancillary_contrasts.csv \
  --output /path/to/new-local-figure-output
```

The `--reference-descriptives` and `--original-plan` checks are optional in
the first two scripts. `verify_run.py` also accepts `--original-plan`,
`--original-main-word`, and `--original-supplement-word` to check their
recorded hashes. The original run used all these checks.

`bounded_models.py` reuses the archived common-prior preparation and fixed
spline formula; `ancillary_checks.py` reuses the archived IPW functions. The
frozen inputs are necessary to verify identical scales and prior estimates.
Fractional-logit uncertainty uses school-cluster CR1 covariance and a t(9)
approximation; the early interval has only four schools. The seven-item
measurement model is a limited-information, descriptive polychoric diagnostic,
not a clustered confirmatory factor analysis or invariance test. None of the
scripts establishes a clinically validated symptom-specific effect.
