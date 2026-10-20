# CTQ neglect analysis

Analysis and figure code for **Childhood neglect and subsequent emotional distress across symptom domains in school-aged youth**.

本仓库保存 SEARCH 学校问卷论文的统计分析、图件生成、核查代码和分析计划。历史代码取自 2026-10-17 复核修订包并保持不变；2026-10-20 增加了事后、范围受限的忽视函数形式和症状内容对比分析代码。仓库不包含参与者问卷、原始数据、录音、汇总结果、论文或图件。

This repository is a code snapshot, not a self-contained data reproduction package. The participant-level questionnaire data and frozen reference inputs remain in the authorized local research package. Repository creation does not grant permission to share those inputs.

## Contents

| Path | Purpose |
| --- | --- |
| `code/build_current_figures.py` | Prepare aggregate figure sources and current figures from an existing local source package; no model fitting |
| `code/export_high_score_support.py` | Export school-level support summaries from authorized local questionnaire data |
| `code/run_model_form_sensitivity.py` | Refit three fixed neglect-function specifications and original-spline leave-one-school-out diagnostics; aggregate outputs only |
| `code/CTQ_MODEL_FORM_PLAN_20261020.md` | Post-result exploratory analysis plan saved before those diagnostic fits |
| `code/specificity_20261020/` | Portable CLI adaptations of the seven-item measurement, bounded-mean, ancillary and verification analyses; Figure S4 renderer; path-scrubbed post-result plan |
| `code/figure_qa_tools/` | Figure alignment, PDF text and collision checks |
| `archive/原临床分析代码/revision_20261012/` | Common-prior content, four-subscale and follow-up/IPW analyses and their dated plans |
| `archive/原临床分析代码/post_review_20261015/` | Grade-restricted common-prior comparisons, spline and influence diagnostics, and the diagnostic addendum |
| `archive/原临床分析代码/revision_common.py` | Shared statistical helpers |
| `archive/原临床分析代码/revision_runtime.py` | Output protection for historical entry points |
| `scripts/check_repository.py` | Standard-library-only code, fingerprint and file-scope checks |
| `SOURCE_MANIFEST.json` | Original relative paths and SHA-256 values for unchanged copied source files |

The archived directory name preserves the relative depth required by the scripts. Do not move these files without updating and checking their path resolution.

## Check the code without data

Python 3.12 is the recorded statistical runtime. The following command needs only the Python standard library and does not import analysis modules, fit models, install software or run bootstrap tests:

```sh
python3 scripts/check_repository.py
```

It checks syntax, the copied-source fingerprints, required local code dependencies and the absence of dataset or binary file types. It does not certify numerical reproduction.

## Install dependencies for local analysis

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

`requirements.txt` lists direct statistical and plotting dependencies without claiming a tested lock file. Historical versions and the limits of the environment record are explained in [docs/SOFTWARE.md](docs/SOFTWARE.md).

## Run analyses with authorized local inputs

The unchanged entry points expect the following inputs at their original relative locations:

- `data/questionnaire_analysis.csv`: restricted, corrected participant-level questionnaire data.
- `archive/原临床分析代码/revision_20261012/reference_results/`: frozen aggregate comparison and scaling inputs.
- `archive/原临床分析代码/post_review_20261015/reference_results/`: frozen aggregate comparison and scaling inputs.
- For the addendum, an explicitly supplied local initial-results directory containing `spline_knots.json` and `function_form_iqr_contrasts.csv`.

These inputs are intentionally absent and ignored by Git. Obtain authorization before making a local working copy. Do not infer permission from possession of a de-identified file. Keep generated outputs outside this repository.

After placing authorized local inputs, run from the repository root and choose a new, previously nonexistent output directory outside the repository:

```sh
python3 archive/原临床分析代码/revision_20261012/content_analysis.py \
  --output /path/to/new-content-output

python3 archive/原临床分析代码/revision_20261012/ipw_analysis.py \
  --output /path/to/new-ipw-output

python3 archive/原临床分析代码/post_review_20261015/post_review_analysis.py \
  --output /path/to/new-post-review-output

python3 archive/原临床分析代码/post_review_20261015/diagnostic_addendum.py \
  --initial-results /path/to/authorized-initial-results \
  --output /path/to/new-addendum-output
```

For the current figure exports, supply the local original source package containing the aggregate tables and unchanged figures:

```sh
python3 code/build_current_figures.py \
  --source-package /path/to/authorized-original-source-package \
  --output /path/to/new-figure-output \
  --qa /path/to/new-figure-qa \
  --audit-scripts code/figure_qa_tools
```

The later model-form diagnostic uses the same corrected 22,265-row questionnaire table and retains the original 4,219-person paired sample. It also requires three authorized **aggregate** inputs: the historical `spline_knots.json`, `manuscript_scaling.json`, and the frozen `fig3_conditional_difference_curves.csv`. These inputs are not in this public code snapshot. Pass their local paths explicitly; the script checks the participant-table SHA-256 and reproduces the frozen original-spline contrasts before writing new aggregate outputs. The `--output` directory must not already exist:

```sh
python3 code/run_model_form_sensitivity.py \
  --data /path/to/authorized/questionnaire_analysis.csv \
  --knots /path/to/authorized/spline_knots.json \
  --scaling /path/to/authorized/manuscript_scaling.json \
  --frozen /path/to/authorized/fig3_conditional_difference_curves.csv \
  --output /path/to/new-model-form-output
```

This exploratory check was specified after the original model discrepancy was known. It does not select a preferred model by significance, validate a clinical threshold, or independently reproduce questionnaire collection and scoring.

The subsequent symptom-content check is also post-result exploratory. It uses
the same fixed 4,219-student paired sample for its four model-form cells,
with school-cluster intervals, and separately reports measurement, early
interval, grade, school-influence and observed-support IPW diagnostics. The
script-by-script input list, commands, portability validation and limits are
in [code/specificity_20261020/README.md](code/specificity_20261020/README.md).
Its local measurement output contains small item-category counts; those
outputs, Figure S4 and its source-data CSVs are deliberately absent here.

The scripts in `code/README_reproduction.md` and `code/README_figures.md` describe the original full local package, including inputs deliberately omitted here. In particular, `code/check_reproduction_inputs.py` checks that full package and will report missing inputs in a fresh code-only clone. Use `scripts/check_repository.py` for this repository's data-free checks.

## Interpretation and scope

The current main symptom-content comparisons use common-prior adjustment. Overall distress retains the previously specified linear-versus-spline comparison; the October 20 checks supplement it without replacing its central 17-versus-10 contrast. Own-prior and earlier-interval content analyses have different adjustment specifications and are not identical replications of the common-prior model.

The historical `post_review_analysis.py` still reads adult aggregate constants to preserve its historical scaling-register output. Those reference inputs and adult rows are not part of the current manuscript's active analyses. See [docs/SCOPE_AND_PROVENANCE.md](docs/SCOPE_AND_PROVENANCE.md) before running or interpreting historical outputs.

This repository does not claim to contain an executable reimplementation of every one of the 26 earlier core-model verification records. Those records describe earlier work in the full local package. The four analysis entry points above and current figure code are the source files supplied in the latest submission package.

## Attribution and access

The four figure QA helpers are unchanged copies of files in `Yuan1z0825/nature-skills` at commit `8630593a086755d2003cbd1a79cd1118c6963202`, distributed under Apache-2.0. Their complete upstream licence is included in [LICENSES/nature-figure-Apache-2.0.txt](LICENSES/nature-figure-Apache-2.0.txt); file-level provenance is documented in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Original source notices are preserved. This licence applies to those helpers; no blanket licence is assigned to the project-specific research code.

The repository is publicly accessible at https://github.com/zhejiushiwodeyiyi/ctq-neglect-analysis following the author's release instruction on 19 October 2026. Participant-level data and frozen result inputs remain excluded. Public code access does not grant access to research data or establish that the omitted inputs may be shared.
