# CTQ neglect analysis

Analysis and figure code for **Childhood neglect and subsequent emotional distress across symptom domains in school-aged youth**.

本仓库保存 SEARCH 学校问卷论文的统计分析、图件生成、核查代码和分析计划。代码取自 2026-10-17 复核修订包，原统计与绘图脚本保持不变。仓库不包含参与者问卷、原始数据、录音、汇总结果、论文或图件；创建仓库时未重新拟合模型。

This repository is a code snapshot, not a self-contained data reproduction package. The participant-level questionnaire data and frozen reference inputs remain in the authorized local research package. Repository creation does not grant permission to share those inputs.

## Contents

| Path | Purpose |
| --- | --- |
| `code/build_current_figures.py` | Prepare aggregate figure sources and current figures from an existing local source package; no model fitting |
| `code/export_high_score_support.py` | Export school-level support summaries from authorized local questionnaire data |
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

The scripts in `code/README_reproduction.md` and `code/README_figures.md` describe the original full local package, including inputs deliberately omitted here. In particular, `code/check_reproduction_inputs.py` checks that full package and will report missing inputs in a fresh code-only clone. Use `scripts/check_repository.py` for this repository's data-free checks.

## Interpretation and scope

The current main symptom-content comparisons use common-prior adjustment. Overall distress retains the previously specified linear-versus-spline comparison. Own-prior and earlier-interval content analyses have different adjustment specifications and are not identical replications of the common-prior model.

The historical `post_review_analysis.py` still reads adult aggregate constants to preserve its historical scaling-register output. Those reference inputs and adult rows are not part of the current manuscript's active analyses. See [docs/SCOPE_AND_PROVENANCE.md](docs/SCOPE_AND_PROVENANCE.md) before running or interpreting historical outputs.

This repository does not claim to contain an executable reimplementation of every one of the 26 earlier core-model verification records. Those records describe earlier work in the full local package. The four analysis entry points above and current figure code are the source files supplied in the latest submission package.

## Attribution and access

The figure QA helpers were supplied with the manuscript package and attributed there to the `nature-figure` skill. Their original code is preserved. No blanket open-source licence is assigned here; this snapshot does not assert sole authorship or relicensing rights for those helpers. Review ownership and applicable terms before public redistribution.

The GitHub repository was initially created as private. Data access, software access and any future public release are separate decisions. A private repository URL alone is not a public code-availability statement.
