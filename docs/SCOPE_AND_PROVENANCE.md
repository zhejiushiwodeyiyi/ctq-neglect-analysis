# Source scope and provenance

This repository was prepared on 19 October 2026 from the code supplied in the `CTQ_复核修订_20261017` manuscript package. The later removal and restoration of author-input markers changed manuscript text only; they did not replace these scientific scripts.

`SOURCE_MANIFEST.json` records hashes of byte-identical copied Python scripts, analysis plans and original code documentation. The version-metadata JSON was adapted only to omit local executable and module paths. The repository adds its own README, dependency list, Git ignore rules and a standard-library-only file check. It does not change statistical formulas, observations, fitted estimates or interpretation.

## Current and historical analysis outputs

Common-prior models supply the current main symptom-content comparisons. Overall-distress linear and spline comparisons retain their existing covariate definitions. Own-prior, early-interval and school-stage results have their explicitly documented roles. The dated local plans record an exploratory sequence and are not public prospective registrations.

The archived `post_review_analysis.py::scale_register()` also reads `external_clinical_scaling.json` and `external_clinical_sample.csv` from its local `reference_results` directory. They contain historical aggregate adult constants, not newly fitted adult models in that function. Those files are deliberately not in this code-only repository. Their historical output must not overwrite the current school-only 24-row scaling register or be presented as an active adult analysis in the manuscript.

The statistical reference folders, model estimates, figure-source CSVs, manuscript files, figures, raw exports, participant data and recordings remain outside this repository. Current figure generation needs the original authorized source package. Reusing code does not verify the underlying questionnaire wording, private-source reconciliation, ethical approval, consent or data-sharing permission.

## Source files with limited roles

- `code/check_reproduction_inputs.py` is the unchanged preflight for the complete local manuscript package. It intentionally expects files absent here and is not the code-only repository check.
- `post_review_20261015/dependencies/content_analysis_20261012.py` is loaded by the post-review entry point and is not a separate command-line workflow.
- `post_review_20261015/build_diagnostic_figure.py` is a historical figure helper; the current submission figure entry point is `code/build_current_figures.py`.
- The four scripts in `code/figure_qa_tools/` are byte-identical to public upstream commit `8630593a086755d2003cbd1a79cd1118c6963202` of `Yuan1z0825/nature-skills`. Before the public release, their Apache-2.0 licence was verified and its full text added under `LICENSES/`, with fixed source links and hashes in `THIRD_PARTY_NOTICES.md`. They are preserved rather than relabelled as original project code. No repository-wide licence has been inferred from that component licence.

## Verification performed for repository preparation

Only code syntax, copied-file hashes, required code paths and upload file scope were checked. No analysis module was imported, no participant data were read into a statistical pipeline, no model was refitted and no bootstrap calculation was run.

## Subsequent 20 October content-comparison release

`code/specificity_20261020/` was added after the code-only repository was first
prepared. Its five Python files are **adaptations**, not entries in the
byte-identical `SOURCE_MANIFEST.json`: the four executed analytical/checking
scripts were converted from machine-specific paths to explicit CLI inputs,
and the Figure S4 renderer now reads aggregate inputs via CLI and imports
the already attributed panel-alignment helper from this repository. The
path-scrubbed analysis plan preserves the original post-result decisions but
is not byte-identical to the local plan. The original locked plan hash and
input checks remain documented in the code and subdirectory README.

The portable scripts were run locally with the authorized original inputs,
writing only to a temporary directory outside the repository. All ten shared
aggregate CSV outputs matched the original executed run field by field,
including every shared numeric value. The Figure S4 source CSVs and 600-dpi
PNG matched the submission-package files byte for byte. `verify_run.py`
passed with the original plan and prior Word hashes checked. This validates
the path adaptation against that run; it does not independently confirm the
clinical construct or authorize public release of the underlying inputs,
figures, or small-cell aggregate outputs. None of those outputs was added to
the repository.
