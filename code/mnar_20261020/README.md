# Supported-domain missing-outcome scenarios

This is one **post-result exploratory** pattern-mixture sensitivity, not a preregistered analysis or a missing-data correction. The dated `ANALYSIS_PLAN.md` was written before its new fits, after earlier manuscript results were available.

The script requires authorized access to the private `questionnaire_analysis.csv` and checks its SHA-256 before analysis. Run it with an explicit input and an output folder outside this repository:

```sh
python code/mnar_20261020/run_mnar_sensitivity.py \
  --source /path/to/authorized/questionnaire_analysis.csv \
  --output /path/to/new/aggregate-results
```

The checked run used 4,519 baseline students in school-by-grade cells with at least one matched wave-3 outcome: 4,219 matched and 300 without a wave-3 record. The script reproduces the complete-case general-distress and common-prior content point estimates before assigning hypothetical bounded conditional means to the 300 unmatched students. Its outputs contain only aggregate scenario curves, selected assumptions, zero-crossing diagnostics, grade counts, and a model audit. No participant records, predicted row-level values, or fitted weights are written.

Scenario outputs are **point estimates under unverifiable assumptions**, with no confidence intervals or P values. They cannot recover whole school-by-grade groups with zero matched outcomes or establish why records are missing. The resulting aggregate files are not part of this public code repository; consult the manuscript supplement and the authorized local source-data package.

`code/build_figureS5_followup_flow.py` renders the supplementary flow graphic from the aggregate `source_data/figS5_followup_flow.csv` in the authorized local manuscript package. The CSV contains counts already reported in the article but is not included in this code-only repository.
