"""Bounded, post-result CTQ neglect functional-form diagnostics.

Reads de-identified participant-level input; writes aggregate tables only.
See CTQ_MODEL_FORM_PLAN_20261020.md, copied from the plan saved before fitting.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
from patsy import build_design_matrices, dmatrix
import scipy
import statsmodels


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_KNOTS = ROOT/"source_data/diagnostic_20261020/spline_knots_frozen.json"
DEFAULT_SCALING = ROOT/"source_data/diagnostic_20261020/manuscript_scaling_frozen.json"
DEFAULT_FROZEN = ROOT/"source_data/fig3_conditional_difference_curves.csv"
PLAN = HERE/"CTQ_MODEL_FORM_PLAN_20261020.md"
EXPECTED_SOURCE_SHA = "5342bd167496989d0ebf083894d24547c50b210356715cdeceb49b76b5ad352b"
SCORES = (17, 25, 32)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve_upstream_code(explicit: Path | None) -> Path:
    candidates = ([explicit] if explicit is not None else []) + [
        HERE/"upstream_analysis/post_review_20261015/post_review_analysis.py",
        ROOT/"archive/原临床分析代码/post_review_20261015/post_review_analysis.py",
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError("The upstream post_review_analysis.py is missing; pass --upstream-code.")


def load_old_module(path: Path):
    spec = importlib.util.spec_from_file_location("ctq_fixed_prior_analysis", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def spline_term(name: str, knots: dict) -> str:
    k = knots[name]
    return (
        f"cr({name}, knots=({k['inner'][0]!r},{k['inner'][1]!r}), "
        f"lower_bound={k['lower']!r}, upper_bound={k['upper']!r}, "
        "constraints='center')"
    )


def build_formula(model_id: str, knots: dict) -> str:
    shared = [spline_term("abuse_z", knots), spline_term("distress_z", knots),
              "C(grade)*sex", "C(site)"]
    if model_id in ("original_spline", "person_quantile_spline"):
        exposure = spline_term("neglect_z", knots)
    elif model_id == "minimum_score_pileup":
        exposure = "above_floor+excess_points"
    else:
        raise ValueError(model_id)
    return "+".join([exposure] + shared)


def difference_vector(info, raw_columns, kept_columns, frame, model_id, score, scales, prior):
    if model_id == "minimum_score_pileup":
        low = frame.copy()
        high = frame.copy()
        low["above_floor"] = 0
        low["excess_points"] = 0.0
        high["above_floor"] = int(score > 10)
        high["excess_points"] = float(score - 10)
        low_x, high_x = build_design_matrices([info], low)[0], build_design_matrices([info], high)[0]
        whole = np.asarray(high_x).mean(axis=0) - np.asarray(low_x).mean(axis=0)
    else:
        mu, sd = scales["neglect"]["mean"], scales["neglect"]["sd"]
        whole = prior.design_diff(info, frame, "neglect_z", (10-mu)/sd, (score-mu)/sd)
    positions = [raw_columns.get_loc(name) for name in kept_columns]
    return np.asarray(whole)[positions]


def fit_one(prior, frame, model_id, knots, scales):
    formula = build_formula(model_id, knots)
    raw_x = dmatrix(formula, frame, return_type="dataframe")
    info = raw_x.design_info
    x, fit = prior.fit_matrix(raw_x, frame.outcome, frame.site, model_id)
    vectors = {score: difference_vector(info, raw_x.columns, x.columns, frame,
                                        model_id, score, scales, prior)
               for score in range(10, 33)}
    result = []
    for score, vector in vectors.items():
        c = prior.clustered_contrast(fit, x, frame.site, vector, model_id)
        result.append(dict(model=model_id, neglect_score=score, reference_score=10,
                           **c))
    assert all(abs(c) < 1e-12 for c in [result[0]["estimate"], result[0]["se_school"]])
    model_info = dict(model=model_id, formula=formula, n=len(frame),
                      n_schools=int(frame.site.nunique()), parameters=int(x.shape[1]),
                      r_squared=float(fit.rsquared))
    return result, model_info, (raw_x, info, x, fit)


def run(args):
    out = args.output.resolve()
    for name in ("data","knots","scaling","frozen"):
        if not getattr(args,name).is_file():
            raise FileNotFoundError(f"Missing --{name} input: {getattr(args,name)}")
    upstream_code = resolve_upstream_code(args.upstream_code)
    assert sha256(args.data) == EXPECTED_SOURCE_SHA, "The participant table differs from the fixed source."
    out.mkdir(parents=True, exist_ok=False)
    prior = load_old_module(upstream_code)
    prior.SOURCE = args.data
    prior.prior.SOURCE = args.data
    prior.prior.SCALING_REFERENCE = args.scaling
    baseline, _, d, scales, scale_error = prior.prior.prepare()
    assert len(baseline) == 8902 and len(d) == 4219 and d.site.nunique() == 10
    assert scale_error < 1e-10
    assert tuple(float(baseline.neglect.quantile(q)) for q in (.25, .5, .75, .9, .95)) == (10.,12.,17.,23.,32.)
    frame = d.copy()
    frame["outcome"] = frame.distress_future_z
    frame["above_floor"] = frame.neglect.gt(10).astype(int)
    frame["excess_points"] = np.maximum(frame.neglect-10., 0.)

    original_knots = json.loads(args.knots.read_text())
    assert np.allclose(original_knots["neglect_z"]["raw_inner"], [21.333333333333332,32.666666666666664])
    person_knots = json.loads(args.knots.read_text())
    mu, sd = scales["neglect"]["mean"], scales["neglect"]["sd"]
    person_knots["neglect_z"]["inner"] = [(score-mu)/sd for score in (12.,23.)]
    person_knots["neglect_z"]["raw_inner"] = [12.,23.]

    all_curves, model_records, fits = [], [], {}
    for model_id, knot_set in [
        ("original_spline", original_knots),
        ("person_quantile_spline", person_knots),
        ("minimum_score_pileup", original_knots),
    ]:
        rows, meta, fit_data = fit_one(prior, frame, model_id, knot_set, scales)
        all_curves.extend(rows)
        model_records.append(meta)
        fits[model_id] = fit_data

    curves = pd.DataFrame(all_curves)
    frozen = pd.read_csv(args.frozen)
    frozen = frozen[frozen.specification.eq("exposure_and_baseline_splines")]
    max_error = 0.
    for score in SCORES:
        new = curves[curves.model.eq("original_spline") & curves.neglect_score.eq(score)].iloc[0]
        old = frozen[frozen.neglect_score.eq(score)].iloc[0]
        for field in ("estimate", "se_school", "school_low", "school_high"):
            max_error = max(max_error, abs(float(new[field])-float(old[field])))
    assert max_error < 1e-10, f"Original model did not reproduce frozen results: {max_error}"

    # School deletion keeps the original whole-sample spline basis and scales fixed.
    raw_x, info, _, _ = fits["original_spline"]
    loo = []
    for index, school in enumerate(sorted(frame.site.unique()), 1):
        include = frame.site.ne(school)
        reduced, reduced_fit = prior.fit_matrix(raw_x.loc[include], frame.loc[include,"outcome"],
                                                frame.loc[include,"site"], f"original_spline_LOSO_{index:02d}")
        for score in SCORES:
            vector = difference_vector(info, raw_x.columns, reduced.columns, frame,
                                       "original_spline", score, scales, prior)
            c = prior.clustered_contrast(reduced_fit, reduced, frame.loc[include,"site"],
                                         vector, f"LOSO_{index:02d}")
            loo.append(dict(omitted_school=f"school_{index:02d}", neglect_score=score,
                            reference_score=10, estimate=c["estimate"],
                            remaining_participants=int(include.sum()),
                            remaining_schools=int(frame.loc[include,"site"].nunique())))
    loo_df = pd.DataFrame(loo)
    assert len(loo_df) == 30 and loo_df.remaining_schools.eq(9).all()

    support = []
    for sample_name, table in [("full_wave2_baseline", baseline),("paired_wave2_wave3", frame)]:
        for score in (10,17,25,32):
            exact = table.neglect.eq(score)
            tail = table.neglect.ge(score)
            support.append(dict(sample=sample_name, n=len(table), score=score,
                                exact_n=int(exact.sum()), exact_percent=100*float(exact.mean()),
                                exact_schools=int(table.loc[exact,"site"].nunique()),
                                at_or_above_n=int(tail.sum()),
                                at_or_above_percent=100*float(tail.mean()),
                                at_or_above_schools=int(table.loc[tail,"site"].nunique())))
    support_df = pd.DataFrame(support)

    selected = curves[curves.neglect_score.isin(SCORES)].copy()
    selected = selected[["model","neglect_score","reference_score","estimate","se_school",
                         "school_low","school_high","schools","school_df","n_rows","parameters"]]
    selected.to_csv(out/"fixed_score_contrasts.csv", index=False)
    curves.to_csv(out/"conditional_difference_grid_10_to_32.csv", index=False)
    loo_df.to_csv(out/"original_spline_leave_one_school_out.csv", index=False)
    support_df.to_csv(out/"score_support_aggregate.csv", index=False)
    pd.DataFrame(model_records).to_csv(out/"model_manifest.csv", index=False)

    # A picture of the already specified score grid; no additional fitting or testing.
    fig, (ax, axn) = plt.subplots(2, 1, figsize=(9.4,6.4), sharex=True,
                                  gridspec_kw={"height_ratios":[3.1,1.0]})
    palette = {"original_spline":"#1565C0", "person_quantile_spline":"#BF360C",
               "minimum_score_pileup":"#00695C"}
    labels = {"original_spline":"Original spline (21.3, 32.7)",
              "person_quantile_spline":"Person-quantile spline (12, 23)",
              "minimum_score_pileup":"Floor-pileup model"}
    for model_id in palette:
        table = curves[curves.model.eq(model_id)].sort_values("neglect_score")
        ax.plot(table.neglect_score, table.estimate, label=labels[model_id],
                color=palette[model_id], linewidth=2.1)
        chosen = table[table.neglect_score.isin(SCORES)]
        ax.scatter(chosen.neglect_score, chosen.estimate, color=palette[model_id], s=20)
    linear = pd.read_csv(args.frozen)
    linear = linear[linear.specification.eq("all_linear") & linear.neglect_score.le(32)].sort_values("neglect_score")
    ax.plot(linear.neglect_score, linear.estimate, color="#555555", ls="--", lw=1.6,
            label="Historical all-linear model (other terms also linear)")
    ax.axhline(0, color="#555555", lw=.6)
    for score in SCORES:
        ax.axvline(score, color="#BDBDBD", lw=.7, ls=":")
    ax.set_ylabel("Conditional W3 distress difference\n(full-W2 baseline SD)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(axis="y", color="#E0E0E0", lw=.5)
    frequencies = frame.neglect.value_counts().reindex(range(10,33),fill_value=0)
    axn.bar(frequencies.index, frequencies.values, color="#78909C", width=.8)
    axn.set_ylabel("Paired n\nat exact score")
    axn.set_xlabel("W2 CTQ neglect composite score")
    axn.set_xlim(9.5,32.5)
    axn.set_xticks([10,12,15,17,20,23,25,28,30,32])
    axn.grid(axis="y", color="#E0E0E0", lw=.5)
    fig.tight_layout()
    for suffix in ("png","pdf","svg"):
        fig.savefig(out/f"functional_form_comparison.{suffix}", dpi=220)
    plt.close(fig)

    audit = dict(
        analysis_date="2026-10-20", post_result_exploratory=True,
        complete_sample_n=len(frame), complete_sample_schools=int(frame.site.nunique()),
        full_reference_n=len(baseline), original_fit_max_abs_difference_from_frozen=max_error,
        input_sha256=sha256(args.data), frozen_sha256=sha256(args.frozen),
        knots_sha256=sha256(args.knots), scaling_sha256=sha256(args.scaling),
        plan_sha256=sha256(PLAN), code_sha256=sha256(Path(__file__)),
        upstream_code_sha256=sha256(upstream_code),
        upstream_content_code_sha256=sha256(upstream_code.parent/"dependencies/content_analysis_20261012.py"),
        upstream_common_code_sha256=sha256(upstream_code.parent.parent/"revision_common.py"),
        versions=dict(python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__,
                      scipy=scipy.__version__, statsmodels=statsmodels.__version__,
                      patsy=patsy.__version__, matplotlib=matplotlib.__version__),
        restrictions=["No manuscript or historical results modified.",
                      "All confidence intervals are school-CR1 t(9) pointwise exploratory intervals.",
                      "No p-value or model selection search.",
                      "No participant-level rows, IDs, or exact-score-by-school cells exported."],
    )
    (out/"analysis_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2)+"\n")
    print(selected.to_string(index=False, float_format=lambda x:f"{x:.5f}"))
    print("LOSO 25 vs 10 range:",
          loo_df.loc[loo_df.neglect_score.eq(25),"estimate"].min(),
          loo_df.loc[loo_df.neglect_score.eq(25),"estimate"].max())
    print("Output:",out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True,
                        help="Authorized local de-identified questionnaire_analysis.csv; never copied to output")
    parser.add_argument("--knots", type=Path, default=DEFAULT_KNOTS)
    parser.add_argument("--scaling", type=Path, default=DEFAULT_SCALING)
    parser.add_argument("--frozen", type=Path, default=DEFAULT_FROZEN)
    parser.add_argument("--upstream-code", type=Path, default=None,
                        help="Existing post_review_analysis.py; defaults to bundled or repository archive copy")
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
