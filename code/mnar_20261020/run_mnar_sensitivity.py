"""Aggregate-only, post-result pattern-mixture/tipping-point sensitivity.

See ANALYSIS_PLAN.md, which was written before this script was executed.
This script never writes row-level data, identifiers, predictions, or weights.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
import patsy
from patsy import dmatrix
from scipy import linalg, optimize
from scipy.special import expit
import statsmodels
import statsmodels.api as sm

EXPECTED_SHA = "5342bd167496989d0ebf083894d24547c50b210356715cdeceb49b76b5ad352b"
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "results"
SUB = [
    "ctq_affective_ignore_score", "ctq_body_ignore_score",
    "ctq_affective_abuse_score", "ctq_body_abuse_score",
]
DOMAINS = ["dep_score", "anxiety_score", "stress_score"]
CONTINUING = [4, 5, 7, 8, 10, 11]
DISTRESS_FORMULA = "neglect_z+abuse_z+distress_z+C(grade)*sex+C(site)"
CONTENT_FORMULA = (
    "neglect_z+abuse_z+dep_score_z+anxiety_score_z+stress_score_z"
    "+C(grade)*sex+C(site)"
)
FROZEN_DISTRESS_BETA = 0.16415519467560802
FROZEN_CONTENT_BETA_PER_SD = 0.027584846938074684
FROZEN_CONTENT_DELTA_17_VS_10 = 0.02864053969648072


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    source_hash = digest(source)
    if source_hash != EXPECTED_SHA:
        raise ValueError("Private source SHA-256 differs from frozen source; stop rather than silently reanalyse")
    m = pd.read_csv(source, low_memory=False)
    if len(m) != 22265 or m.duplicated(["participant_code", "wave"]).any():
        raise ValueError("Unexpected participant/wave structure")
    m["sex"] = m.gender.astype("string").str.strip().str.lower().map(
        {"female": 1.0, "male": 0.0, "女": 1.0, "男": 0.0}
    )
    m["site"] = m.school_code.astype(str)
    m["neglect"] = m[SUB[:2]].sum(axis=1, min_count=2)
    m["abuse"] = m[SUB[2:]].sum(axis=1, min_count=2)
    m["distress"] = m[DOMAINS].mean(axis=1)
    m["anhedonic"] = m[["dass3", "dass10", "dass16"]].sum(axis=1, min_count=3) / 3
    m["autonomic"] = m[["dass2", "dass4", "dass7", "dass19"]].sum(axis=1, min_count=4) / 4

    baseline = m.loc[
        m.wave.eq(2) & m.ctq_source_complete & m.dass_source_complete
    ].set_index("participant_code")
    future = m.loc[m.wave.eq(3) & m.dass_source_complete].set_index("participant_code")
    if len(baseline) != 8902 or baseline.sex.isna().any():
        raise ValueError("Expected 8,902 fully measured W2 baseline records with coded sex")
    scales = {}
    for name in ["neglect", "abuse", "distress", *DOMAINS]:
        scales[name] = {
            "mean": float(baseline[name].mean()),
            "sd": float(baseline[name].std(ddof=0)),
        }
        baseline[name + "_z"] = (baseline[name] - scales[name]["mean"]) / scales[name]["sd"]

    eligible = baseline.loc[baseline.grade.isin(CONTINUING)].copy()
    eligible["observed"] = eligible.index.isin(future.index)
    if len(eligible) != 5941 or int(eligible.observed.sum()) != 4219:
        raise ValueError("Unexpected continuing-grade/follow-up counts")
    by_cell = eligible.groupby(["grade", "site"], observed=True)["observed"].sum()
    supported_cells = set(by_cell[by_cell.gt(0)].index)
    eligible["supported"] = [
        (grade, site) in supported_cells for grade, site in zip(eligible.grade, eligible.site)
    ]
    supported = eligible.loc[eligible.supported].copy()
    if len(supported) != 4519 or int(supported.observed.sum()) != 4219 or supported.site.nunique() != 10:
        raise ValueError("Expected 4,519 supported records (4,219 observed, 300 unmatched) in 10 schools")
    observed_ids = supported.index[supported.observed]
    w2_site = supported.loc[observed_ids, "site"].to_numpy()
    w3_site = future.loc[observed_ids, "site"].to_numpy()
    mismatch = w2_site != w3_site
    mismatch_pairs = pd.DataFrame({"w2_site": w2_site[mismatch], "w3_site": w3_site[mismatch]})
    discrepancy = mismatch_pairs.groupby(["w2_site", "w3_site"]).size().to_dict()
    # Individual-ID linkage is real, but a source-code change is NOT evidence of a
    # global school-code mapping. The original models use baseline/W2 school.
    if discrepancy != {("site_01", "site_03"): 228}:
        raise ValueError(f"Unexpected W2/W3 site-code discrepancies: {discrepancy}")
    for name in ["distress", "anhedonic", "autonomic"]:
        supported.loc[observed_ids, name + "_future"] = future.loc[observed_ids, name].to_numpy(float)
    supported["content_gap_future"] = supported["anhedonic_future"] - supported["autonomic_future"]
    observed = supported.loc[observed_ids].copy()
    if (not observed["distress_future"].between(0, 42).all() or
            not observed["anhedonic_future"].between(0, 3).all() or
            not observed["autonomic_future"].between(0, 3).all()):
        raise ValueError("Observed DASS outcome outside instrument range")
    # Private row-level frames remain in memory and are never serialized.
    metadata = {
        "source_sha256": source_hash,
        "baseline_n": len(baseline),
        "continuing_grade_baseline_n": len(eligible),
        "zero_W3_support_continuing_grade_n": int((~eligible.supported).sum()),
        "supported_baseline_n": len(supported),
        "supported_matched_n": len(observed),
        "supported_unmatched_n": int((~supported.observed).sum()),
        "supported_schools": int(supported.site.nunique()),
        "paired_w2_w3_site_code_mismatch_n": int(mismatch.sum()),
        "paired_site_code_mismatch_aggregate": [
            {"w2_site": left, "w3_site": right, "n": int(n)}
            for (left, right), n in discrepancy.items()
        ],
        "full_baseline_scaling": scales,
    }
    return supported, observed, metadata


def ranked_design(supported: pd.DataFrame, observed: pd.DataFrame, formula: str):
    x_full = dmatrix(formula, supported, return_type="dataframe")
    x_obs = x_full.loc[observed.index]
    a = x_obs.to_numpy(float)
    _, r, piv = linalg.qr(a, mode="economic", pivoting=True)
    tol = max(a.shape) * np.finfo(float).eps * np.abs(r).max()
    rank = int((np.abs(np.diag(r)) > tol).sum())
    keep = sorted(piv[:rank])
    if "neglect_z" not in x_full.columns[keep]:
        raise ValueError("Exposure coefficient aliased")
    x_full = x_full.iloc[:, keep]
    x_obs = x_full.loc[observed.index]
    if np.linalg.matrix_rank(x_full.to_numpy(float)) != x_full.shape[1]:
        raise ValueError("Supported-domain model design lost rank")
    return x_full, x_obs


def reference_glm(x_obs: pd.DataFrame, y: pd.Series, outcome_max: float):
    fit = sm.GLM(
        y.to_numpy(float) / outcome_max,
        x_obs,
        family=sm.families.Binomial(link=sm.families.links.Logit()),
    ).fit(maxiter=250, tol=1e-10)
    if not fit.converged or not np.isfinite(fit.params.to_numpy()).all():
        raise RuntimeError("Reference fractional-logit model failed to converge")
    return fit


def coef_projection(x: pd.DataFrame, target: str = "neglect_z") -> np.ndarray:
    a = x.to_numpy(float)
    basis = np.zeros(a.shape[1])
    basis[x.columns.get_loc(target)] = 1.0
    # beta_target = projection dot y; avoids refitting identical design at every scenario.
    return a @ np.linalg.solve(a.T @ a, basis)


def crossing(curve_delta: np.ndarray, curve_beta: np.ndarray, evaluator):
    # The negative departure nearest zero is the meaningful adverse tipping point.
    hits = []
    for i in range(len(curve_delta) - 1):
        if curve_beta[i] == 0:
            hits.append(float(curve_delta[i]))
        elif curve_beta[i] * curve_beta[i + 1] < 0:
            hits.append(float(optimize.brentq(evaluator, curve_delta[i], curve_delta[i + 1])))
    negative = [v for v in hits if v < 0]
    return max(negative) if negative else None


def run(source: Path, output: Path) -> None:
    plan = HERE / "ANALYSIS_PLAN.md"
    if not plan.exists():
        raise FileNotFoundError("Dated analysis plan must exist before sensitivity fits")
    supported, observed, info = prepare(source)
    missing = supported.loc[~supported.observed]
    obs_mask = supported.observed.to_numpy(bool)
    mis_mask = ~obs_mask
    z_missing = missing.neglect_z.to_numpy(float)
    centered_z = z_missing - z_missing.mean()
    info["unmatched_neglect_mean_z"] = float(z_missing.mean())
    info["unmatched_neglect_z_p10_p50_p90"] = [float(v) for v in np.quantile(z_missing, [.1, .5, .9])]

    xd_full, xd_obs = ranked_design(supported, observed, DISTRESS_FORMULA)
    xc_full, xc_obs = ranked_design(supported, observed, CONTENT_FORMULA)
    yd_obs = (observed.distress_future - info["full_baseline_scaling"]["distress"]["mean"]) / info["full_baseline_scaling"]["distress"]["sd"]
    yc_obs = observed.content_gap_future
    anchor_d = float(sm.OLS(yd_obs, xd_obs).fit().params["neglect_z"])
    anchor_c = float(sm.OLS(yc_obs, xc_obs).fit().params["neglect_z"])
    if abs(anchor_d - FROZEN_DISTRESS_BETA) > 1e-10:
        raise ValueError(f"Distress anchor mismatch: {anchor_d} versus {FROZEN_DISTRESS_BETA}")
    if abs(anchor_c - FROZEN_CONTENT_BETA_PER_SD) > 1e-10:
        raise ValueError(f"Content anchor mismatch: {anchor_c} versus {FROZEN_CONTENT_BETA_PER_SD}")
    factor = 7 / info["full_baseline_scaling"]["neglect"]["sd"]
    if abs(anchor_c * factor - FROZEN_CONTENT_DELTA_17_VS_10) > 1e-10:
        raise ValueError("Content 17-vs-10 scaling mismatch")

    # Baseline-covariate-dependent reference means; these are counterfactual predictions,
    # not recovered W3 observations. Each observed outcome remains exactly as recorded.
    fd = reference_glm(xd_obs, observed.distress_future, 42.0)
    fa = reference_glm(xc_obs, observed.anhedonic_future, 3.0)
    fu = reference_glm(xc_obs, observed.autonomic_future, 3.0)
    pd_mis = 42 * expit(xd_full.loc[missing.index].to_numpy(float) @ fd.params.to_numpy(float))
    pa_mis = 3 * expit(xc_full.loc[missing.index].to_numpy(float) @ fa.params.to_numpy(float))
    pu_mis = 3 * expit(xc_full.loc[missing.index].to_numpy(float) @ fu.params.to_numpy(float))
    if not (np.isfinite(pd_mis).all() and np.isfinite(pa_mis).all() and np.isfinite(pu_mis).all()):
        raise ValueError("Nonfinite counterfactual reference mean")
    if not (np.all((pd_mis >= 0) & (pd_mis <= 42)) and
            np.all((pa_mis >= 0) & (pa_mis <= 3)) and
            np.all((pu_mis >= 0) & (pu_mis <= 3))):
        raise ValueError("Reference mean exceeds DASS bounds")

    hd = coef_projection(xd_full)
    hc = coef_projection(xc_full)
    distress_mean = info["full_baseline_scaling"]["distress"]["mean"]
    distress_sd = info["full_baseline_scaling"]["distress"]["sd"]
    observed_d_contribution = float(hd[obs_mask] @ yd_obs.to_numpy(float))
    observed_c_contribution = float(hc[obs_mask] @ yc_obs.to_numpy(float))

    def evaluate_distress(delta: float):
        expected = np.clip(pd_mis + delta * centered_z, 0, 42)
        beta = observed_d_contribution + float(hd[mis_mask] @ ((expected - distress_mean) / distress_sd))
        return beta, float(expected.mean()), int(np.count_nonzero((expected <= 0) | (expected >= 42)))

    def evaluate_content(delta: float):
        a = np.clip(pa_mis + .5 * delta * centered_z, 0, 3)
        u = np.clip(pu_mis - .5 * delta * centered_z, 0, 3)
        beta = observed_c_contribution + float(hc[mis_mask] @ (a - u))
        return beta * factor, float((a - u).mean()), int(np.count_nonzero((a <= 0) | (a >= 3) | (u <= 0) | (u >= 3)))

    # Independent direct OLS refits verify the algebraic projection used for the grid.
    projection_errors = []
    for delta in [0., -4., 4.]:
        y_all = np.empty(len(supported))
        y_all[obs_mask] = yd_obs.to_numpy(float)
        y_all[mis_mask] = (
            np.clip(pd_mis + delta * centered_z, 0, 42) - distress_mean
        ) / distress_sd
        explicit = float(sm.OLS(y_all, xd_full).fit().params["neglect_z"])
        projection_errors.append(abs(explicit - evaluate_distress(delta)[0]))
    for delta in [0., -.25, .25]:
        y_all = np.empty(len(supported))
        y_all[obs_mask] = yc_obs.to_numpy(float)
        a = np.clip(pa_mis + .5 * delta * centered_z, 0, 3)
        u = np.clip(pu_mis - .5 * delta * centered_z, 0, 3)
        y_all[mis_mask] = a - u
        explicit = float(sm.OLS(y_all, xc_full).fit().params["neglect_z"]) * factor
        projection_errors.append(abs(explicit - evaluate_content(delta)[0]))
    max_projection_error = float(max(projection_errors))
    if max_projection_error > 1e-10:
        raise ValueError("Scenario coefficient projection failed independent OLS refit check")

    d_grid = np.arange(-42, 42.0001, .5)
    c_grid = np.arange(-3, 3.0001, .025)
    records = []
    for endpoint, grid, evaluator, unit in [
        ("general_distress", d_grid, evaluate_distress, "W3 distress SD per W2 neglect SD"),
        ("content_A_minus_U", c_grid, evaluate_content, "W3 item-mean 17-vs-10 neglect contrast"),
    ]:
        for delta in grid:
            beta, missing_mean, n_at_bound = evaluator(float(delta))
            records.append({
                "endpoint": endpoint, "delta_assumed_per_neglect_SD": float(delta),
                "delta_units": ("raw DASS distress points" if endpoint == "general_distress"
                                else "raw 0–3 content-gap item-mean units"),
                "scenario_point_estimate": beta,
                "estimate_units": unit,
                "hypothetical_missing_mean": missing_mean,
                "hypothetical_missing_n_at_any_bound": n_at_bound,
                "supported_baseline_n": len(supported),
                "observed_n": len(observed), "unmatched_n": len(missing),
            })
    curve = pd.DataFrame(records)
    baseline_d, baseline_d_missing_mean, _ = evaluate_distress(0)
    baseline_c, baseline_c_missing_mean, _ = evaluate_content(0)
    tip_d = crossing(d_grid, curve.loc[curve.endpoint.eq("general_distress"), "scenario_point_estimate"].to_numpy(float),
                     lambda v: evaluate_distress(v)[0])
    tip_c = crossing(c_grid, curve.loc[curve.endpoint.eq("content_A_minus_U"), "scenario_point_estimate"].to_numpy(float),
                     lambda v: evaluate_content(v)[0])
    tipping_rows = []
    for endpoint, tip, evaluator in [
        ("general_distress", tip_d, evaluate_distress),
        ("content_A_minus_U", tip_c, evaluate_content),
    ]:
        if tip is not None:
            estimate, hypothetical_mean, n_at_bound = evaluator(tip)
            tipping_rows.append({
                "endpoint": endpoint,
                "negative_delta_point_zero_tipping": tip,
                "scenario_point_estimate_at_tipping": estimate,
                "hypothetical_missing_mean_at_tipping": hypothetical_mean,
                "hypothetical_missing_n_at_any_bound_at_tipping": n_at_bound,
                "hypothetical_missing_fraction_at_any_bound_at_tipping": n_at_bound / len(missing),
            })

    selected = []
    for endpoint, values, evaluator, observed_anchor, reference_anchor in [
        ("general_distress", [-8., -4., 0., 4., 8.], evaluate_distress, anchor_d, baseline_d),
        ("content_A_minus_U", [-.5, -.25, 0., .25, .5], evaluate_content, anchor_c * factor, baseline_c),
    ]:
        for delta in values:
            beta, missing_mean, n_at_bound = evaluator(delta)
            selected.append({
                "endpoint": endpoint, "delta_assumed_per_neglect_SD": delta,
                "scenario_point_estimate": beta,
                "difference_from_complete_case_anchor": beta - observed_anchor,
                "difference_from_reference_delta_0": beta - reference_anchor,
                "hypothetical_missing_mean": missing_mean,
                "hypothetical_missing_n_at_any_bound": n_at_bound,
            })

    # Aggregate only: no participant IDs or row-level predicted or imputed values.
    by_grade = (
        supported.groupby("grade", observed=True)["observed"].agg(["count", "sum"])
        .rename(columns={"count": "supported_baseline_n", "sum": "matched_n"})
        .reset_index()
    )
    by_grade["unmatched_n"] = by_grade.supported_baseline_n - by_grade.matched_n
    by_grade["matched_fraction"] = by_grade.matched_n / by_grade.supported_baseline_n
    info.update({
        "analysis_status": "post_result_exploratory_not_preregistered",
        "target": "W2 continuing-grade school-by-grade cells with >=1 matched W3 outcome only",
        "what_delta_0_means": "bounded fractional-logit conditional-mean reference among supported-cell unmatched records",
        "source_boundary": "no equivalence assumed between W2-only site_09 and W3-only site_03",
        "complete_case_anchor": {
            "distress_beta_per_W2_neglect_SD": anchor_d,
            "content_gap_beta_per_W2_neglect_SD": anchor_c,
            "content_gap_delta_score17_minus10": anchor_c * factor,
            "matches_frozen_to_1e_10": True,
        },
        "delta_zero_reference": {
            "distress_beta_per_W2_neglect_SD": baseline_d,
            "content_gap_delta_score17_minus10": baseline_c,
            "hypothetical_missing_distress_mean": baseline_d_missing_mean,
            "hypothetical_missing_content_gap_mean": baseline_c_missing_mean,
        },
        "negative_delta_point_estimate_zero_tipping": {
            "distress_raw_points_per_neglect_SD": tip_d,
            "content_gap_0to3_units_per_neglect_SD": tip_c,
            "interpretation": "counterfactual point estimate only, not a significance/p-value tipping point",
        },
        "reference_glm": {
            "distress_converged": bool(fd.converged),
            "anhedonic_converged": bool(fa.converged),
            "autonomic_converged": bool(fu.converged),
            "distress_design_columns": xd_full.shape[1],
            "content_design_columns": xc_full.shape[1],
            "missing_reference_distress_mean": float(pd_mis.mean()),
            "missing_reference_anhedonic_mean": float(pa_mis.mean()),
            "missing_reference_autonomic_mean": float(pu_mis.mean()),
        },
        "max_projection_vs_direct_OLS_absolute_error": max_projection_error,
        "software": {
            "python": sys.version.split()[0], "numpy": np.__version__,
            "pandas": pd.__version__, "patsy": patsy.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "inferential_note": "No CI or p value is calculated from pseudo-completed outcomes; counterfactual means are deterministic and have unquantified model uncertainty.",
    })
    output.mkdir(parents=True, exist_ok=True)
    curve.to_csv(output / "scenario_curve.csv", index=False, float_format="%.12g")
    pd.DataFrame(selected).to_csv(output / "selected_scenarios.csv", index=False, float_format="%.12g")
    pd.DataFrame(tipping_rows).to_csv(output / "tipping_scenarios.csv", index=False, float_format="%.12g")
    by_grade.to_csv(output / "supported_grade_counts.csv", index=False)
    (output / "diagnostics.json").write_text(json.dumps(info, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "complete_case_anchor": info["complete_case_anchor"],
        "delta_zero_reference": info["delta_zero_reference"],
        "negative_delta_point_estimate_zero_tipping": info["negative_delta_point_estimate_zero_tipping"],
        "output": str(output),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="Private questionnaire_analysis.csv with frozen SHA-256; read in place, never exported")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.source, args.output)
