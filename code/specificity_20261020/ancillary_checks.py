"""Fixed ancillary checks for the 2026-10-20 CTQ content-specificity protocol.

Reads the private participant table in place. Writes only aggregate results.
This is post-result exploratory work; none of the fits is an independent test.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
import patsy
from patsy import dmatrix
import scipy
from scipy import linalg, stats
import sklearn
import statsmodels
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

threadpool_limits(limits=1)

REPO = Path(__file__).resolve().parents[2]
DEFAULT_IPW_CODE = REPO / "archive/原临床分析代码/revision_20261012/ipw_analysis.py"
OUT: Path
SOURCE: Path
OLD_IPW: Path
OLD_IPW_SUMMARY: Path
SECONDARY: Path
SUB = ["ctq_affective_ignore_score", "ctq_body_ignore_score", "ctq_affective_abuse_score", "ctq_body_abuse_score"]
DOMAINS = ["dep_score", "anxiety_score", "stress_score"]
EXPECTED_SHA = "5342bd167496989d0ebf083894d24547c50b210356715cdeceb49b76b5ad352b"


def describe_contents(m: pd.DataFrame) -> pd.DataFrame:
    m = m.copy()
    m["neglect"] = m[SUB[:2]].sum(axis=1, min_count=2)
    m["abuse"] = m[SUB[2:]].sum(axis=1, min_count=2)
    m["anhedonic"] = m[["dass3", "dass10", "dass16"]].sum(axis=1, min_count=3) / 3
    m["autonomic"] = m[["dass2", "dass4", "dass7", "dass19"]].sum(axis=1, min_count=4) / 4
    m["site"] = m.school_code.astype(str)
    m["sex"] = m.gender.astype("string").str.strip().str.lower().map(
        {"female": 1.0, "male": 0.0, "女": 1.0, "男": 0.0}
    )
    assert m.sex.notna().all()
    return m


def pair(m: pd.DataFrame, first: int, last: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float]:
    baseline = m.loc[m.wave.eq(first) & m.ctq_source_complete & m.dass_source_complete].set_index("participant_code")
    future = m.loc[m.wave.eq(last) & m.dass_source_complete].set_index("participant_code")
    ids = baseline.index.intersection(future.index)
    d = baseline.loc[ids].copy()
    d["anhedonic_future"] = future.loc[ids, "anhedonic"]
    d["autonomic_future"] = future.loc[ids, "autonomic"]
    neglect_sd = float(baseline.neglect.std(ddof=0))
    for c in ["neglect", "abuse"] + DOMAINS:
        sd = float(baseline[c].std(ddof=0))
        d[c + "_z"] = (d[c] - float(baseline[c].mean())) / sd
    assert d[["anhedonic_future", "autonomic_future", "neglect_z", "abuse_z"] + [c + "_z" for c in DOMAINS]].notna().all().all()
    return baseline, future, d, neglect_sd


def stacked_fit(d: pd.DataFrame, neglect_sd: float, weights: pd.Series | None, label: str) -> dict:
    frames = []
    for content, is_anhedonic in [("autonomic", 0), ("anhedonic", 1)]:
        f = d[["grade", "sex", "site", "neglect_z", "abuse_z"] + [c + "_z" for c in DOMAINS]].copy()
        f["outcome"] = d[content + "_future"]
        f["content"] = is_anhedonic
        f["weight"] = 1.0 if weights is None else weights.reindex(d.index).to_numpy(float)
        frames.append(f)
    stack = pd.concat(frames, ignore_index=True)
    assert stack.weight.notna().all() and stack.outcome.between(0, 3).all()
    formula = "content*(neglect_z+abuse_z+dep_score_z+anxiety_score_z+stress_score_z+C(grade)*sex+C(site))"
    x = dmatrix(formula, stack, return_type="dataframe")
    # QR only removes nuisance columns, if any; the planned contrast is protected.
    a = x.to_numpy(float)
    _, r, piv = linalg.qr(a, mode="economic", pivoting=True)
    tol = max(a.shape) * np.finfo(float).eps * abs(r).max()
    rank = int((abs(np.diag(r)) > tol).sum())
    target = "content:neglect_z"
    keep = sorted(piv[:rank])
    assert target in x.columns[keep], "Planned contrast aliased"
    x = x.iloc[:, keep]
    assert np.linalg.matrix_rank(x.to_numpy()) == x.shape[1]
    g = stack.site.nunique()
    fit = sm.WLS(stack.outcome, x, weights=stack.weight).fit(
        cov_type="cluster", cov_kwds={"groups": stack.site, "use_correction": True}, use_t=True
    )
    factor = 7 / neglect_sd
    beta_sd = float(fit.params[target])
    se_sd = float(fit.bse[target])
    beta = beta_sd * factor
    se = se_sd * factor
    crit = float(stats.t.ppf(0.975, g - 1))
    return {
        "check": label, "n": len(d), "stacked_rows": len(stack), "schools": g,
        "matrix_rank": x.shape[1], "baseline_neglect_sd_ddof0": neglect_sd,
        "neglect_contrast": "17 minus 10", "units": "DASS item mean 0 to 3",
        "delta_anhedonic_minus_autonomic": beta,
        "school_cr1_se": se,
        "school_cr1_t_low": beta - crit * se,
        "school_cr1_t_high": beta + crit * se,
        "school_t_df": g - 1,
        "school_t_p": float(2 * stats.t.sf(abs(beta / se), g - 1)),
        "weighted": weights is not None,
        "weight_repeated_for_two_rows": weights is not None,
    }


def separate_slope_identity(d: pd.DataFrame, neglect_sd: float, weights: pd.Series | None,
                            stacked_delta: float) -> float:
    """Cross-check fully interacted stacked point estimate with two separate fits."""
    x = dmatrix(
        "neglect_z+abuse_z+dep_score_z+anxiety_score_z+stress_score_z+C(grade)*sex+C(site)",
        d, return_type="dataframe",
    )
    a = x.to_numpy(float)
    _, r, piv = linalg.qr(a, mode="economic", pivoting=True)
    tol = max(a.shape) * np.finfo(float).eps * abs(r).max()
    rank = int((abs(np.diag(r)) > tol).sum())
    keep = sorted(piv[:rank])
    assert "neglect_z" in x.columns[keep]
    x = x.iloc[:, keep]
    ww = np.ones(len(d)) if weights is None else weights.reindex(d.index).to_numpy(float)
    slopes = {}
    for name in ["anhedonic", "autonomic"]:
        slopes[name] = float(sm.WLS(d[name + "_future"], x, weights=ww).fit().params["neglect_z"])
    direct = (slopes["anhedonic"] - slopes["autonomic"]) * 7 / neglect_sd
    discrepancy = abs(direct - stacked_delta)
    assert discrepancy < 1e-10, "Stacked and separate content point estimates disagree"
    return discrepancy


def old_ipw_module():
    spec = importlib.util.spec_from_file_location("ctq_existing_ipw", OLD_IPW)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def ipw_weights(m: pd.DataFrame, d: pd.DataFrame) -> tuple[pd.Series, dict, pd.DataFrame]:
    old = old_ipw_module()
    a = m[m.wave.eq(2)].set_index("participant_code").copy()
    assert len(a) == 8902
    for c in ["neglect", "abuse"] + DOMAINS:
        a[c + "_z"] = (a[c] - a[c].mean()) / a[c].std(ddof=0)
    b_ids = m.loc[m.wave.eq(3), "participant_code"]
    a["followup"] = a.index.isin(b_ids)
    eligible = a[a.grade.isin([4, 5, 7, 8, 10, 11])].dropna(subset=["sex"]).copy()
    retained = eligible.index[eligible.followup]
    assert len(eligible) == 5941 and len(retained) == len(d) == 4219
    assert set(retained) == set(d.index)
    XP = dmatrix(old.PS_FORMULA, eligible, return_type="dataframe")
    p, fit_info = old.fit_ps(XP, eligible.followup, revised=True)
    weights, info = old.weight_set(p, retained, float(eligible.followup.mean()))
    w = weights["stabilized_p99"].reindex(d.index)
    assert w.notna().all()
    previous = json.loads(OLD_IPW_SUMMARY.read_text())
    assert abs(info["p99_cap"] - previous["revised_weights"]["p99_cap"]) < 1e-5
    assert abs(info["stabilized_p99"]["ess"] - previous["revised_weights"]["stabilized_p99"]["ess"]) < 1e-4
    assert abs(w.sum() - previous["revised_weights"]["stabilized_p99"]["sum"]) < 1e-3
    # Existing support rule: discard grade-school cells with zero observed follow-up.
    cells = eligible.groupby(["grade", "site"], observed=True).followup.agg(["sum", "count"])
    zeros = set(cells[cells["sum"].eq(0)].index.tolist())
    supported = eligible.loc[[ (grade, site) not in zeros for grade, site in zip(eligible.grade, eligible.site) ]]
    assert len(supported) == 4519 and set(supported.index[supported.followup]) == set(d.index)
    bal = old.balance(supported, eligible.loc[retained], {"stabilized_p99": weights["stabilized_p99"]})
    vals = bal[bal.specification.eq("stabilized_p99")].smd.abs().dropna()
    weight_info = {
        "observation_model": old.PS_FORMULA,
        "propensity_fit_iterations": fit_info["n_iter"],
        "propensity_fit_warnings": fit_info["warnings"],
        "n_baseline_continuing": len(eligible), "n_supported_target": len(supported),
        "n_paired": len(d), "n_paired_schools": d.site.nunique(),
        "n_zero_observation_cells": len(zeros),
        "p99_cap": info["p99_cap"], "n_truncated": info["n_truncated"],
        "weight_min": info["stabilized_p99"]["min"],
        "weight_max": info["stabilized_p99"]["max"],
        "weight_ess": info["stabilized_p99"]["ess"],
        "supported_balance_max_abs_smd": float(vals.max()),
        "supported_balance_mean_abs_smd": float(vals.mean()),
        "supported_balance_n_indicators": len(vals),
        "supported_balance_n_abs_smd_over_0_1": int((vals > .1).sum()),
    }
    return w, weight_info, bal


def secondary_existing(neglect_sd: float) -> dict:
    rows = pd.read_csv(SECONDARY)
    hit = rows.loc[(rows.units == "raw_0to3") & (rows.exposure == "neglect")]
    assert len(hit) == 1
    r = hit.iloc[0]
    factor = 7 / neglect_sd
    return {
        "check": "grade_ge7_existing_common_prior_linear_OLS_rescaled",
        "n": int(r.n_participants), "schools": int(r.schools),
        "baseline_neglect_sd_ddof0": neglect_sd,
        "neglect_contrast": "17 minus 10", "units": "DASS item mean 0 to 3",
        "delta_anhedonic_minus_autonomic": float(r.beta * factor),
        "school_cr1_se": float(r.se_school * factor),
        "school_cr1_t_low": float(r.school_low * factor),
        "school_cr1_t_high": float(r.school_high * factor),
        "school_t_df": int(r.school_df),
        "school_t_p": float(r.p_school),
        "source_sha256": hashlib.sha256(SECONDARY.read_bytes()).hexdigest(), "new_fit": False,
    }


def main():
    parser = argparse.ArgumentParser(description="Post-result content-comparison ancillary checks")
    parser.add_argument("--data", type=Path, required=True, help="Authorized participant table")
    parser.add_argument("--ipw-code", type=Path, default=DEFAULT_IPW_CODE)
    parser.add_argument("--ipw-summary", type=Path, required=True, help="Frozen aggregate IPW summary JSON")
    parser.add_argument("--secondary", type=Path, required=True, help="Frozen grade-restricted common-prior table")
    parser.add_argument("--output", type=Path, required=True, help="New output directory outside this repository")
    opts = parser.parse_args()
    global SOURCE, OLD_IPW, OLD_IPW_SUMMARY, SECONDARY, OUT
    SOURCE, OLD_IPW, OLD_IPW_SUMMARY, SECONDARY = (
        opts.data, opts.ipw_code, opts.ipw_summary, opts.secondary)
    for source_file in (SOURCE, OLD_IPW, OLD_IPW_SUMMARY, SECONDARY):
        if not source_file.is_file():
            raise FileNotFoundError(source_file)
    OUT = opts.output.resolve()
    if OUT.is_relative_to(REPO):
        raise ValueError("Keep generated outputs outside the public code repository")
    OUT.mkdir(parents=True, exist_ok=False)
    sha = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert sha == EXPECTED_SHA, "Private source data changed; do not silently rerun."
    m = describe_contents(pd.read_csv(SOURCE, low_memory=False))
    assert len(m) == 22265 and not m.duplicated(["participant_code", "wave"]).any()
    b1, f2, d12, sd1 = pair(m, 1, 2)
    b2, f3, d23, sd2 = pair(m, 2, 3)
    assert len(d12) == 2022 and d12.site.nunique() == 4
    assert len(d23) == 4219 and d23.site.nunique() == 10
    assert abs(sd2 - 6.741979) < .001
    results = [stacked_fit(d12, sd1, None, "wave1_to_2_common_prior_linear_OLS"),
               secondary_existing(sd2)]
    w, weight_info, balance = ipw_weights(m, d23)
    results.append(stacked_fit(d23, sd2, w, "wave2_to_3_IPW_p99_common_prior_linear_WLS"))
    # The unweighted match is a required identity check against frozen Table S14.
    unweighted = stacked_fit(d23, sd2, None, "wave2_to_3_unweighted_identity_check")
    assert abs(unweighted["delta_anhedonic_minus_autonomic"] - .02864053969648072) < 1e-10
    assert abs(unweighted["school_cr1_t_low"] - .0038424127270906307) < 1e-10
    separate_identity_errors = {
        "wave1_to_2": separate_slope_identity(d12, sd1, None, results[0]["delta_anhedonic_minus_autonomic"]),
        "wave2_to_3_ipw": separate_slope_identity(d23, sd2, w, results[2]["delta_anhedonic_minus_autonomic"]),
        "wave2_to_3_unweighted": separate_slope_identity(d23, sd2, None, unweighted["delta_anhedonic_minus_autonomic"]),
    }
    pd.DataFrame(results + [unweighted]).to_csv(OUT / "ancillary_contrasts.csv", index=False)
    balance.to_csv(OUT / "ipw_supported_target_balance.csv", index=False)
    metadata = {
        "source_sha256": sha, "source_rows": len(m), "baseline_wave1_n": len(b1),
        "baseline_wave2_n": len(b2), "wave1_to_2_paired_n": len(d12),
        "wave2_to_3_paired_n": len(d23), "ipw": weight_info,
        "wave1_to_2_exposure_support": {
            "paired_min": float(d12.neglect.min()), "paired_max": float(d12.neglect.max()),
            "paired_n_score_10": int(d12.neglect.eq(10).sum()),
            "paired_n_score_17": int(d12.neglect.eq(17).sum()),
            "paired_q25": float(d12.neglect.quantile(.25)),
            "paired_q75": float(d12.neglect.quantile(.75)),
        },
        "unweighted_identity_check": unweighted,
        "stacked_vs_separate_slope_identity_abs_errors": separate_identity_errors,
        "interval_caveats": [
            "Four schools make wave1-to-2 CR1 t(3) interval particularly unstable.",
            "These are post-result exploratory checks, not independent validation.",
            "Wave1-to-2 and wave2-to-3 have different eligible schools and students.",
            "Fitted IPW probabilities are treated as fixed in the school-CR1 interval.",
            "IPW addresses measured follow-up selection only within observed support; one school and grade-school cells with no observed outcomes cannot be recovered.",
        ],
        "software": {"python": sys.version.split()[0], "numpy": np.__version__,
                     "pandas": pd.__version__, "patsy": patsy.__version__,
                     "scipy": scipy.__version__, "statsmodels": statsmodels.__version__,
                     "sklearn": sklearn.__version__},
    }
    (OUT / "ancillary_audit.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(pd.DataFrame(results + [unweighted]).to_string(index=False))


if __name__ == "__main__":
    main()
