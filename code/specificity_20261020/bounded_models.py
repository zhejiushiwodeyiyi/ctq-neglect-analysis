"""Locked post-result CTQ content comparison: OLS anchors and bounded means.

Only aggregate estimates and diagnostics are written. The participant table remains
at its private source location. Run with --help for required paths.
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
from patsy import build_design_matrices, dmatrix
import scipy
from scipy import linalg, stats
import statsmodels
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

threadpool_limits(limits=1)

REPO = Path(__file__).resolve().parents[2]
DEFAULT_PRIOR_CODE = REPO / 'archive/原临床分析代码/post_review_20261015/post_review_analysis.py'
EXPECTED_SHA = '5342bd167496989d0ebf083894d24547c50b210356715cdeceb49b76b5ad352b'
EXPECTED_PLAN_SHA = '645da2b59fdab730fa604924cbf908b80066e5dfeda0440e2e3d733e0f208cda'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_prior(prior_code: Path, source: Path, scaling: Path):
    spec = importlib.util.spec_from_file_location('ctq_prior_for_bounded', prior_code)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.prior.SOURCE = source
    module.prior.SCALING_REFERENCE = scaling
    return module


def stack_content(d: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for content, flag in [('anhedonic_content', 1), ('autonomic_content', 0)]:
        part = d.copy()
        part['outcome_raw'] = part[content + '_future']
        part['outcome_fraction'] = part['outcome_raw'] / 3.0
        part['content'] = flag
        frames.append(part)
    return pd.concat(frames, ignore_index=True)


def rank_keep(x: pd.DataFrame):
    a = x.to_numpy(float)
    _, rr, pivot = linalg.qr(a, mode='economic', pivoting=True)
    tolerance = max(a.shape) * np.finfo(float).eps * np.abs(rr).max()
    rank = int((np.abs(np.diag(rr)) > tolerance).sum())
    kept = sorted(pivot[:rank])
    protected = [i for i, name in enumerate(x.columns) if 'neglect_z' in name]
    if not set(protected).issubset(kept):
        raise ValueError('A neglect exposure basis is aliased after rank reduction')
    return kept, rank


def scenario_matrices(info, kept, d: pd.DataFrame, mu: float, sd: float):
    out = {}
    for content, flag in [('A', 1), ('U', 0)]:
        for score in [10.0, 17.0]:
            frame = d.copy()
            frame['content'] = flag
            frame['neglect_z'] = (score - mu) / sd
            out[(content, score)] = np.asarray(build_design_matrices([info], frame)[0])[:, kept]
    return out


def marginal(fit, scenarios, bounded: bool):
    beta = fit.params.to_numpy()
    predictions, gradients = {}, {}
    for key, x in scenarios.items():
        linear = x @ beta
        if bounded:
            mu = scipy.special.expit(linear)
            predictions[key] = 3.0 * mu
            gradients[key] = 3.0 * (x * (mu * (1.0 - mu))[:, None]).mean(axis=0)
        else:
            predictions[key] = linear
            gradients[key] = x.mean(axis=0)
    a = float(predictions[('A', 17.0)].mean() - predictions[('A', 10.0)].mean())
    u = float(predictions[('U', 17.0)].mean() - predictions[('U', 10.0)].mean())
    gradient = (gradients[('A', 17.0)] - gradients[('A', 10.0)] -
                gradients[('U', 17.0)] + gradients[('U', 10.0)])
    variance = float(gradient @ fit.cov_params().to_numpy() @ gradient)
    if variance < -1e-12:
        raise ValueError('Negative delta-method variance')
    standard_error = float(np.sqrt(max(variance, 0.0)))
    return dict(delta=a-u, a_difference=a, u_difference=u,
                se_school=standard_error,
                predicted_at_10_A=float(predictions[('A', 10.0)].mean()),
                predicted_at_17_A=float(predictions[('A', 17.0)].mean()),
                predicted_at_10_U=float(predictions[('U', 10.0)].mean()),
                predicted_at_17_U=float(predictions[('U', 17.0)].mean()))


def finish_row(row, school_count: int):
    critical = float(stats.t.ppf(.975, school_count - 1))
    row['school_low'] = row['delta'] - critical * row['se_school']
    row['school_high'] = row['delta'] + critical * row['se_school']
    row['schools'] = school_count
    row['school_df'] = school_count - 1
    row['p_t_approx'] = float(2 * stats.t.sf(abs(row['delta'] / row['se_school']), school_count - 1))
    return row


def fit_fractional(x, f, start_params=None):
    model = sm.GLM(f['outcome_fraction'].to_numpy(float), x,
                   family=sm.families.Binomial(link=sm.families.links.Logit()))
    fit = model.fit(maxiter=250, tol=1e-9, start_params=start_params,
                    cov_type='cluster', cov_kwds={'groups': f['site'], 'use_correction': True},
                    use_t=True)
    if not fit.converged or not np.isfinite(fit.params).all():
        raise RuntimeError('Fractional-logit fit failed to converge with finite coefficients')
    return fit


def run():
    args = argparse.ArgumentParser(description='Reproduce the locked, post-result content comparison')
    args.add_argument('--data', type=Path, required=True, help='Authorized questionnaire_analysis.csv')
    args.add_argument('--prior-code', type=Path, default=DEFAULT_PRIOR_CODE)
    args.add_argument('--knots', type=Path, required=True, help='Frozen spline_knots.json')
    args.add_argument('--scaling', type=Path, required=True, help='Frozen manuscript_scaling.json')
    args.add_argument('--frozen', type=Path, required=True, help='Frozen Table S14 source table')
    args.add_argument('--original-plan', type=Path, help='Optional byte-identical original plan for hash check')
    args.add_argument('--output', type=Path, required=True, help='New output directory outside this repository')
    opts = args.parse_args()
    source, prior_code, knots_path, scaling, frozen_path = (
        opts.data, opts.prior_code, opts.knots, opts.scaling, opts.frozen)
    if sha(source) != EXPECTED_SHA:
        raise ValueError('Private input does not match locked SHA-256')
    if opts.original_plan is not None and sha(opts.original_plan) != EXPECTED_PLAN_SHA:
        raise ValueError('Locked analysis protocol has changed')
    for path in [prior_code, knots_path, scaling, frozen_path]:
        if not path.is_file():
            raise FileNotFoundError(path)
    output = opts.output.resolve()
    if output.is_relative_to(REPO):
        raise ValueError('Keep generated outputs outside the public code repository')
    output.mkdir(parents=True, exist_ok=False)
    prior = load_prior(prior_code, source, scaling)
    baseline, _, d, scales, scale_error = prior.prior.prepare()
    if len(d) != 4219 or d.site.nunique() != 10 or len(baseline) != 8902:
        raise ValueError('Paired sample does not match locked protocol')
    f = stack_content(d)
    if not f.outcome_fraction.between(0, 1).all():
        raise ValueError('Outcome fraction outside [0,1]')
    knots = json.loads(knots_path.read_text())
    frozen = pd.read_csv(frozen_path)
    result_rows, diagnostics, loso_rows = [], [], []
    for form in ['linear', 'natural_cubic_spline_3df']:
        formula = prior.formula_for(form, 'content', knots)
        raw_x = dmatrix(formula, f, return_type='dataframe')
        info = raw_x.design_info
        kept, rank = rank_keep(raw_x)
        x = raw_x.iloc[:, kept]
        scenarios = scenario_matrices(info, kept, d,
                                      scales['neglect']['mean'], scales['neglect']['sd'])
        ols = sm.OLS(f.outcome_raw.to_numpy(float), x).fit(
            cov_type='cluster', cov_kwds={'groups': f.site, 'use_correction': True}, use_t=True)
        ols_row = finish_row(marginal(ols, scenarios, bounded=False), 10)
        old = frozen[(frozen.system.eq('content')) &
                     (frozen.outcome_units.eq('raw_0to3')) &
                     (frozen.function_form.eq(form)) &
                     (frozen.exposure.eq('neglect'))].iloc[0]
        old_errors = {name: abs(ols_row[dest] - float(old[source]))
                      for name, dest, source in [('estimate', 'delta', 'estimate'),
                                                 ('se', 'se_school', 'se_school'),
                                                 ('lower', 'school_low', 'school_low'),
                                                 ('upper', 'school_high', 'school_high')]}
        if max(old_errors.values()) > 1e-9:
            raise ValueError(f'Frozen OLS reproduction failed for {form}: {old_errors}')
        result_rows.append(dict(model='OLS', form=form, n_participants=len(d),
                                n_rows=len(f), parameters=x.shape[1],
                                **ols_row))
        bounded = fit_fractional(x, f)
        bounded_row = finish_row(marginal(bounded, scenarios, bounded=True), 10)
        result_rows.append(dict(model='fractional_logit', form=form,
                                n_participants=len(d), n_rows=len(f), parameters=x.shape[1],
                                **bounded_row))
        fitted_raw = 3.0 * bounded.predict(x)
        condition_number = float(np.linalg.cond(x.to_numpy(float)))
        diagnostics.append(dict(form=form, formula=formula, design_columns_raw=raw_x.shape[1],
                                parameters=x.shape[1], design_rank=rank,
                                condition_number=condition_number,
                                converged=bool(bounded.converged),
                                iterations=int(bounded.fit_history['iteration']),
                                cov_type=bounded.cov_type,
                                cov_correction='statsmodels cluster use_correction=True; t_9 interval',
                                min_prediction=float(fitted_raw.min()),
                                max_prediction=float(fitted_raw.max()),
                                below_zero_n=int((fitted_raw < 0).sum()),
                                above_three_n=int((fitted_raw > 3).sum()),
                                near_zero_n=int((fitted_raw < 1e-6).sum()),
                                near_three_n=int((fitted_raw > 3-1e-6).sum()),
                                ols_reproduction_abs_errors=old_errors))
        for site in sorted(d.site.unique()):
            keep_person = d.site.ne(site)
            small_d = d.loc[keep_person].copy()
            keep_rows = f.site.ne(site)
            small_f = f.loc[keep_rows].copy()
            small_raw = raw_x.loc[keep_rows]
            small_kept, small_rank = rank_keep(small_raw)
            small_x = small_raw.iloc[:, small_kept]
            small_scenarios = scenario_matrices(info, small_kept, small_d,
                                                scales['neglect']['mean'], scales['neglect']['sd'])
            small_fit = fit_fractional(small_x, small_f)
            point = marginal(small_fit, small_scenarios, bounded=True)
            loso_rows.append(dict(form=form, removed_school=site,
                                  removed_participants=int((~keep_person).sum()),
                                  remaining_participants=len(small_d),
                                  remaining_schools=int(small_d.site.nunique()),
                                  n_rows=len(small_f), parameters=small_x.shape[1],
                                  design_rank=small_rank,
                                  converged=bool(small_fit.converged),
                                  delta=point['delta'],
                                  a_difference=point['a_difference'],
                                  u_difference=point['u_difference']))
        print(f'{form}: OLS {ols_row["delta"]:.6f}, fractional {bounded_row["delta"]:.6f}', flush=True)
    pd.DataFrame(result_rows).to_csv(output/'four_model_grid.csv', index=False)
    pd.DataFrame(loso_rows).to_csv(output/'fractional_leave_one_school_out.csv', index=False)
    provenance = dict(status='post-result exploratory, not prespecified before seeing cohort results',
                      original_protocol_sha256=EXPECTED_PLAN_SHA,
                      original_protocol_verified=opts.original_plan is not None,
                      source_sha256=sha(source), scaling_sha256=sha(scaling),
                      knots_sha256=sha(knots_path), frozen_table_sha256=sha(frozen_path),
                      script_sha256=sha(Path(__file__)),
                      upstream_code_sha256=sha(prior_code),
                      n_baseline=len(baseline), n_paired=len(d), schools=d.site.nunique(),
                      scale_max_abs_error=scale_error,
                      outcome='seven fixed DASS items, mean per item 0–3; fractional response = mean/3',
                      estimand='anhedonic minus autonomic marginal change from neglect 10 to 17',
                      reference_population='same observed paired covariate distribution for both contents',
                      libraries=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                                     patsy=patsy.__version__, statsmodels=statsmodels.__version__),
                      diagnostics=diagnostics,
                      note='Leave-one-school-out refits use full-sample scaling, fixed knots, and fixed Patsy design info; rank reduction removes only newly aliased columns.')
    (output/'diagnostics.json').write_text(json.dumps(provenance, indent=2, ensure_ascii=False)+'\n')


if __name__ == '__main__':
    run()
