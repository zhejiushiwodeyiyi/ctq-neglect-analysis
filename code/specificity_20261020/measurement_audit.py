"""Fixed seven-item DASS content measurement audit for the CTQ specificity plan.

Only aggregate files are written.  The source participant table is read in memory;
its rows and identifiers are never copied to the output directory.  The ordinal
diagnostic is limited-information: pairwise four-category polychoric correlations
followed by a fixed positive-loading two-factor least-squares model.  It does not
provide a formal clustered CFA test or measurement-invariance test.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize_scalar, least_squares
from scipy.stats import multivariate_normal, norm

REPO = Path(__file__).resolve().parents[2]
EXPECTED_SHA = '5342bd167496989d0ebf083894d24547c50b210356715cdeceb49b76b5ad352b'
SOURCE: Path
OUT: Path
CONTENTS = {'anhedonic_content': (3, 10, 16), 'autonomic_content': (2, 4, 7, 19)}
ITEMS = (3, 10, 16, 2, 4, 7, 19)
COLS = tuple(f'dass{i}' for i in ITEMS)
PAIRS = tuple((i, j) for i in range(7) for j in range(i + 1, 7))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def prepare() -> dict[str, pd.DataFrame]:
    raw = pd.read_csv(SOURCE, low_memory=False)
    assert len(raw) == 22265
    assert not raw.duplicated(['participant_code', 'wave']).any()
    b = raw.loc[(raw.wave == 2) & raw.ctq_source_complete & raw.dass_source_complete].set_index('participant_code')
    f = raw.loc[(raw.wave == 3) & raw.dass_source_complete].set_index('participant_code')
    ids = b.index.intersection(f.index)
    b, f = b.loc[ids].copy(), f.loc[ids].copy()
    assert len(b) == len(f) == 4219 and b.school_code.nunique() == 10
    assert b.grade.isin([4, 5, 7, 8, 10, 11]).all()
    assert b[list(COLS)].notna().all().all() and f[list(COLS)].notna().all().all()
    domain_items = {'dep_score': (3, 5, 10, 13, 16, 17, 21),
                    'anxiety_score': (2, 4, 7, 9, 15, 19, 20),
                    'stress_score': (1, 6, 8, 11, 12, 14, 18)}
    for frame in (b, f):
        assert frame[[f'dass{i}' for i in range(1, 22)]].isin([0, 1, 2, 3]).all().all()
        for domain, items in domain_items.items():
            calc = 2 * frame[[f'dass{i}' for i in items]].sum(axis=1)
            assert np.array_equal(calc.to_numpy(), frame[domain].to_numpy())
    masks = {'All_paired': np.ones(len(b), dtype=bool),
             'Primary_grades_4_5': b.grade.isin([4, 5]).to_numpy(),
             'Secondary_grade_ge7': b.grade.ge(7).to_numpy()}
    frames = {}
    for wave, frame in ((2, b), (3, f)):
        for stage, mask in masks.items():
            frames[f'wave{wave}__{stage}'] = frame.iloc[np.flatnonzero(mask)]
    assert len(frames['wave3__Primary_grades_4_5']) == 1813
    assert len(frames['wave3__Secondary_grade_ge7']) == 2406
    return frames


def threshold_vector(x: np.ndarray) -> np.ndarray:
    counts = np.bincount(x.astype(int), minlength=4)
    if (counts == 0).any():
        raise ValueError('Empty item category; four-category threshold unavailable')
    return norm.ppf(np.cumsum(counts[:3]) / counts.sum())


def rectangle_probabilities(t1: np.ndarray, t2: np.ndarray, rho: float) -> np.ndarray:
    lo1 = np.r_[-np.inf, t1]
    hi1 = np.r_[t1, np.inf]
    lo2 = np.r_[-np.inf, t2]
    hi2 = np.r_[t2, np.inf]
    cov = np.array([[1., rho], [rho, 1.]])
    p = np.empty((4, 4))
    for a in range(4):
        for b in range(4):
            p[a, b] = multivariate_normal.cdf([hi1[a], hi2[b]], cov=cov,
                                                lower_limit=[lo1[a], lo2[b]],
                                                abseps=1e-10, releps=1e-10)
    return p


def polychoric_pair(x: np.ndarray, y: np.ndarray, t1: np.ndarray, t2: np.ndarray) -> tuple[float, int, float]:
    cells = np.zeros((4, 4), dtype=int)
    np.add.at(cells, (x.astype(int), y.astype(int)), 1)
    def objective(rho: float) -> float:
        p = rectangle_probabilities(t1, t2, rho)
        if not np.isfinite(p).all() or np.any(p[cells > 0] <= 0):
            return np.inf
        return -float(np.sum(cells[cells > 0] * np.log(p[cells > 0])))
    opt = minimize_scalar(objective, method='bounded', bounds=(-.97, .97),
                          options={'xatol': 1e-5, 'maxiter': 100})
    if not opt.success or not np.isfinite(opt.fun):
        raise RuntimeError(f'Polychoric MLE failed: {opt.message}')
    expected = rectangle_probabilities(t1, t2, float(opt.x)) * len(x)
    return float(opt.x), int((cells == 0).sum()), float(expected.min())


def fixed_two_factor(poly: np.ndarray) -> dict:
    """Unit-factor-variance, positive primary loadings, zero cross-loadings."""
    def model(par: np.ndarray) -> np.ndarray:
        loads, corr = par[:7], par[7]
        fac = np.array([0, 0, 0, 1, 1, 1, 1])
        pred = np.eye(7)
        for i, j in PAIRS:
            pred[i, j] = pred[j, i] = loads[i] * loads[j] * (1 if fac[i] == fac[j] else corr)
        return pred
    def resid(par: np.ndarray) -> np.ndarray:
        pred = model(par)
        return np.array([pred[i, j] - poly[i, j] for i, j in PAIRS])
    bounds = (np.r_[np.repeat(.001, 7), -.97], np.r_[np.repeat(.999, 7), .97])
    fit = least_squares(resid, np.r_[np.repeat(.6, 7), .6], bounds=bounds,
                        xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=5000)
    residual = resid(fit.x)
    return {'converged': bool(fit.success), 'message': fit.message,
            'loadings': fit.x[:7].tolist(), 'factor_correlation': float(fit.x[7]),
            'rms_offdiagonal_residual': float(np.sqrt(np.mean(residual ** 2))),
            'max_abs_offdiagonal_residual': float(np.max(np.abs(residual))),
            'parameter_at_bound': bool(np.any(fit.x <= bounds[0] + .01) or np.any(fit.x >= bounds[1] - .01)),
            'n_offdiagonal_correlations': len(PAIRS), 'n_free_parameters': 8}


def main() -> None:
    parser = argparse.ArgumentParser(description='Seven-item DASS measurement diagnostics')
    parser.add_argument('--data', type=Path, required=True, help='Authorized participant table')
    parser.add_argument('--reference-descriptives', type=Path,
                        help='Optional frozen measurement_grade_content_descriptives.csv')
    parser.add_argument('--output', type=Path, required=True,
                        help='New local output directory; contains small item counts')
    opts = parser.parse_args()
    global SOURCE, OUT
    SOURCE, OUT = opts.data, opts.output.resolve()
    if sha256_file(SOURCE) != EXPECTED_SHA:
        raise ValueError('Private input does not match the locked SHA-256')
    if OUT.is_relative_to(REPO):
        raise ValueError('Keep small-cell outputs outside the public code repository')
    OUT.mkdir(parents=True, exist_ok=False)
    frames = prepare()
    frequencies, content_rows, threshold_rows, loading_rows, factor_rows, pair_rows = [], [], [], [], [], []
    for label, frame in frames.items():
        wave, stage = label.split('__')
        values = frame[list(COLS)].to_numpy(dtype=float)
        n = len(values)
        assert np.isin(values, [0, 1, 2, 3]).all()
        thresholds = []
        for j, item in enumerate(ITEMS):
            counts = np.bincount(values[:, j].astype(int), minlength=4)
            for category, count in enumerate(counts):
                frequencies.append(dict(wave=int(wave[4:]), stage=stage, item=item,
                                        category=category, count=int(count), denominator=n,
                                        percent=100 * float(count) / n, missing_n=0))
            try:
                t = threshold_vector(values[:, j])
            except ValueError:
                t = np.repeat(np.nan, 3)
            thresholds.append(t)
            for cut, threshold in enumerate(t, start=1):
                threshold_rows.append(dict(wave=int(wave[4:]), stage=stage, item=item,
                                           threshold_at_or_below=cut - 1,
                                           threshold=float(threshold)))
        for content, items in CONTENTS.items():
            ix = [ITEMS.index(item) for item in items]
            arr = values[:, ix]
            score = arr.mean(axis=1)
            covariance = np.cov(arr, rowvar=False, ddof=1)
            k = len(ix)
            alpha = k / (k - 1) * (1 - np.trace(covariance) / covariance.sum())
            r = np.corrcoef(arr, rowvar=False)
            content_rows.append(dict(wave=int(wave[4:]), stage=stage, content=content,
                                     n=n, schools=frame.school_code.nunique(), items=k,
                                     mean=float(score.mean()), sd=float(score.std(ddof=1)),
                                     zero_n=int((score == 0).sum()),
                                     zero_percent=100 * float((score == 0).mean()),
                                     alpha_raw=float(alpha),
                                     mean_interitem_r=float(r[np.triu_indices(k, 1)].mean()),
                                     between_content_r=float(np.corrcoef(values[:, :3].mean(axis=1),
                                                                          values[:, 3:].mean(axis=1))[0, 1]),
                                     missing_any_item_n=int(np.isnan(arr).any(axis=1).sum())))
        if np.isnan(thresholds).any():
            factor_rows.append(dict(wave=int(wave[4:]), stage=stage, n=n,
                                    status='not_estimable_empty_category'))
            continue
        poly = np.eye(7)
        for i, j in PAIRS:
            rho, empty_cells, minimum_expected = polychoric_pair(
                values[:, i], values[:, j], thresholds[i], thresholds[j])
            poly[i, j] = poly[j, i] = rho
            pair_rows.append(dict(wave=int(wave[4:]), stage=stage, item_i=ITEMS[i],
                                  item_j=ITEMS[j], polychoric_r=rho,
                                  observed_zero_cells=empty_cells,
                                  minimum_expected_cell=minimum_expected))
        fac = fixed_two_factor(poly)
        min_eigenvalue = float(np.linalg.eigvalsh(poly).min())
        factor_rows.append(dict(wave=int(wave[4:]), stage=stage, n=n,
                                status='estimated' if fac['converged'] else 'failed',
                                factor_r=fac['factor_correlation'],
                                polychoric_min_eigenvalue=min_eigenvalue,
                                rms_offdiagonal_residual=fac['rms_offdiagonal_residual'],
                                max_abs_offdiagonal_residual=fac['max_abs_offdiagonal_residual'],
                                parameter_at_bound=fac['parameter_at_bound'],
                                pairs_with_observed_zero_cells=sum(row['observed_zero_cells'] > 0 for row in pair_rows if row['wave'] == int(wave[4:]) and row['stage'] == stage),
                                pairs_with_min_expected_below5=sum(row['minimum_expected_cell'] < 5 for row in pair_rows if row['wave'] == int(wave[4:]) and row['stage'] == stage)))
        for j, item in enumerate(ITEMS):
            loading_rows.append(dict(wave=int(wave[4:]), stage=stage, item=item,
                                     factor='anhedonic_content' if j < 3 else 'autonomic_content',
                                     loading=fac['loadings'][j]))
    content_frame = pd.DataFrame(content_rows)
    assert abs(float(content_frame.loc[(content_frame.wave == 3) &
                                        (content_frame.stage == 'All_paired') &
                                        (content_frame.content == 'anhedonic_content'),
                                        'zero_percent'].iloc[0]) - 69.542546) < 1e-5
    assert abs(float(content_frame.loc[(content_frame.wave == 3) &
                                        (content_frame.stage == 'All_paired') &
                                        (content_frame.content == 'autonomic_content'),
                                        'zero_percent'].iloc[0]) - 61.057123) < 1e-5
    if opts.reference_descriptives is not None:
        old = pd.read_csv(opts.reference_descriptives)
        for row in old.itertuples(index=False):
            new = content_frame[(content_frame.wave == row.visit) &
                                (content_frame.stage == row.stage_at_baseline) &
                                (content_frame.content == row.content)]
            assert len(new) == 1
            for key in ('mean', 'sd', 'zero_percent', 'alpha_raw', 'between_content_r'):
                assert abs(float(new[key].iloc[0]) - float(getattr(row, key))) < 1e-12
    pd.DataFrame(frequencies).to_csv(OUT / 'item_category_frequencies.csv', index=False)
    content_frame.to_csv(OUT / 'content_descriptives.csv', index=False)
    pd.DataFrame(threshold_rows).to_csv(OUT / 'ordinal_thresholds.csv', index=False)
    pd.DataFrame(loading_rows).to_csv(OUT / 'ordinal_two_factor_loadings.csv', index=False)
    pd.DataFrame(factor_rows).to_csv(OUT / 'ordinal_two_factor_diagnostics.csv', index=False)
    pd.DataFrame(pair_rows).to_csv(OUT / 'polychoric_pair_diagnostics.csv', index=False)
    result = {'source_sha256': sha256_file(SOURCE), 'code_sha256': sha256_file(Path(__file__)),
              'analysis_status': 'post_result_exploratory', 'paired_n': 4219, 'schools': 10,
              'item_sets': CONTENTS, 'software': {'python': platform.python_version(),
                                                   'numpy': np.__version__, 'pandas': pd.__version__,
                                                   'scipy': scipy.__version__},
              'method': 'pairwise four-category polychoric MLE; fixed two-factor unweighted least squares; descriptive only',
              'no_participant_level_outputs': True,
              'reference_descriptives_checked': opts.reference_descriptives is not None,
              'publication_caveat': 'Item-category frequencies may include small cells; do not publish without data-custodian review.'}
    (OUT / 'audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(pd.DataFrame(content_rows).to_string(index=False))
    print(pd.DataFrame(factor_rows).to_string(index=False))


if __name__ == '__main__':
    main()
