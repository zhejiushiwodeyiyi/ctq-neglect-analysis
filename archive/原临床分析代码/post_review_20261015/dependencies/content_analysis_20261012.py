"""Fixed post-review content/subscale sensitivities; see CONTENT_PLAN.md.

No source files or historical model ledgers are modified. Outputs are aggregates.
"""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import sys
sys.dont_write_bytecode=True

import numpy as np
import pandas as pd
import patsy
import scipy
from scipy import stats
import statsmodels
import statsmodels.api as sm
from patsy import dmatrix
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / 'reproduction_output/revision_20261012_content'
SOURCE = ROOT / 'data/questionnaire_analysis.csv'
EVIDENCE = HERE / 'reference_results'
SCALING_REFERENCE = EVIDENCE / 'manuscript_scaling.json'
sys.path.insert(0, str(HERE.parent.parent))
from revision_common import SUB, bh, full_rank, sex_code

threadpool_limits(limits=1)
ITEMS = {'anhedonic_content': [3, 10, 16], 'autonomic_content': [2, 4, 7, 19]}
DOMAINS = ['dep_score', 'anxiety_score', 'stress_score']
LABELS = {'ctq_affective_ignore_score': 'Emotional neglect',
          'ctq_body_ignore_score': 'Physical neglect',
          'ctq_affective_abuse_score': 'Emotional abuse',
          'ctq_body_abuse_score': 'Physical abuse'}


def school_fit(x, y, school):
    x = full_rank(x)
    assert np.linalg.matrix_rank(x.to_numpy()) == x.shape[1]
    fit = sm.OLS(y, x).fit(cov_type='cluster',
                          cov_kwds={'groups': school, 'use_correction': True}, use_t=True)
    groups = school.nunique()
    beta, se = float(fit.params.ctq_z), float(fit.bse.ctq_z)
    critical = stats.t.ppf(.975, groups - 1)
    ci = fit.conf_int().loc['ctq_z']
    assert np.max(np.abs(ci - [beta - critical * se, beta + critical * se])) < 1e-10
    return x, dict(beta=beta, se_school=se, school_low=float(ci.iloc[0]),
                   school_high=float(ci.iloc[1]), p_school=float(fit.pvalues.ctq_z),
                   schools=int(groups), school_df=int(groups - 1),
                   matrix_rank=int(x.shape[1]), matrix_columns=int(x.shape[1]))


def wild_matrix(x, y, school):
    a = x.to_numpy(float)
    yy = np.asarray(y, float)
    target = x.columns.get_loc('ctq_z')
    inverse = np.linalg.pinv(a)
    influence = inverse[target]
    codes = pd.Categorical(school).codes
    groups = int(school.nunique())
    n, k = a.shape
    correction = groups / (groups - 1) * (n - 1) / (n - k)

    def tstat(v):
        coef = inverse @ v
        error = v - a @ coef
        cluster_score = np.bincount(codes, weights=influence * error, minlength=groups)
        se = np.sqrt(np.sum(cluster_score ** 2) * correction)
        return float(coef[target] / se)

    observed = tstat(yy)
    null = np.delete(a, target, axis=1)
    fitted = null @ np.linalg.lstsq(null, yy, rcond=None)[0]
    residual = yy - fitted
    exceed = 0
    for signs in itertools.product([-1., 1.], repeat=groups):
        bootstrap_t = tstat(fitted + residual * np.asarray(signs)[codes])
        exceed += abs(bootstrap_t) >= abs(observed) - 1e-10
    return dict(p_wild=float(exceed / 2 ** groups), wild_patterns=int(2 ** groups),
                observed_t=observed, bootstrap_exceedance_count=int(exceed))


def prepare():
    m = pd.read_csv(SOURCE, low_memory=False)
    assert len(m) == 22265
    assert not m.duplicated(['participant_code', 'wave']).any()
    m['sex'] = sex_code(m.gender)
    m['site'] = m.school_code
    for content, items in ITEMS.items():
        columns = [f'dass{i}' for i in items]
        m[content] = m[columns].sum(axis=1, min_count=len(items)) / len(items)
    m['neglect'] = m[SUB[:2]].sum(axis=1, min_count=2)
    m['abuse'] = m[SUB[2:]].sum(axis=1, min_count=2)
    m['distress'] = m[DOMAINS].mean(axis=1)
    baseline = m[m.wave.eq(2) & m.ctq_source_complete & m.dass_source_complete].set_index('participant_code')
    future = m[m.wave.eq(3) & m.dass_source_complete].set_index('participant_code')
    ids = baseline.index.intersection(future.index)
    d = baseline.loc[ids].copy()
    assert len(baseline) == 8902 and len(d) == 4219
    assert d.site.nunique() == 10
    assert d.grade.isin([4, 5, 7, 8, 10, 11]).all()
    columns = ['sex', 'grade', 'site', 'neglect', 'abuse', 'distress'] + DOMAINS + SUB + list(ITEMS)
    assert d[columns].notna().all().all()
    assert future.loc[ids, DOMAINS + list(ITEMS)].notna().all().all()
    assert baseline[[f'dass{i}' for i in range(1, 22)]].isin([0, 1, 2, 3]).all().all()
    assert future.loc[ids, [f'dass{i}' for i in range(1, 22)]].isin([0, 1, 2, 3]).all().all()
    domain_items = {'dep_score': [3, 5, 10, 13, 16, 17, 21],
                    'anxiety_score': [2, 4, 7, 9, 15, 19, 20],
                    'stress_score': [1, 6, 8, 11, 12, 14, 18]}
    for domain, items in domain_items.items():
        for frame in [baseline, future.loc[ids]]:
            recomputed = 2 * frame[[f'dass{i}' for i in items]].sum(axis=1, min_count=len(items))
            assert np.max(np.abs(recomputed - frame[domain])) == 0
    assert baseline[SUB].ge(5).all().all() and baseline[SUB].le(25).all().all()
    scales = {}
    for c in ['neglect', 'abuse', 'distress'] + DOMAINS + SUB + list(ITEMS):
        mu, sd = float(baseline[c].mean()), float(baseline[c].std(ddof=0))
        scales[c] = dict(mean=mu, sd=sd, baseline_n=int(baseline[c].notna().sum()))
        d[c + '_z'] = (d[c] - mu) / sd
        d[c + '_future'] = future.loc[ids, c]
        d[c + '_future_z'] = (future.loc[ids, c] - mu) / sd
    old_scales = json.loads((SCALING_REFERENCE).read_text())
    scale_error = max(abs(scales[c][field] - old_scales[c][field])
                      for c in old_scales for field in ['mean', 'sd'])
    assert scale_error < 1e-10
    return baseline, future, d, scales, scale_error


def content_models(d):
    contrasts, components, identity_checks = [], [], []
    for units in ['standardized', 'raw_0to3']:
        for exposure in ['neglect', 'abuse']:
            frames, slopes = [], {}
            other = 'abuse' if exposure == 'neglect' else 'neglect'
            for content in ITEMS:
                f = d[['grade', 'sex', 'site']].copy()
                f['ctq_z'], f['other'] = d[exposure + '_z'], d[other + '_z']
                for domain in DOMAINS:
                    f[domain + '_z'] = d[domain + '_z']
                f['outcome'] = d[content + ('_future_z' if units == 'standardized' else '_future')]
                formula = 'ctq_z+other+dep_score_z+anxiety_score_z+stress_score_z+C(grade)*sex+C(site)'
                x, row = school_fit(dmatrix(formula, f, return_type='dataframe'), f.outcome, f.site)
                slopes[content] = row['beta']
                components.append(dict(units=units, exposure=exposure, content=content,
                                       baseline_adjustment='common_three_DASS_domains',
                                       n_participants=len(f), **row))
                f = f.rename(columns={'ctq_z': 'exposure'})
                f['content'] = int(content == 'anhedonic_content')
                frames.append(f)
            f = pd.concat(frames, ignore_index=True)
            x = dmatrix('content*(exposure+other+dep_score_z+anxiety_score_z+stress_score_z+C(grade)*sex+C(site))',
                        f, return_type='dataframe').rename(columns={'content:exposure': 'ctq_z'})
            x, row = school_fit(x, f.outcome, f.site)
            bootstrap = wild_matrix(x, f.outcome, f.site)
            assert abs(bootstrap['observed_t'] - row['beta'] / row['se_school']) < 1e-8
            error = abs(row['beta'] - (slopes['anhedonic_content'] - slopes['autonomic_content']))
            assert error < 1e-10
            identity_checks.append(dict(units=units, exposure=exposure, largest_absolute_error=error))
            contrasts.append(dict(units=units, exposure=exposure, contrast='anhedonic_minus_autonomic',
                                  baseline_adjustment='common_three_DASS_domains', n_participants=len(d),
                                  stacked_rows=len(f), **row, **bootstrap))
    result = pd.DataFrame(contrasts)
    result['q_wild_2'] = result.groupby('units').p_wild.transform(bh)
    result.to_csv(OUT / 'content_revision_common_baseline_contrasts.csv', index=False)
    pd.DataFrame(components).to_csv(OUT / 'content_revision_common_baseline_slopes.csv', index=False)
    return result, identity_checks


def preserve_comparison(new):
    rows = []
    for units, filename, field in [('standardized', 'source_clinical_content_contrast.csv', 'beta_difference'),
                                    ('raw_0to3', 'content_raw_unit_contrast.csv', 'beta_difference_raw')]:
        old = pd.read_csv(EVIDENCE / filename)
        old = old[old.transition.eq('2-3')]
        for row in old.to_dict('records'):
            rows.append(dict(units=units, exposure=row['exposure'],
                             baseline_adjustment='different_own_prior_content',
                             beta=row[field], school_low=row['school_low'], school_high=row['school_high'],
                             p_wild=row['p_wild']))
    for row in new.to_dict('records'):
        rows.append({key: row[key] for key in ['units', 'exposure', 'baseline_adjustment',
                                              'beta', 'school_low', 'school_high', 'p_wild']})
    pd.DataFrame(rows).to_csv(OUT / 'content_revision_original_vs_common_baseline.csv', index=False)


def subscale_models(d, scales):
    rows = []
    for target in SUB:
        f = d[['grade', 'sex', 'site']].copy()
        f['ctq_z'] = d[target + '_z']
        others = [c for c in SUB if c != target]
        for c in others:
            f[c + '_z'] = d[c + '_z']
        f['prior'] = d.distress_z
        f['outcome'] = d.distress_future_z
        formula = 'ctq_z+prior+C(grade)*sex+C(site)+' + '+'.join(c + '_z' for c in others)
        x, row = school_fit(dmatrix(formula, f, return_type='dataframe'), f.outcome, f.site)
        bootstrap = wild_matrix(x, f.outcome, f.site)
        assert abs(bootstrap['observed_t'] - row['beta'] / row['se_school']) < 1e-8
        factor = scales['distress']['sd'] / scales[target]['sd'] * 5
        rows.append(dict(subscale=target, subscale_label=LABELS[target], n_participants=len(f),
                         baseline_subscale_sd=scales[target]['sd'], **row, **bootstrap,
                         raw_distress_per_5points=row['beta'] * factor,
                         raw_distress_per_5points_school_low=row['school_low'] * factor,
                         raw_distress_per_5points_school_high=row['school_high'] * factor))
    result = pd.DataFrame(rows)
    result['q_wild_4'] = bh(result.p_wild)
    result.to_csv(OUT / 'subscale_revision_mutually_adjusted_estimates.csv', index=False)
    return result


def content_describe(frame, stage, visit):
    rows = []
    between = float(frame[list(ITEMS)].corr().iloc[0, 1])
    for content, items in ITEMS.items():
        array = frame[[f'dass{i}' for i in items]].to_numpy(float)
        assert np.isfinite(array).all()
        covariance = np.cov(array, rowvar=False, ddof=1)
        k = len(items)
        alpha = k / (k - 1) * (1 - np.trace(covariance) / covariance.sum())
        correlations = np.corrcoef(array, rowvar=False)
        score = array.mean(axis=1)
        rows.append(dict(sample='paired_wave2_to_wave3', stage_at_baseline=stage, visit=visit,
                         content=content, n=len(frame), schools=frame.site.nunique(), items=k,
                         mean=score.mean(), sd=score.std(ddof=1), zero_n=int((score == 0).sum()),
                         zero_percent=100 * (score == 0).mean(), alpha_raw=alpha,
                         mean_interitem_r=correlations[np.triu_indices(k, 1)].mean(),
                         between_content_r=between))
    return rows


def grade_describe(baseline, future, d):
    measurement, ctq = [], []
    for stage, condition in [('Primary_grades_4_5', lambda f: f.grade.isin([4, 5])),
                              ('Secondary_grade_ge7', lambda f: f.grade.ge(7))]:
        paired = d[condition(d)]
        future_subset = future.loc[paired.index].copy()
        assert len(paired) == (1813 if stage == 'Primary_grades_4_5' else 2406)
        measurement.extend(content_describe(paired, stage, 2))
        measurement.extend(content_describe(future_subset, stage, 3))
        for sample, frame in [('full_eligible_wave2', baseline[condition(baseline)]),
                              ('paired_wave2_to_wave3_baseline', paired)]:
            for subscale in SUB:
                score = frame[subscale].dropna()
                ctq.append(dict(sample=sample, stage_at_baseline=stage, subscale=subscale,
                                subscale_label=LABELS[subscale], n=len(score), schools=frame.site.nunique(),
                                mean=score.mean(), sd=score.std(ddof=1), median=score.median(),
                                q25=score.quantile(.25), q75=score.quantile(.75), minimum=score.min(),
                                maximum=score.max(), floor_score5_n=int((score == 5).sum()),
                                floor_score5_percent=100 * (score == 5).mean()))
    pd.DataFrame(measurement).to_csv(OUT / 'measurement_grade_content_descriptives.csv', index=False)
    pd.DataFrame(ctq).to_csv(OUT / 'measurement_grade_CTQ_descriptives.csv', index=False)
    return pd.DataFrame(measurement), pd.DataFrame(ctq)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, future, d, scales, scaling_error = prepare()
    print('Prepared full wave-2 baseline n=8902; paired n=4219, schools=10.', flush=True)
    content, checks = content_models(d)
    old_components = pd.read_csv(EVIDENCE / 'clinical_sensitivity.csv')
    new_components = pd.read_csv(OUT / 'content_revision_common_baseline_slopes.csv')
    component_checks = []
    for row in new_components[new_components.units.eq('standardized')].to_dict('records'):
        old = old_components[(old_components.specification == 'prior_all_symptoms') &
                             (old_components.exposure == row['exposure']) &
                             (old_components.endpoint == row['content'])].iloc[0]
        error = max(abs(row['beta'] - old.beta), abs(row['school_low'] - old.school_low),
                    abs(row['school_high'] - old.school_high))
        assert error < 1e-10
        component_checks.append(dict(exposure=row['exposure'], content=row['content'],
                                     largest_estimate_or_interval_error=error))
    preserve_comparison(content)
    print('Four common-baseline direct contrasts completed.', flush=True)
    subscale = subscale_models(d, scales)
    print('Four mutually adjusted CTQ subscales completed.', flush=True)
    measurement, ctq = grade_describe(baseline, future, d)
    audit = dict(plan='CONTENT_PLAN.md', plan_sha256=hashlib.sha256((HERE / 'CONTENT_PLAN.md').read_bytes()).hexdigest(),
                 post_review_exploratory=True, source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                 source_row_count=22265, baseline_n=len(baseline), n_participants=len(d), schools=d.site.nunique(),
                 baseline_scaling=scales, scaling_ddof=0, original_scaling_max_error=scaling_error,
                 grade_strata_n={'primary_4_5': 1813, 'secondary_ge7': 2406},
                 direct_contrast_identity_checks=checks, identity_checks_passed=True,
                 existing_all_prior_component_checks=component_checks,
                 DASS_domain_item_reconstruction_passed=True, CTQ_baseline_subscale_range_checks_passed=True,
                 bootstrap='Restricted null, school Rademacher wild bootstrap, CR1 studentized, two-sided; full enumeration',
                 bootstrap_patterns=1024, school_CI='CR1 standard errors with t(9) critical value',
                 multiplicity={'content': 'BH across neglect and abuse, separately within each units representation',
                               'subscales': 'BH across the four mutually adjusted CTQ-subscale tests'},
                 versions={'python': sys.version.split()[0], 'numpy': np.__version__, 'pandas': pd.__version__,
                           'scipy': scipy.__version__, 'statsmodels': statsmodels.__version__, 'patsy': patsy.__version__},
                 limitations=['New post-review exploratory sensitivity analyses of the existing cohort.',
                              'Unit representations are not independent replications.',
                              'BH within these families does not correct all previous manuscript exploration.',
                              'Raw Cronbach alpha and item correlations do not establish construct or child-age validity.',
                              'Chronological age and actual administration materials require author confirmation.',
                              'School bootstrap is bootstrap inference, not an exact randomization test.'])
    (OUT / 'content_revision_audit.json').write_text(json.dumps(audit, indent=2))
    print(content.to_string(index=False))
    print(subscale.to_string(index=False))
    print(measurement.to_string(index=False))


if __name__ == '__main__':
    from revision_runtime import configure_output
    OUT = configure_output(ROOT, 'Re-run the fixed October 12 content/subscale analyses.')
    main()
