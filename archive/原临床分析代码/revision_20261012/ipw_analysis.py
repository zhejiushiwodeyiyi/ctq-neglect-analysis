"""Post-review, aggregate-only observation-weight diagnostics and sensitivities.

See IPW_PLAN.md. Reads the deidentified working copy, never writes source files
or individual participant-level predictions/weights.
"""
from pathlib import Path
import hashlib
import json
import platform
import sys
sys.dont_write_bytecode=True
import warnings

import numpy as np
import pandas as pd
import patsy
from patsy import dmatrix
import scipy
from scipy import linalg, optimize, special, stats
import sklearn
from sklearn.linear_model import LogisticRegression
import statsmodels
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

threadpool_limits(limits=1)
HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'reproduction_output/revision_20261012_ipw'
SOURCE = ROOT / 'data/questionnaire_analysis.csv'
EVIDENCE = HERE / 'reference_results'
SUB = ['ctq_affective_ignore_score', 'ctq_body_ignore_score',
       'ctq_affective_abuse_score', 'ctq_body_abuse_score']
DOMAINS = {'dep_score': [3, 5, 10, 13, 16, 17, 21],
           'anxiety_score': [2, 4, 7, 9, 15, 19, 20],
           'stress_score': [1, 6, 8, 11, 12, 14, 18]}
ENDPOINTS = ['distress', 'dep_score', 'anxiety_score', 'stress_score',
             'anhedonic_content', 'autonomic_content', 'negative_self_content']
PS_FORMULA = ('C(grade)*sex + C(site) + abuse_z + neglect_z + '
              'dep_score_z + anxiety_score_z + stress_score_z')


def clean_json(x):
    if isinstance(x, dict):
        return {str(k): clean_json(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean_json(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def describe(v):
    x = np.asarray(v, dtype=float)
    result = {'n': len(x), 'mean': x.mean(), 'sd_ddof0': x.std(ddof=0),
              'min': x.min(), 'max': x.max()}
    result.update({f'p{100*q:g}': np.quantile(x, q)
                   for q in [.01, .05, .25, .5, .75, .95, .99]})
    return clean_json(result)


def full_rank(X):
    others = [c for c in X.columns if c != 'ctq_z']
    a = X[others].to_numpy(float)
    _, r, piv = linalg.qr(a, mode='economic', pivoting=True)
    tol = max(a.shape) * np.finfo(float).eps * abs(r).max()
    rank = int((abs(np.diag(r)) > tol).sum())
    keep = [others[j] for j in sorted(piv[:rank])]
    resid = X.ctq_z.to_numpy() - a[:, piv[:rank]] @ np.linalg.lstsq(
        a[:, piv[:rank]], X.ctq_z, rcond=None)[0]
    assert np.linalg.norm(resid) >= 1e-8
    return X[keep + ['ctq_z']]


def fit_ps(X, y, revised=False):
    kwargs = {'C': 1, 'max_iter': 2000}
    if revised:
        kwargs.update(solver='lbfgs', fit_intercept=True, max_iter=10000,
                      tol=1e-10, l1_ratio=0.0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        model = LogisticRegression(**kwargs).fit(X, y)
    p = pd.Series(model.predict_proba(X)[:, 1], index=X.index)
    pars = model.get_params()
    theta = np.r_[model.coef_[0], model.intercept_[0]]
    a = np.c_[X.to_numpy(float), np.ones(len(X))]
    yy = y.to_numpy(float)

    def loss_grad(v):
        eta = a @ v
        penalty = .5 * (v[:-1] @ v[:-1]) / pars['C']
        loss = np.logaddexp(0., eta).sum() - yy @ eta + penalty
        grad = a.T @ (special.expit(eta) - yy)
        grad[:-1] += v[:-1] / pars['C']
        return loss, grad

    loss, grad = loss_grad(theta)
    info = {'parameters': pars, 'n_iter': model.n_iter_.tolist(),
            'warnings': [{'category': w.category.__name__, 'message': str(w.message)}
                         for w in caught],
            'convergence_warning': any(w.category.__name__ == 'ConvergenceWarning'
                                       for w in caught),
            'max_iter_reached': bool(np.max(model.n_iter_) >= pars['max_iter']),
            'penalized_objective_sum_scale': loss,
            'gradient_abs_max_sum_scale': abs(grad).max(),
            'gradient_abs_max_mean_scale': abs(grad).max()/len(X),
            'patsy_intercept_coefficient': float(model.coef_[0][X.columns.get_loc('Intercept')]),
            'sklearn_intercept': float(model.intercept_[0])}
    if revised:
        cross = optimize.minimize(loss_grad, np.zeros(a.shape[1]), jac=True,
                                  method='L-BFGS-B', options={
                                      'maxiter': 10000, 'ftol': 1e-14,
                                      'gtol': 1e-10, 'maxls': 100})
        pcross = special.expit(a @ cross.x)
        info['independent_optimization'] = {
            'method': 'SciPy L-BFGS-B analytic penalized sum-log-likelihood; zero start',
            'success': cross.success, 'status': cross.status,
            'message': str(cross.message), 'n_iter': cross.nit,
            'gradient_abs_max_sum_scale': abs(cross.jac).max(),
            'objective_difference': float(loss - cross.fun),
            'prediction_max_abs_difference': float(abs(p.to_numpy()-pcross).max())}
        assert not info['convergence_warning'] and not info['max_iter_reached']
        assert cross.success
        assert abs(p.to_numpy()-pcross).max() < 1e-5
    return p, clean_json(info)


def weight_set(p, retained, retention):
    # Series alignment is explicit; random identifiers remain internal only.
    pp = p.reindex(retained)
    assert pp.notna().all() and pp.between(0, 1, inclusive='neither').all()
    unstable = 1/pp
    stable = retention/pp
    cap = float(stable.quantile(.99))
    capped = stable.clip(upper=cap)
    info = {'unstabilized': describe(unstable), 'stabilized_raw': describe(stable),
            'stabilized_p99': describe(capped), 'p99_cap': cap,
            'n_truncated': int((stable > cap).sum()),
            'fraction_truncated': float((stable > cap).mean())}
    for name, w in [('unstabilized', unstable), ('stabilized_raw', stable),
                    ('stabilized_p99', capped)]:
        s1 = float(np.sum(w.to_numpy()))
        s2 = float(np.dot(w.to_numpy(), w.to_numpy()))
        # Independent loop avoids accidentally calculating a Series with wrong index.
        s1_loop = sum(float(z) for z in w)
        s2_loop = sum(float(z)**2 for z in w)
        assert np.isclose(s1, s1_loop, rtol=1e-12)
        assert np.isclose(s2, s2_loop, rtol=1e-12)
        info[name].update(sum=s1, ess=s1*s1/s2)
    return {'unstabilized': unstable, 'stabilized_raw': stable,
            'stabilized_p99': capped}, clean_json(info)


def balance(eligible, paired, weights):
    features = {'sex_female': (eligible.sex, paired.sex, 'binary')}
    for c in ['grade', 'neglect', 'abuse', 'dep_score', 'anxiety_score', 'stress_score']:
        features[c] = (eligible[c], paired[c], 'continuous')
    for c in ['grade', 'site']:
        for level in sorted(eligible[c].unique()):
            features[f'{c}_{level}'] = ((eligible[c] == level).astype(float),
                                        (paired[c] == level).astype(float), 'binary')
    for grade, site in eligible[['grade','site']].drop_duplicates().sort_values(['grade','site']).itertuples(index=False):
        features[f'grade_{grade}_site_{site}'] = (
            ((eligible.grade == grade) & (eligible.site == site)).astype(float),
            ((paired.grade == grade) & (paired.site == site)).astype(float), 'binary')
    features['high_school_grades_10_11'] = (eligible.grade.isin([10,11]).astype(float),
                                            paired.grade.isin([10,11]).astype(float), 'binary')
    rows = []
    for name, (target, observed, kind) in features.items():
        target_mean = float(target.mean())
        denominator = (np.sqrt(target_mean * (1-target_mean)) if kind == 'binary'
                       else float(target.std(ddof=0)))
        for specification, w in {'unweighted': pd.Series(1., index=paired.index), **weights}.items():
            ww = w.reindex(paired.index)
            assert ww.notna().all()
            observed_mean = float(np.average(observed.to_numpy(), weights=ww.to_numpy()))
            rows.append({'variable': name, 'type': kind, 'specification': specification,
                         'target_mean': target_mean, 'observed_mean': observed_mean,
                         'target_sd_denominator': denominator,
                         'smd': (observed_mean-target_mean)/denominator if denominator > 0 else np.nan})
    return pd.DataFrame(rows)


def regress(paired, weights):
    rows = []
    manual_cov_differences = []
    for endpoint in ENDPOINTS:
        for exposure in ['neglect', 'abuse']:
            f = paired.copy()
            f['outcome'] = f[endpoint+'_future_z']
            f['ctq_z'] = f[exposure+'_z']
            f['other'] = f[('neglect' if exposure == 'abuse' else 'abuse')+'_z']
            f['prior'] = f[endpoint+'_z']
            for specification in ['all_prior_unweighted', 'all_prior_raw_weighted',
                                  'all_prior_p99_weighted'] + (
                                      ['own_prior_unweighted'] if endpoint == 'distress' else []):
                nuisance = ('prior' if specification == 'own_prior_unweighted'
                            else 'dep_score_z + anxiety_score_z + stress_score_z')
                X = full_rank(dmatrix('ctq_z + other + C(grade)*sex + C(site) + '+nuisance,
                                      f, return_type='dataframe'))
                ff = f.loc[X.index]
                if specification == 'all_prior_raw_weighted':
                    w = weights['stabilized_raw'].reindex(ff.index)
                elif specification == 'all_prior_p99_weighted':
                    w = weights['stabilized_p99'].reindex(ff.index)
                else:
                    w = pd.Series(1., index=ff.index)
                assert len(ff) == len(paired) and w.notna().all()
                fit = sm.WLS(ff.outcome, X, weights=w).fit()
                robust = fit.get_robustcov_results(cov_type='cluster', groups=ff.site,
                                                   use_correction=True, use_t=True)
                j = X.columns.get_loc('ctq_z')
                n, k = X.shape
                G = ff.site.nunique()
                beta, se = float(fit.params.ctq_z), float(robust.bse[j])
                crit = float(stats.t.ppf(.975, G-1))
                # Direct weighted CR1 sandwich, independent of statsmodels robust routine.
                xa = X.to_numpy(float); wa = w.to_numpy(float)
                bread = np.linalg.inv(xa.T @ (wa[:,None]*xa))
                resid = ff.outcome.to_numpy() - xa @ fit.params.to_numpy()
                scores = np.vstack([xa[ff.site.to_numpy()==g].T @
                                    (wa[ff.site.to_numpy()==g]*resid[ff.site.to_numpy()==g])
                                    for g in sorted(ff.site.unique())])
                cov = bread @ scores.T @ scores @ bread * G/(G-1)*(n-1)/(n-k)
                d = float(abs(np.sqrt(cov[j,j])-se))
                manual_cov_differences.append(d)
                assert d < 1e-10
                if specification == 'all_prior_raw_weighted':
                    unstable_fit = sm.WLS(ff.outcome, X, weights=weights['unstabilized'].reindex(ff.index)).fit()
                    assert np.max(abs(unstable_fit.params.to_numpy()-fit.params.to_numpy())) < 1e-10
                rows.append({'endpoint': endpoint, 'exposure': exposure,
                             'specification': specification, 'n': n, 'schools': G,
                             'rank': k, 'beta': beta, 'school_se': se,
                             'school_low': beta-crit*se, 'school_high': beta+crit*se,
                             'p_school': float(2*stats.t.sf(abs(beta/se), G-1)),
                             'interval_df': G-1, 'interval_method': 'school CR1 t interval; fitted weights treated as fixed'})
    return pd.DataFrame(rows), max(manual_cov_differences)


def main():
    file_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    m = pd.read_csv(SOURCE, low_memory=False)
    assert not m.duplicated(['participant_code','wave']).any()
    mismatch = {}
    for c, ii in DOMAINS.items():
        mismatch[c] = int((2*m[[f'dass{i}' for i in ii]].sum(axis=1) != m[c]).sum())
    assert not any(mismatch.values())
    m['neglect'] = m[SUB[:2]].sum(axis=1,min_count=2)
    m['abuse'] = m[SUB[2:]].sum(axis=1,min_count=2)
    m['distress'] = m[list(DOMAINS)].mean(axis=1)
    m['anhedonic_content'] = m[['dass3','dass10','dass16']].mean(axis=1)
    m['autonomic_content'] = m[['dass2','dass4','dass7','dass19']].mean(axis=1)
    m['negative_self_content'] = m[['dass17','dass21']].mean(axis=1)
    m['site'] = m.school_code.astype(str)
    m['sex'] = m.gender.astype('string').str.strip().str.lower().map(
        {'female':1.,'male':0.,'女':1.,'男':0.})
    a = m[m.wave==2].set_index('participant_code').copy()
    b = m[m.wave==3].set_index('participant_code').copy()
    a['followup'] = a.index.isin(b.index)
    scales = {}
    for c in ['neglect','abuse']+ENDPOINTS:
        mean = float(a[c].mean()); sd = float(a[c].std(ddof=0))
        a[c+'_z'] = (a[c]-mean)/sd
        scales[c] = {'mean':mean,'sd_ddof0':sd}
    eligible = a[a.grade.isin([4,5,7,8,10,11])].dropna(subset=['sex']).copy()
    retained = eligible.index[eligible.followup]
    paired = eligible.loc[retained].copy()
    for c in ENDPOINTS:
        paired[c+'_future_z'] = (b.loc[retained,c]-scales[c]['mean'])/scales[c]['sd_ddof0']
    assert len(a) == 8902 and len(eligible) == 5941 and len(paired) == 4219
    XP = dmatrix(PS_FORMULA, eligible, return_type='dataframe')
    assert len(XP) == len(eligible) and XP.index.equals(eligible.index)
    retention = float(eligible.followup.mean())
    original_p, original_info = fit_ps(XP, eligible.followup, revised=False)
    revised_p, revised_info = fit_ps(XP, eligible.followup, revised=True)
    original_weights, original_weight_info = weight_set(original_p, retained, retention)
    revised_weights, revised_weight_info = weight_set(revised_p, retained, retention)
    frozen = json.loads((EVIDENCE/'clinical_structure_summary.json').read_text())
    freeze_comparison = {}
    for key, new in [('weight_cap_p99', original_weight_info['p99_cap']),
                     ('weight_ess', original_weight_info['stabilized_p99']['ess']),
                     ('weight_min', original_weight_info['stabilized_p99']['min']),
                     ('weight_max', original_weight_info['stabilized_p99']['max'])]:
        old = frozen[key]
        freeze_comparison[key] = {'frozen':old,'current_historical_spec_refit':new,
                                  'difference':new-old,'matches_at_1e_6': bool(abs(new-old)<1e-6)}
    # Only aggregate cells are persisted; identities and individual predictions stay internal.
    cell_rows=[]
    for (grade, site), f in eligible.groupby(['grade','site'], observed=True):
        obs=f[f.followup]; ids=obs.index
        row={'grade':grade,'school':site,'baseline_n':len(f),'observed_n':int(f.followup.sum()),
             'observed_rate':float(f.followup.mean()),
             'historical_spec_mean_predicted':float(original_p.loc[f.index].mean()),
             'revised_mean_predicted':float(revised_p.loc[f.index].mean()),
             'revised_prediction_min':float(revised_p.loc[f.index].min()),
             'revised_prediction_max':float(revised_p.loc[f.index].max())}
        for name,w in revised_weights.items():
            row[name+'_weighted_n']=float(w.reindex(ids).sum())
        cell_rows.append(row)
    cells=pd.DataFrame(cell_rows)
    cells.to_csv(OUT/'ipw_grade_school_diagnostics.csv',index=False)
    grade_rows=[]
    for grade,f in eligible.groupby('grade',observed=True):
        ids=f.index[f.followup]
        row={'grade':grade,'baseline_n':len(f),'observed_n':int(f.followup.sum()),
             'observed_rate':float(f.followup.mean()),'baseline_share':len(f)/len(eligible),
             'observed_share':len(ids)/len(paired),
             'revised_mean_predicted':float(revised_p.loc[f.index].mean())}
        for name,w in revised_weights.items():
            row[name+'_weighted_share']=float(w.reindex(ids).sum()/w.sum())
        grade_rows.append(row)
    pd.DataFrame(grade_rows).to_csv(OUT/'ipw_grade_composition.csv',index=False)
    both_weights={**{'original_spec_'+k:v for k,v in original_weights.items()},**revised_weights}
    bal=balance(eligible,paired,both_weights)
    bal.to_csv(OUT/'ipw_balance.csv',index=False)
    balance_summary={}
    for label,f in bal.groupby('specification',observed=True):
        vals=f.smd.dropna().abs()
        balance_summary[label]={'n_indicators':len(vals),'max_abs_smd':float(vals.max()),
                                'mean_abs_smd':float(vals.mean()),
                                'n_abs_smd_above_0_1':int((vals>.1).sum()),
                                'share_abs_smd_above_0_1':float((vals>.1).mean())}
    outcome,cov_difference=regress(paired,revised_weights)
    outcome.to_csv(OUT/'ipw_outcome_sensitivity.csv',index=False)
    probability_summary={}
    for label,p in [('historical_spec_refit',original_p),('revised',revised_p)]:
        probability_summary[label]={'eligible':describe(p),'observed':describe(p.loc[retained]),
                                    'unobserved':describe(p.loc[eligible.index[~eligible.followup]])}
    highschool=bal[bal.variable=='high_school_grades_10_11'].set_index('specification')
    zero_cells=cells[cells.observed_n==0]
    unsupported_mask = pd.Series(False, index=eligible.index)
    for row in zero_cells.itertuples(index=False):
        unsupported_mask |= (eligible.grade == row.grade) & (eligible.site == row.school)
    supported = eligible.loc[~unsupported_mask].copy()
    assert len(supported) == 4519 and supported.index[supported.followup].equals(paired.index)
    supported_bal = balance(supported, paired, both_weights)
    supported_bal.to_csv(OUT/'ipw_supported_target_balance.csv',index=False)
    supported_balance_summary={}
    for label,f in supported_bal.groupby('specification',observed=True):
        vals=f.smd.dropna().abs()
        supported_balance_summary[label]={'n_indicators':len(vals),
                                          'max_abs_smd':float(vals.max()),
                                          'mean_abs_smd':float(vals.mean()),
                                          'n_abs_smd_above_0_1':int((vals>.1).sum()),
                                          'share_abs_smd_above_0_1':float((vals>.1).mean())}
    supported_highschool = supported_bal[supported_bal.variable=='high_school_grades_10_11'].set_index('specification')
    target_comparison=[]
    for grade in sorted(eligible.grade.unique()):
        row={'grade':grade,'original_target_n':int((eligible.grade==grade).sum()),
             'supported_target_n':int((supported.grade==grade).sum()),
             'observed_n':int((paired.grade==grade).sum()),
             'original_target_share':float((eligible.grade==grade).mean()),
             'supported_target_share':float((supported.grade==grade).mean()),
             'observed_share':float((paired.grade==grade).mean())}
        for name,w in revised_weights.items():
            row[name+'_weighted_share']=float(w[paired.grade==grade].sum()/w.sum())
        target_comparison.append(row)
    pd.DataFrame(target_comparison).to_csv(OUT/'ipw_target_composition_comparison.csv',index=False)
    result={'analysis_date':'2026-10-12','status':'post-review exploratory sensitivity',
            'source_sha256':file_hash,
            'versions':{'python':platform.python_version(),'numpy':np.__version__,
                        'pandas':pd.__version__,'scipy':scipy.__version__,
                        'statsmodels':statsmodels.__version__,'sklearn':sklearn.__version__,
                        'patsy':patsy.__version__},
            'counts':{'wave2':len(a),'eligible_continuing_grades_mapped_sex':len(eligible),
                      'paired':len(paired),'baseline_eligible_schools':eligible.site.nunique(),
                      'paired_schools':paired.site.nunique(),'retention_fraction':retention,
                      'zero_observation_grade_school_cells':len(zero_cells),
                      'baseline_in_zero_observation_cells':int(zero_cells.baseline_n.sum()),
                      'source_complete_flags_baseline':{
                          c:eligible[c].value_counts(dropna=False).to_dict()
                          for c in ['ctq_source_complete','dass_source_complete']}},
            'dass_reconstruction_mismatches':mismatch,'scaling':scales,
            'observation_model_formula':PS_FORMULA,
            'observation_matrix':{'n':len(XP),'columns':XP.shape[1],
                                  'rank':np.linalg.matrix_rank(XP.to_numpy(float)),
                                  'condition_number':np.linalg.cond(XP.to_numpy(float))},
            'historical_spec_refit':original_info,'revised_observation_fit':revised_info,
            'original_vs_revised_prediction_max_abs_difference':float(abs(original_p-revised_p).max()),
            'probability_distribution':probability_summary,
            'historical_spec_refit_weights':original_weight_info,
            'revised_weights':revised_weight_info,
            'frozen_weight_summary_comparison':freeze_comparison,
            'frozen_reproduced':all(v['matches_at_1e_6'] for v in freeze_comparison.values()),
            'balance_summary':balance_summary,
            'supported_target_diagnostics':{
                'description':'Post-diagnostic descriptive amendment; same fitted probabilities and weights, no new model or change to outcome sample',
                'n':len(supported),'schools':supported.site.nunique(),
                'retention_fraction':float(supported.followup.mean()),
                'excluded_zero_observation_cells':zero_cells[['grade','school','baseline_n','observed_n']].to_dict(orient='records'),
                'school_absence_checks':[
                    {'school':school,
                     'all_wave2_records':int((a.site==school).sum()),
                     'all_wave2_records_with_any_wave3_match':int(a.loc[a.site==school,'followup'].sum()),
                     'wave3_records_in_same_school':int((b.site==school).sum())}
                    for school in sorted(zero_cells.school.unique())],
                'balance_summary':supported_balance_summary,
                'high_school_share':{'baseline':float(supported_highschool.target_mean.iloc[0]),
                                     **{k:float(row.observed_mean) for k,row in supported_highschool.iterrows()}}},
            'high_school_share':{'baseline':float(highschool.target_mean.iloc[0]),
                                  **{k:float(row.observed_mean) for k,row in highschool.iterrows()}},
            'independent_checks':{'manual_weight_sum_and_ess_agree':True,
                                  'explicit_weight_index_alignment_passed':True,
                                  'unstabilized_stabilized_wls_invariance_passed':True,
                                  'manual_cr1_vs_statsmodels_se_max_abs_difference':cov_difference,
                                  'source_hash_unchanged':hashlib.sha256(SOURCE.read_bytes()).hexdigest()==file_hash},
            'outcome_rows':outcome.to_dict(orient='records'),
            'interpretation_boundary':'Fixed-weight school-CR1 intervals do not account for observation-model estimation; measured-selection sensitivity only, not MAR verification or a causal analysis; no support for entirely absent grades/cells.'}
    (OUT/'ipw_summary.json').write_text(json.dumps(clean_json(result),indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(clean_json({'counts':result['counts'],'frozen_reproduced':result['frozen_reproduced'],
                                'revised_weights':result['revised_weights'],
                                'balance_summary':balance_summary,
                                'supported_target_diagnostics':result['supported_target_diagnostics'],
                                'high_school_share':result['high_school_share'],
                                'distress':outcome[outcome.endpoint=='distress'].to_dict(orient='records')}),indent=2))


if __name__=='__main__':
    sys.path.insert(0, str(HERE.parent))
    from revision_runtime import configure_output
    OUT = configure_output(ROOT, 'Re-run fixed October 12 observation-weight diagnostics.')
    main()
