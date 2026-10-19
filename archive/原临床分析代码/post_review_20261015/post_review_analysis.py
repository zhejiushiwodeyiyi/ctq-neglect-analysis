"""Fixed 15 October 2026 post-review analyses; outputs are aggregate only."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
sys.dont_write_bytecode=True

import numpy as np
import pandas as pd
import patsy
from patsy import dmatrix, build_design_matrices
import scipy
from scipy import stats, linalg
import statsmodels
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
OUT = HERE / 'reproduction_output'
ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / 'data/questionnaire_analysis.csv'
EVIDENCE = HERE / 'reference_results'
sys.path.insert(0,str(HERE.parent))
spec = importlib.util.spec_from_file_location('prior_content', HERE/'dependencies/content_analysis_20261012.py')
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
prior.ROOT = ROOT
prior.SOURCE = SOURCE
prior.EVIDENCE = EVIDENCE
prior.SCALING_REFERENCE = EVIDENCE / 'manuscript_scaling.json'
threadpool_limits(limits=1)
ITEMS = prior.ITEMS
DOMAINS = prior.DOMAINS
EXPOSURES = ['neglect', 'abuse']
SAVED_MODELS = []


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fit_matrix(x, y, school, label):
    a=x.to_numpy(float)
    _,rr,pivot=linalg.qr(a,mode='economic',pivoting=True)
    tolerance=max(a.shape)*np.finfo(float).eps*np.abs(rr).max()
    rank=int((np.abs(np.diag(rr))>tolerance).sum())
    kept=sorted(pivot[:rank])
    protected=[i for i,c in enumerate(x.columns) if 'neglect_z' in c or 'abuse_z' in c]
    assert set(protected).issubset(kept), 'An exposure basis is aliased; do not silently change the estimand.'
    x=x.iloc[:,kept]
    assert np.linalg.matrix_rank(x.to_numpy()) == x.shape[1]
    fitted = sm.OLS(np.asarray(y), x).fit(
        cov_type='cluster', cov_kwds={'groups': school, 'use_correction': True}, use_t=True)
    SAVED_MODELS.append(dict(model_id=label, n_rows=len(x), parameters=x.shape[1],
                             schools=int(pd.Series(school).nunique())))
    return x, fitted


def secondary_models(d):
    second = d[d.grade.ge(7)].copy()
    assert len(second) == 2406 and second.site.nunique() == 7
    stage_out = OUT / 'secondary_implementation'
    stage_out.mkdir(exist_ok=True)
    prior.OUT = stage_out
    result, checks = prior.content_models(second)
    result['sample'] = 'baseline_grade_ge7'
    result['reference'] = 'full_eligible_wave2_n8902_ddof0'
    result['q_wild_4_new_contrasts'] = prior.bh(result.p_wild)
    result = result.rename(columns={'q_wild_2': 'q_wild_2_within_units'})
    result.to_csv(OUT / 'secondary_common_prior_contrasts.csv', index=False)
    components = pd.read_csv(stage_out / 'content_revision_common_baseline_slopes.csv')
    components['sample'] = 'baseline_grade_ge7'
    components.to_csv(OUT / 'secondary_common_prior_components.csv', index=False)
    assert result.wild_patterns.eq(128).all()
    assert np.allclose(result.p_wild * 128, result.bootstrap_exceedance_count)
    for r in result.to_dict('records'):
        SAVED_MODELS.append(dict(model_id=f"secondary_direct_{r['units']}_{r['exposure']}",
                                 n_rows=4812, schools=7, parameters=r['matrix_columns']))
    for r in components.to_dict('records'):
        SAVED_MODELS.append(dict(model_id=f"secondary_component_{r['units']}_{r['exposure']}_{r['content']}",
                                 n_rows=2406, schools=7, parameters=r['matrix_columns']))
    return result, checks


def scale_register():
    m = pd.read_csv(SOURCE, low_memory=False)
    m['neglect'] = m[prior.SUB[:2]].sum(axis=1, min_count=2)
    m['abuse'] = m[prior.SUB[2:]].sum(axis=1, min_count=2)
    m['distress'] = m[DOMAINS].mean(axis=1)
    for name, numbers in ITEMS.items():
        columns = [f'dass{i}' for i in numbers]
        m[name] = m[columns].sum(axis=1, min_count=len(columns)) / len(columns)
    rows = []
    school_refs = {}
    for first, last in [(1, 2), (2, 3)]:
        a = m[m.wave.eq(first) & m.ctq_source_complete & m.dass_source_complete]
        school_refs[first] = a
        for variable in EXPOSURES + ['distress'] + DOMAINS + prior.SUB + list(ITEMS):
            v = a[variable].dropna()
            rows.append(dict(population='school', interval=f'{first}-{last}',
                             reference=f'full_eligible_complete_wave{first}_baseline',
                             reference_population_n=len(a), variable_n=len(v), variable=variable,
                             mean=v.mean(), sd=v.std(ddof=0), ddof=0,
                             provenance='recomputed_from_corrected_deidentified_table',
                             notes='same reference retained in paired, restricted and follow-up values'))
    adult = json.loads((EVIDENCE / 'external_clinical_scaling.json').read_text())
    adult_sample = pd.read_csv(EVIDENCE / 'external_clinical_sample.csv').set_index('cohort')
    for cohort, variables in adult.items():
        for variable, scale in variables.items():
            rows.append(dict(population='adult', interval='baseline-to-two-year',
                             reference=cohort + '_available_baseline_before_complete_case_selection',
                             reference_population_n=int(adult_sample.loc[cohort, 'n']),
                             variable_n=np.nan, variable=variable, mean=scale['mean'], sd=scale['sd'],
                             ddof=0, provenance='existing_historical_aggregate_scaling',
                             notes='variable-specific nonmissing N not recorded in supplied aggregate; adult input not independently reread'))
    full = pd.DataFrame(rows)
    full.to_csv(OUT / 'scaling_register_full.csv', index=False)
    keep = (full.population.eq('school') & full.variable.isin(EXPOSURES + ['distress'] + list(ITEMS))) | (
        full.population.eq('adult') & full.variable.isin(['neglect', 'abuse', 'prior']))
    full[keep].to_csv(OUT / 'scaling_register_manuscript.csv', index=False)
    return school_refs, full


def formula_for(form, system, knots):
    continuous = ['neglect_z', 'abuse_z'] + (['distress_z'] if system == 'distress' else [c+'_z' for c in DOMAINS])
    terms = []
    for variable in continuous:
        if form == 'linear':
            terms.append(variable)
        else:
            k = knots[variable]
            terms.append(f"cr({variable}, knots=({k['inner'][0]!r}, {k['inner'][1]!r}), lower_bound={k['lower']!r}, upper_bound={k['upper']!r}, constraints='center')")
    rhs = '+'.join(terms) + '+C(grade)*sex+C(site)'
    return rhs if system == 'distress' else 'content*(' + rhs + ')'


def clustered_contrast(fitted, x, school, row, label):
    v = np.asarray(row, dtype=float)
    estimate = float(v @ fitted.params)
    variance = float(v @ fitted.cov_params().to_numpy() @ v)
    se = np.sqrt(max(0., variance))
    groups = int(pd.Series(school).nunique())
    critical = stats.t.ppf(.975, groups - 1)
    return dict(model_id=label, estimate=estimate, se_school=se,
                school_low=estimate-critical*se, school_high=estimate+critical*se,
                schools=groups, school_df=groups-1, n_rows=len(x), parameters=x.shape[1])


def design_diff(info, frame, variable, lo, hi):
    lower, upper = frame.copy(), frame.copy()
    lower[variable], upper[variable] = lo, hi
    xl, xu = build_design_matrices([info], lower)[0], build_design_matrices([info], upper)[0]
    return np.asarray(xu).mean(axis=0) - np.asarray(xl).mean(axis=0)


def influence_summary(x, y, school, participant, vectors, system, form):
    a, yy = x.to_numpy(float), np.asarray(y, float)
    n, p = a.shape
    pinv = np.linalg.pinv(a)
    xtxi = pinv @ pinv.T
    beta = pinv @ yy
    error = yy - a @ beta
    leverage = np.einsum('ij,ji->i', a, pinv)
    rss = float(error @ error)
    group_codes = pd.Categorical(participant).codes
    unique_codes = np.unique(group_codes)
    group_rows = [np.flatnonzero(group_codes == c) for c in unique_codes]
    max_leverage = np.array([leverage[rows].max() for rows in group_rows])
    rows_out = []
    severe = np.zeros(len(group_rows), dtype=bool)
    cfit = sm.OLS(yy, x).fit(cov_type='cluster', cov_kwds={'groups':school, 'use_correction':True}, use_t=True)
    for exposure, v in vectors.items():
        v = np.asarray(v, float)
        denominator = np.sqrt(max(float(v @ xtxi @ v), 0.))
        changes, dfbetas = [], []
        for rr in group_rows:
            ai, ei = a[rr], error[rr]
            block = np.eye(len(rr)) - ai @ xtxi @ ai.T
            adjusted = np.linalg.solve(block, ei)
            delta = xtxi @ ai.T @ adjusted
            change = float(v @ delta)
            rss_deleted = max(0., rss - float(ei @ adjusted))
            sigma_deleted = np.sqrt(rss_deleted / (n - p - len(rr)))
            dfbeta = change / (sigma_deleted * denominator)
            changes.append(change)
            dfbetas.append(dfbeta)
        changes, dfbetas = np.asarray(changes), np.asarray(dfbetas)
        original = float(v @ beta)
        cluster_se = np.sqrt(float(v @ cfit.cov_params().to_numpy() @ v))
        severe |= (max_leverage > .5) | (np.abs(dfbetas) > 1.)
        rows_out.append(dict(system=system, form=form, target=exposure,
                             participant_n=len(group_rows), outcome_rows=n, parameters=p,
                             deletion_unit='participant_all_content_rows' if system=='content' else 'participant',
                             original_estimate=original, se_school=cluster_se,
                             max_row_leverage=float(leverage.max()),
                             leverage_p99=float(np.quantile(max_leverage,.99)),
                             leverage_gt_2p_over_n_n=int((max_leverage>2*p/n).sum()),
                             max_abs_group_dfbeta=float(np.abs(dfbetas).max()),
                             abs_group_dfbeta_p99=float(np.quantile(np.abs(dfbetas),.99)),
                             abs_group_dfbeta_gt_2_over_sqrt_participant_n=int((np.abs(dfbetas)>2/np.sqrt(len(group_rows))).sum()),
                             max_abs_single_participant_delete_change=float(np.abs(changes).max()),
                             max_delete_change_over_school_se=float(np.abs(changes).max()/cluster_se),
                             minimum_single_delete_estimate=float((original-changes).min()),
                             maximum_single_delete_estimate=float((original-changes).max())))
    severe_codes = set(unique_codes[severe])
    keep = np.array([c not in severe_codes for c in group_codes])
    for row in rows_out:
        row['severe_participant_union_n'] = int(severe.sum())
    return rows_out, keep


def function_form(baseline, d, scales):
    knots = {}
    variables = EXPOSURES + ['distress'] + DOMAINS
    for variable in variables:
        values = np.sort(baseline[variable].dropna().unique())
        standardized = (values - scales[variable]['mean']) / scales[variable]['sd']
        lower, upper = float(standardized.min()), float(standardized.max())
        inner = np.quantile(standardized, [1/3, 2/3]).tolist()
        assert lower < inner[0] < inner[1] < upper
        knots[variable+'_z'] = dict(lower=lower, inner=inner, upper=upper,
                                    unique_value_n=len(values),
                                    raw_lower=float(values.min()), raw_inner=np.quantile(values,[1/3,2/3]).tolist(), raw_upper=float(values.max()),
                                    rule='terciles_of_sorted_unique_full_baseline_values')
    (OUT/'spline_knots.json').write_text(json.dumps(knots, indent=2))
    rows, predictions, influences, severe_rows = [], [], [], []
    context = {}
    for form in ['linear', 'natural_cubic_spline_3df']:
        for system, units in [('distress','standardized'), ('content','standardized'), ('content','raw_0to3')]:
            base = d.copy()
            if system == 'distress':
                f = base.copy(); f['outcome'] = f.distress_future_z
                f['participant'] = np.arange(len(base))
            else:
                fs=[]
                for content in ITEMS:
                    ff=base.copy()
                    ff['outcome']=ff[content+('_future_z' if units=='standardized' else '_future')]
                    ff['content']=int(content=='anhedonic_content')
                    ff['participant']=np.arange(len(base))
                    fs.append(ff)
                f=pd.concat(fs,ignore_index=True)
            formula = formula_for(form, system, knots)
            raw_x = dmatrix(formula, f, return_type='dataframe')
            info = raw_x.design_info
            label=f'{system}_{units}_{form}'
            x, fitted = fit_matrix(raw_x, f.outcome, f.site, label)
            kept = [raw_x.columns.get_loc(c) for c in x.columns]
            vectors = {}
            for exposure in EXPOSURES:
                lo, hi = baseline[exposure].quantile([.25,.75])
                lo_z, hi_z = (lo-scales[exposure]['mean'])/scales[exposure]['sd'], (hi-scales[exposure]['mean'])/scales[exposure]['sd']
                if system=='distress':
                    contrast=design_diff(info, f, exposure+'_z',lo_z,hi_z)[kept]
                else:
                    aa=f[f.content.eq(1)]; bb=f[f.content.eq(0)]
                    contrast=(design_diff(info,aa,exposure+'_z',lo_z,hi_z)-design_diff(info,bb,exposure+'_z',lo_z,hi_z))[kept]
                vectors[exposure]=contrast
                rows.append(dict(system=system, outcome_units=units, function_form=form, exposure=exposure,
                                 estimand='conditional_Q25_to_Q75_difference' if system=='distress' else 'anhedonic_minus_autonomic_conditional_Q25_to_Q75_difference',
                                 raw_exposure_q25=lo, raw_exposure_q75=hi,
                                 exposure_sd_steps=hi_z-lo_z,
                                 **clustered_contrast(fitted,x,f.site,contrast,label)))
            predicted=np.asarray(fitted.predict(x))
            if system=='distress':
                raw_pred=predicted*scales['distress']['sd']+scales['distress']['mean']
                predictions.append(dict(model_id=label,outcome='distress',n=len(f),
                                        observed_raw_min=float(f.distress_future.min()),observed_raw_max=float(f.distress_future.max()),
                                        fitted_raw_min=float(raw_pred.min()),fitted_raw_max=float(raw_pred.max()),
                                        fitted_below_scale_n=int((raw_pred<0).sum()),fitted_above_scale_n=int((raw_pred>42).sum()),
                                        scale_lower=0,scale_upper=42))
            else:
                for content, flag in [('anhedonic_content',1),('autonomic_content',0)]:
                    mask=f.content.eq(flag).to_numpy()
                    raw_pred=(predicted[mask]*scales[content]['sd']+scales[content]['mean']) if units=='standardized' else predicted[mask]
                    predictions.append(dict(model_id=label,outcome=content,n=int(mask.sum()),
                                            observed_raw_min=float(f.loc[mask,content+'_future'].min()),observed_raw_max=float(f.loc[mask,content+'_future'].max()),
                                            fitted_raw_min=float(raw_pred.min()),fitted_raw_max=float(raw_pred.max()),
                                            fitted_below_scale_n=int((raw_pred<0).sum()),fitted_above_scale_n=int((raw_pred>3).sum()),
                                            scale_lower=0,scale_upper=3))
            influence, keep=influence_summary(x,f.outcome,f.site,f.participant,vectors,system,form)
            for row in influence: row['outcome_units']=units
            influences.extend(influence)
            deleted=int(f.loc[~keep,'participant'].nunique())
            if deleted:
                xx, fit2=fit_matrix(x.loc[keep],f.loc[keep,'outcome'],f.loc[keep,'site'],label+'_severe_deletion')
                for exposure,contrast in vectors.items():
                    assert list(xx.columns)==list(x.columns)
                    severe_rows.append(dict(model_id=label,exposure=exposure,deleted_participants=deleted,
                                             **clustered_contrast(fit2,xx,f.loc[keep,'site'],contrast,label+'_severe_deletion')))
            if form=='linear' and system=='content':
                context[units]=(x,fitted,f,info,kept)
    pd.DataFrame(rows).to_csv(OUT/'function_form_iqr_contrasts.csv',index=False)
    pd.DataFrame(predictions).to_csv(OUT/'function_form_prediction_bounds.csv',index=False)
    pd.DataFrame(influences).to_csv(OUT/'influence_aggregate.csv',index=False)
    pd.DataFrame(severe_rows).to_csv(OUT/'severe_influence_sensitivity.csv',index=False)
    return pd.DataFrame(rows),pd.DataFrame(influences),context


def covariance_identities(context, scales):
    x,fit,f,_,_=context['raw_0to3']
    xs,fits,_,_,_=context['standardized']
    rows=[]
    stored=pd.read_csv(EVIDENCE/'content_revision_common_baseline_contrasts.csv')
    for exposure in EXPOSURES:
        u=np.zeros(x.shape[1]);v=u.copy()
        u[x.columns.get_loc(exposure+'_z')]=1
        v[x.columns.get_loc(exposure+'_z')]=1
        v[x.columns.get_loc('content:'+exposure+'_z')]=1
        components=np.vstack([v,u])
        betas=components @ fit.params.to_numpy()
        covariance=components @ fit.cov_params().to_numpy() @ components.T
        s=np.array([scales[c]['sd'] for c in ITEMS])
        standardized=float(betas[0]/s[0]-betas[1]/s[1])
        standardized_var=float(np.array([1/s[0],-1/s[1]])@covariance@np.array([1/s[0],-1/s[1]]))
        actual=float(fits.params['content:'+exposure+'_z'])
        actual_var=float(fits.cov_params().loc['content:'+exposure+'_z','content:'+exposure+'_z'])
        raw=float(betas[0]-betas[1])
        oldstd=stored[(stored.exposure==exposure)&stored.units.eq('standardized')].iloc[0]
        oldraw=stored[(stored.exposure==exposure)&stored.units.eq('raw_0to3')].iloc[0]
        error=max(abs(standardized-actual),abs(standardized_var-actual_var),abs(standardized-oldstd.beta),abs(raw-oldraw.beta))
        assert error<1e-10
        rows.append(dict(exposure=exposure,raw_anhedonic_slope=betas[0],raw_autonomic_slope=betas[1],
                         anhedonic_reference_sd=s[0],autonomic_reference_sd=s[1],
                         raw_anhedonic_variance=covariance[0,0],raw_autonomic_variance=covariance[1,1],
                         raw_joint_covariance=covariance[0,1],raw_contrast=raw,
                         standardized_contrast=standardized,standardized_variance_from_joint_covariance=standardized_var,
                         standardized_se=np.sqrt(standardized_var),maximum_identity_error=error))
    pd.DataFrame(rows).to_csv(OUT/'common_prior_joint_covariance_identities.csv',index=False)
    # Stored own-prior identities for both intervals use their actual baseline-wave references.
    slopes=pd.read_csv(EVIDENCE/'content_raw_unit_slopes.csv')
    direct=pd.read_csv(EVIDENCE/'source_clinical_content_contrast.csv')
    identity=[]
    for (interval,exposure),group in slopes.groupby(['transition','exposure']):
        ss=group.set_index('content')
        std=ss.raw_beta/ss.baseline_content_sd
        implied=float(std['anhedonic_content']-std['autonomic_content'])
        old=direct[(direct.transition==interval)&(direct.exposure==exposure)].iloc[0]
        err=abs(implied-old.beta_difference)
        assert err<1e-10
        identity.append(dict(interval=interval,exposure=exposure,adjustment='different_own_prior_content',
                             implied_standardized_difference=implied,stored_standardized_difference=old.beta_difference,
                             absolute_identity_error=err))
    pd.DataFrame(identity).to_csv(OUT/'interval_own_prior_scaling_identities.csv',index=False)
    return rows,identity


def subscale_contrast(d):
    # This optional contrast tests specified per-subscale-SD conditional slopes, not contributions.
    f=d.copy()
    formula='+'.join(c+'_z' for c in prior.SUB)+'+distress_z+C(grade)*sex+C(site)'
    x,fitted=fit_matrix(dmatrix(formula,f,return_type='dataframe'),f.distress_future_z,f.site,'subscale_EN_minus_PN')
    contrast=np.zeros(x.shape[1])
    en,pn=prior.SUB[:2]
    contrast[x.columns.get_loc(en+'_z')]=1
    contrast[x.columns.get_loc(pn+'_z')]=-1
    row=clustered_contrast(fitted,x,f.site,contrast,'subscale_EN_minus_PN')
    # Reparameterize so ctq_z is exactly b_EN-b_PN; preserve the fitted column space.
    f['ctq_z']=f[en+'_z'];f['sum_en_pn']=f[en+'_z']+f[pn+'_z']
    wild_formula='ctq_z+sum_en_pn+'+'+'.join(c+'_z' for c in prior.SUB[2:])+'+distress_z+C(grade)*sex+C(site)'
    xx,check=prior.school_fit(dmatrix(wild_formula,f,return_type='dataframe'),f.distress_future_z,f.site)
    SAVED_MODELS.append(dict(model_id='subscale_EN_minus_PN_reparameterization',n_rows=len(xx),schools=int(f.site.nunique()),parameters=xx.shape[1]))
    bootstrap=prior.wild_matrix(xx,f.distress_future_z,f.site)
    assert abs(row['estimate']-check['beta'])<1e-10
    pd.DataFrame([dict(contrast='EN_minus_PN_per_own_full_baseline_SD',n_participants=len(f),
                       **row,**bootstrap)]).to_csv(OUT/'subscale_EN_minus_PN_contrast.csv',index=False)
    return dict(**row,**bootstrap)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    baseline,future,d,scales,scale_error=prior.prepare()
    assert abs(scales['anhedonic_content']['sd']-.529394860062179)<1e-12
    assert abs(scales['autonomic_content']['sd']-.48290090067819813)<1e-12
    print('Inputs prepared and fixed scaling confirmed.',flush=True)
    secondary,second_checks=secondary_models(d)
    print(secondary[['units','exposure','beta','school_low','school_high','p_wild','q_wild_4_new_contrasts']].to_string(index=False),flush=True)
    refs,scaling=scale_register()
    print('Scaling register prepared.',flush=True)
    functions,influence,contexts=function_form(baseline,d,scales)
    print(functions[['system','outcome_units','function_form','exposure','estimate','school_low','school_high']].to_string(index=False),flush=True)
    identities,ownidentities=covariance_identities(contexts,scales)
    subscale=subscale_contrast(d)
    print('Influence and algebraic checks completed.',flush=True)
    pd.DataFrame(SAVED_MODELS).to_csv(OUT/'new_model_manifest.csv',index=False)
    audit=dict(plan='POST_REVIEW_PLAN_20261015.md',analysis_date='2026-10-15',post_review_not_preregistered=True,
               source_sha256=sha(SOURCE),source_path=str(SOURCE),plan_sha256=sha(HERE/'POST_REVIEW_PLAN_20261015.md'),
               prior_source_sha256=sha(HERE/'dependencies/content_analysis_20261012.py'),script_sha256=sha(__file__),
               dependency_sha256={str(p.relative_to(ROOT)):sha(p) for p in list((HERE/'dependencies').glob('*.py'))+list(EVIDENCE.glob('*'))},
               versions=dict(python=sys.version.split()[0],numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__,patsy=patsy.__version__),
               baseline_n=len(baseline),paired_n=len(d),paired_schools=d.site.nunique(),
               secondary_n=2406,secondary_schools=7,secondary_wild_patterns=128,
               full_reference_scaling_max_abs_error=scale_error,
               max_secondary_stacking_error=max(c['largest_absolute_error'] for c in second_checks),
               max_covariance_identity_error=max(c['maximum_identity_error'] for c in identities),
               max_own_prior_scaling_identity_error=max(c['absolute_identity_error'] for c in ownidentities),
               observed_outcome_fit_manifest_rows=len(SAVED_MODELS),
               severe_influence_max_union_n=int(influence.severe_participant_union_n.max()),
               adult_variable_specific_reference_n='unavailable_in_aggregate_records',
               boundaries=['No private raw-questionnaire source adjudication repeated.',
                           'No adult individual inputs reread; adult scales are frozen historical aggregates.',
                           'No participant identifiers or individual observations exported.',
                           'Diagnostic spline intervals are school CR1 with t(G-1), not new wild-bootstrap significance tests.',
                           'Post-review analyses are not independent replication; diagnostic contrasts are conditional rather than causal.'])
    (OUT/'post_review_analysis_audit.json').write_text(json.dumps(audit,indent=2))
    (OUT/'analysis_summary.json').write_text(json.dumps(dict(secondary=secondary.to_dict('records'),function_form=functions.to_dict('records'),subscale_contrast=subscale),indent=2))
    print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':
    from revision_runtime import configure_output
    OUT = configure_output(ROOT, 'Re-run fixed October 15 post-review analyses.')
    main()
