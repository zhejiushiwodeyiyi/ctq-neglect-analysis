"""Bounded post-result diagnostic addendum; no changes to the initial result set."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import patsy
from patsy import dmatrix
import scipy
import statsmodels

HERE=Path(__file__).resolve().parent
ROOT=Path(__file__).resolve().parents[3]
ORIGINAL=ROOT/'paper/evidence/results/revision_20261015'
OUT=ROOT/'reproduction_output/diagnostic_addendum_20261015'
PLAN=HERE/'POST_REVIEW_DIAGNOSTIC_ADDENDUM_20261015.md'
spec=importlib.util.spec_from_file_location('first_postreview',HERE/'post_review_analysis.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    OUT.mkdir(exist_ok=True)
    frozen={str(p.resolve()):sha(p) for p in ORIGINAL.rglob('*') if p.is_file()}
    for p in [HERE/'post_review_analysis.py',HERE/'POST_REVIEW_PLAN_20261015.md']:
        frozen[str(p.resolve())]=sha(p)
    (OUT/'original_outputs_frozen_before.json').write_text(json.dumps(frozen,indent=2))
    baseline,_,d,scales,_=base.prior.prepare()
    assert len(baseline)==8902 and len(d)==4219 and d.site.nunique()==10
    knots=json.loads((ORIGINAL/'spline_knots.json').read_text())
    old=pd.read_csv(ORIGINAL/'function_form_iqr_contrasts.csv')
    frame=d.copy();frame['outcome']=frame.distress_future_z
    q25,q75,p95=map(float,baseline.neglect.quantile([.25,.75,.95]))
    grid=np.arange(10,int(np.floor(p95))+1,dtype=int)
    assert q25==10 and q75==17 and len(grid)>1
    iqr,curves,models=[],[],[]
    reconstruction_errors=[]
    specs=[('exposure_splines_baseline_linear',True,False),
           ('exposure_linear_baseline_spline',False,True),
           ('exposure_and_baseline_splines',True,True)]
    for label,exposure_spline,prior_spline in specs:
        terms=[]
        for variable in ['neglect_z','abuse_z','distress_z']:
            flexible=exposure_spline if variable!='distress_z' else prior_spline
            if flexible:
                k=knots[variable]
                terms.append(f"cr({variable}, knots=({k['inner'][0]!r},{k['inner'][1]!r}), lower_bound={k['lower']!r}, upper_bound={k['upper']!r}, constraints='center')")
            else:terms.append(variable)
        rawx=dmatrix('+'.join(terms)+'+C(grade)*sex+C(site)',frame,return_type='dataframe')
        info=rawx.design_info
        x,fit=base.fit_matrix(rawx,frame.outcome,frame.site,label)
        kept=[rawx.columns.get_loc(c) for c in x.columns]
        models.append(dict(model_id=label,n_participants=len(frame),schools=10,parameters=x.shape[1],
                           purpose='reconstruction_for_curve' if exposure_spline and prior_spline else 'new_mixed_diagnostic_specification'))
        for exposure in ['neglect','abuse']:
            lower,upper=map(float,baseline[exposure].quantile([.25,.75]))
            lo=(lower-scales[exposure]['mean'])/scales[exposure]['sd']
            hi=(upper-scales[exposure]['mean'])/scales[exposure]['sd']
            v=base.design_diff(info,frame,exposure+'_z',lo,hi)[kept]
            row=dict(specification=label,exposure=exposure,raw_q25=lower,raw_q75=upper,
                     **base.clustered_contrast(fit,x,frame.site,v,label))
            iqr.append(row)
            if exposure_spline and prior_spline:
                prev=old[(old.system=='distress')&(old.function_form=='natural_cubic_spline_3df')&(old.exposure==exposure)].iloc[0]
                err=max(abs(row[c]-prev[c]) for c in ['estimate','se_school','school_low','school_high'])
                assert err<1e-10,err
                reconstruction_errors.append(dict(exposure=exposure,maximum_error=err))
        ref=(10-scales['neglect']['mean'])/scales['neglect']['sd']
        for score in grid:
            target=(score-scales['neglect']['mean'])/scales['neglect']['sd']
            v=base.design_diff(info,frame,'neglect_z',ref,target)[kept]
            row=base.clustered_contrast(fit,x,frame.site,v,label)
            curves.append(dict(specification=label,neglect_score=int(score),reference_score=10,
                               outcome='conditional_general_distress_difference_in_full_W2_SD',
                               n_participants=len(frame),interval_type='pointwise_school_CR1_t9',**row))
    linear=old[(old.system=='distress')&(old.function_form=='linear')&(old.exposure=='neglect')].iloc[0]
    for score in grid:
        factor=(score-10)/(q75-q25)
        curves.append(dict(specification='all_linear',neglect_score=int(score),reference_score=10,
                           outcome='conditional_general_distress_difference_in_full_W2_SD',n_participants=len(frame),
                           interval_type='pointwise_school_CR1_t9',model_id='frozen_linear_IQR_algebra',
                           estimate=float(linear.estimate*factor),se_school=float(linear.se_school*factor),
                           school_low=float(linear.school_low*factor),school_high=float(linear.school_high*factor),
                           schools=10,school_df=9,n_rows=4219,parameters=int(linear.parameters)))
    curve=pd.DataFrame(curves)
    assert (curve[curve.neglect_score.eq(10)][['estimate','se_school','school_low','school_high']].to_numpy()==0).all()
    curve.sort_values(['specification','neglect_score']).to_csv(OUT/'neglect_conditional_difference_curves.csv',index=False)
    pd.DataFrame(iqr).to_csv(OUT/'mixed_specification_iqr_contrasts.csv',index=False)
    pd.DataFrame(models).to_csv(OUT/'addendum_model_manifest.csv',index=False)
    distribution=baseline.neglect.value_counts().sort_index()
    pd.DataFrame({'neglect_score':distribution.index.astype(int),'n':distribution.to_numpy(),
                  'percent':100*distribution.to_numpy()/len(baseline)}).to_csv(OUT/'neglect_full_W2_aggregate_distribution.csv',index=False)
    assert all(sha(Path(rel))==expected for rel,expected in frozen.items())
    output=dict(date='2026-10-15',plan_sha256=sha(PLAN),script_sha256=sha(__file__),
                source_sha256=sha(base.SOURCE),knots_sha256=sha(ORIGINAL/'spline_knots.json'),
                written_after_first_diagnostic_results=True,original_outputs_unchanged=True,
                n_new_mixed_model_fits=2,n_necessary_model_reconstructions=1,n_linear_refits=0,
                reconstruction_errors=reconstruction_errors,full_W2_n=8902,paired_n=4219,schools=10,
                neglect_q25=q25,neglect_q75=q75,neglect_p95=p95,integer_grid_min=10,integer_grid_max=int(grid.max()),
                versions=dict(numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__,patsy=patsy.__version__),
                curve_intervals='pointwise school CR1 with t(9), not simultaneous; fixed-reference prediction design-vector covariance',
                participant_identifiers_or_observations_exported=False,
                no_p_value_or_model_selection_search=True,
                output_sha256={p.name:sha(p) for p in OUT.glob('*') if p.is_file() and p.name!='diagnostic_addendum_audit.json'})
    (OUT/'diagnostic_addendum_audit.json').write_text(json.dumps(output,indent=2))
    print(pd.DataFrame(iqr)[['specification','exposure','estimate','school_low','school_high']].to_string(index=False),flush=True)
    print(json.dumps(output,indent=2),flush=True)

if __name__=='__main__':
    import argparse
    sys.path.insert(0, str(HERE.parent))
    from revision_runtime import ensure_external_output
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--initial-results',type=Path,default=ORIGINAL)
    args=parser.parse_args()
    OUT=ensure_external_output(ROOT,args.output)
    ORIGINAL=args.initial_results.resolve()
    if not (ORIGINAL/'spline_knots.json').is_file():parser.error('Initial results lack spline_knots.json')
    main()
