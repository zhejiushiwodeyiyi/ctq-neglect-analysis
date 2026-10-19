"""Plot frozen aggregate conditional contrasts; never read participant data."""
from pathlib import Path
import hashlib
import json
import shutil
import os

ROOT=Path(__file__).resolve().parents[3]
import argparse,sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from revision_runtime import ensure_external_output
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--source',type=Path,default=ROOT/'paper/evidence/results/diagnostic_addendum_20261015')
args=parser.parse_args()
OUT=ensure_external_output(ROOT,args.output)
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'matplotlib_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

SOURCE=args.source.resolve()
curves=pd.read_csv(SOURCE/'neglect_conditional_difference_curves.csv')
freq=pd.read_csv(SOURCE/'neglect_full_W2_aggregate_distribution.csv')
assert len(curves)==92 and int(freq.n.sum())==8902
assert curves.n_participants.eq(4219).all()
assert curves.schools.eq(10).all() and curves.school_df.eq(9).all()
at10=curves[curves.neglect_score.eq(10)]
assert np.max(np.abs(at10[['estimate','school_low','school_high']].to_numpy()))<1e-12

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.5,'axes.labelsize':10,
                     'axes.titlesize':10,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none',
                     'axes.spines.top':False,'axes.spines.right':False})
fig,(ax,bx)=plt.subplots(2,1,figsize=(7.2,5.7),gridspec_kw={'height_ratios':[2.25,1]},
                         constrained_layout=True)
specs=[('all_linear','Linear exposures and baseline','#0072B2','-',1.7),
       ('exposure_linear_baseline_spline','Linear exposures; spline baseline','#009E73','--',1.2),
       ('exposure_splines_baseline_linear','Spline exposures; linear baseline','#CC79A7',':',1.7),
       ('exposure_and_baseline_splines','Spline exposures and baseline','#D55E00','-',1.9)]
for key,label,color,linestyle,width in specs:
    d=curves[curves.specification.eq(key)].sort_values('neglect_score')
    assert len(d)==23 and d.neglect_score.min()==10 and d.neglect_score.max()==32
    if key in ['all_linear','exposure_and_baseline_splines']:
        ax.fill_between(d.neglect_score,d.school_low,d.school_high,color=color,alpha=.11,linewidth=0)
    ax.plot(d.neglect_score,d.estimate,color=color,linestyle=linestyle,linewidth=width,label=label)
ax.axhline(0,color='#666666',linewidth=.7,zorder=0)
ax.axvline(17,color='#888888',linewidth=.7,linestyle=(0,(2,3)),zorder=0)
ax.text(17.15,.55,'Q75 = 17',transform=ax.get_xaxis_transform(),va='center',ha='left',fontsize=8,color='#666666')
ax.set(xlim=(10,32),xticks=[10,15,20,25,30,32],xlabel='Baseline CTQ neglect score (curve range 10 to P95 = 32)',
       ylabel='Conditional general-distress difference\nfrom neglect score 10 (baseline SD)')
ax.legend(loc='upper left',frameon=False,borderaxespad=.5,handlelength=2.8)
ax.set_title('a  Model-dependent conditional differences',loc='left',fontweight='bold',pad=8)

bx.bar(freq.neglect_score,freq.percent,width=.88,color='#71889B',edgecolor='white',linewidth=.3)
bx.axvline(17,color='#777777',linewidth=.7,linestyle=(0,(2,3)))
bx.axvline(32,color='#777777',linewidth=.7,linestyle=(0,(5,3)))
bx.text(.98,.93,'Full eligible wave 2, N = 8,902\nScore 10: 37.4%',transform=bx.transAxes,
        va='top',ha='right',fontsize=8.5)
bx.set(xlim=(9.4,45.5),ylim=(0,42),xticks=[10,17,20,30,32,40,45],yticks=[0,10,20,30,40],
       xlabel='Baseline CTQ neglect score (full observed distribution)',ylabel='Students (%)')
bx.set_title('b  Reference-population score distribution',loc='left',fontweight='bold',pad=8)

for ext in ['png','pdf','svg']:
    fig.savefig(OUT/f'FigureS4_function_form.{ext}',dpi=450,facecolor='white',bbox_inches='tight')
plt.close(fig)
shutil.copy2(SOURCE/'neglect_conditional_difference_curves.csv',OUT/'CTQ_FunctionForm_SourceData_20261015.csv')
audit={'source_files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [
    SOURCE/'neglect_conditional_difference_curves.csv',SOURCE/'neglect_full_W2_aggregate_distribution.csv']},
    'curve_rows':len(curves),'curves':4,'model_n':4219,'schools':10,'reference_population_n':8902,
    'range':[10,32],'reference_score':10,'bands':'pointwise school CR1 t(9) bands for all-linear and full-spline only',
    'interpretation':'conditional model contrasts, not a causal dose response or a validated threshold',
    'matplotlib_version':matplotlib.__version__,
    'outputs':{f'FigureS4_function_form.{ext}':hashlib.sha256((OUT/f'FigureS4_function_form.{ext}').read_bytes()).hexdigest() for ext in ['png','pdf','svg']}}
(OUT/'figureS4_source_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps({k:v for k,v in audit.items() if k not in ['outputs','source_files']},indent=2))
