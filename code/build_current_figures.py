"""Prepare current aggregate figure sources and render Figures 2 and 3.

No model fitting is performed. All plotted estimates and uncertainty limits are
read from existing result tables. Descriptive distributions are copied from
the existing aggregate source data. Requires Python, NumPy, pandas, matplotlib.
"""
from pathlib import Path
import argparse, csv, hashlib, io, json, shutil, sys, os
from types import MethodType
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parent/'matplotlib_cache'))
os.environ.setdefault('XDG_CACHE_HOME', str(Path(__file__).resolve().parent/'font_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
from matplotlib.offsetbox import DrawingArea
from matplotlib.artist import allow_rasterization

P = argparse.ArgumentParser(description=__doc__)
P.add_argument('--source-package', type=Path, required=True)
P.add_argument('--output', type=Path, required=True)
P.add_argument('--qa', type=Path, required=True)
P.add_argument('--audit-scripts', type=Path, required=True)
P.add_argument('--only-figure2', action='store_true',
               help='Prepare all aggregate source mappings, but render only Figure 2 into a staging directory')
A = P.parse_args()
sys.path.insert(0, str(A.audit_scripts))
from audit_panel_alignment import require_matplotlib_panel_alignment
SRC = A.source_package/'paper/evidence/source_data'
RES = A.source_package/'paper/evidence/results'
OUT, QA = A.output, A.qa
SD, FD = OUT/'source_data', OUT/'figures'
for folder in [SD, FD, QA, OUT/'code']:
    folder.mkdir(exist_ok=True, parents=True)

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
REG=[]
def source(relative, target, figure, panel, note='', predicate=None, preserve_csv_rows=False):
    p=A.source_package/relative
    d=pd.read_csv(p)
    original_n=len(d)
    if predicate is not None: d=d.loc[predicate(d)].copy()
    q=SD/target
    if predicate is None: shutil.copy2(p,q)
    elif preserve_csv_rows:
        # Keep the selected historical number strings and empty cells verbatim.
        # csv.reader.line_num also handles quoted fields containing newlines.
        raw=p.read_bytes().decode('utf-8')
        lines=raw.splitlines(keepends=True)
        reader=csv.reader(io.StringIO(raw,newline=''))
        next(reader)
        end=reader.line_num
        selected=[''.join(lines[:end])]
        selected_indices=set(d.index)
        for index,_ in enumerate(reader):
            start,end=end,reader.line_num
            if index in selected_indices:selected.append(''.join(lines[start:end]))
        q.write_bytes(''.join(selected).encode('utf-8'))
    else: d.to_csv(q,index=False)
    REG.append(dict(figure=figure,panel=panel,current_source_file=target,
                    source_package_relative_path=relative,source_sha256=sha(p),
                    current_sha256=sha(q),source_rows=original_n,current_rows=len(d),
                    operation='byte-identical copy' if predicate is None else 'explicit row selection; no numerical calculation',
                    selection_and_role=note))
    return d

def sd(name,target,figure,panel,note='',predicate=None,preserve_csv_rows=False):
    return source('paper/evidence/source_data/'+name,target,figure,panel,note,predicate,preserve_csv_rows)
def rr(rel,target,figure,panel,note=''):
    return source('paper/evidence/results/'+rel,target,figure,panel,note)

# Unchanged figures: exact source exports and stable source-to-panel mappings.
for name,panel in [('fig1a_neglect_distress_hexbin.csv','a'),('fig1a_distress_quantiles_by_neglect_band.csv','a'),
                   ('fig1a_neglect_marginal.csv','a'),('fig1a_distress_marginal.csv','a')]:
    sd(name,name,'Figure1',panel,'Unchanged original descriptive distribution; no fit')
sd('clinical_sensitivity.csv','fig1b_broad_linear_associations.csv','Figure1','b',
   'endpoint in distress/dep_score/anxiety_score/stress_score and specification prior_same_symptom',
   lambda d:d.endpoint.isin(['distress','dep_score','anxiety_score','stress_score']) & d.specification.eq('prior_same_symptom'))
sd('broad_reporting_estimates.csv','fig1c_unadjusted_adjusted_distress.csv','Figure1','c')
sd('manuscript_main_estimates.csv','fig1c_main_distress_units.csv','Figure1','c',
   'sample=All only; retired Discovery/Validation rows remain archived; original CSV values retained',
   lambda d:d['sample'].eq('All'),preserve_csv_rows=True)
for name,panel in [('figS1a_paired_distress_hexbin.csv','a'),('figS1b_distress_ecdf.csv','b'),
                   ('clinical_followup_descriptives.csv','c')]:sd(name,name,'FigureS1',panel)
sd('figS2_school_omission.csv','figS2_school_omission.csv','FigureS2','all')
sd('clinical_leave_one_school_out.csv','figS2_leave_one_school_out.csv','FigureS2','all')
sd('manuscript_main_estimates.csv','figS2_reference_main_estimates.csv','FigureS2','all',
   'sample=All only; retired Discovery/Validation rows remain archived; original CSV values retained',
   lambda d:d['sample'].eq('All'),preserve_csv_rows=True)
curves=sd('neglect_conditional_difference_curves.csv','figS3_conditional_difference_curves.csv','FigureS3','a',
   'Four existing specifications, 10–32 versus 10; pointwise CR1/t9 intervals')
freq=sd('neglect_full_W2_aggregate_distribution.csv','figS3_full_baseline_neglect_distribution.csv','FigureS3','b',
   'Full baseline N=8,902; not paired-sample support counts')

# Figure 2: retention, selection profile, and the two broad-distress intervals.
ret=sd('clinical_retention_by_grade.csv','fig2a_followup_by_grade.csv','Figure2','a')
smd=sd('fig2d_followup_profile_smd.csv','fig2b_followup_profile_smd.csv','Figure2','b')
temp=sd('source_clinical_temporal_extension.csv','fig2c_distress_intervals.csv','Figure2','c',
        'endpoint=distress; both exposures and both intervals retained',lambda d:d.endpoint.eq('distress'))

# Figure 3: unchanged descriptive content panels, with common-prior contrasts.
joint=sd('fig3a_content_joint.csv','fig3a_content_joint.csv','Figure3','a')
marg=sd('fig3a_content_marginals.csv','fig3a_content_marginals.csv','Figure3','a')
items=sd('fig3b_item_responses.csv','fig3b_item_responses.csv','Figure3','b')
full=sd('content_revision_common_baseline_contrasts.csv','fig3cd_full_common_prior_contrasts.csv','Figure3','c,d',
        'N=4,219; both units and both exposures; original within-unit two-exposure BH families')
grade=rr('revision_20261015/secondary_common_prior_contrasts.csv','fig3cd_secondary_common_prior_contrasts.csv','Figure3','c,d',
         'N=2,406; both units and exposures; original four-test BH family retained')
comp=sd('content_revision_common_baseline_slopes.csv','fig3e_full_common_prior_components.csv','Figure3','e',
        'Full-sample standardized components only; no restricted or own-prior component mixing',lambda d:d.units.eq('standardized'))
measure=sd('content_measurement_descriptives.csv','fig3f_content_measurement_descriptives.csv','Figure3','f')

assert joint.students.sum()==4219 and len(full)==len(grade)==4 and len(comp)==4
assert len(temp)==4 and ret['size'].sum()==8902 and ret['sum'].sum()==4219
assert np.allclose(items[[f'response_{i}_percent' for i in range(4)]].sum(axis=1),100)
assert set(grade['schools'])=={7} and set(full['schools'])=={10}
for ex in ['neglect','abuse']:
    c=comp.loc[comp.exposure.eq(ex)].set_index('content')
    diff=c.loc['anhedonic_content','beta']-c.loc['autonomic_content','beta']
    assert np.isclose(diff,float(full.loc[full.units.eq('standardized')&full.exposure.eq(ex),'beta'].iloc[0]),atol=1e-12)

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
    'font.size':7,'axes.labelsize':7,'axes.titlesize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,
    'legend.fontsize':6.2,'axes.linewidth':.7,'axes.spines.top':False,'axes.spines.right':False,
    'legend.frameon':False,'pdf.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':True,
    'xtick.major.size':2.5,'ytick.major.size':2.5,'axes.labelpad':2.5})
BLUE,RED,VIOLET,TEAL,GREY,BLACK='#0F4D92','#B64342','#9A4D8E','#42949E','#767676','#272727'
EX={'neglect':(BLUE,'o','Neglect'),'abuse':(RED,'s','Non-sexual abuse')}
CO={'anhedonic_content':(VIOLET,'o','Anhedonic content'),'autonomic_content':(TEAL,'^','Autonomic content')}
counts=LinearSegmentedColormap.from_list('content_counts',['#E3EAF4','#AFC8E4','#3775BA',BLUE,'#082A52'])

def axes(fig,x,y,w,h):
    W,H=fig.get_size_inches()*25.4
    return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def label(fig,x,y,letter,title):
    W,H=fig.get_size_inches()*25.4
    fig.text(x/W,1-y/H,letter,fontsize=8,fontweight='bold',ha='left',va='top')
    fig.text((x+4)/W,1-y/H-.0008,title,fontsize=7,ha='left',va='top')
def fmt(v,_=None):
    return ('0' if abs(v)<1e-10 else f'{v:.3f}'.rstrip('0').rstrip('.')).replace('-','−')
def pval(v):return '<0.001' if v<.001 else f'{v:.3f}'
def hit(df,**kw):
    mask=np.ones(len(df),bool)
    for k,v in kw.items():mask &= df[k].eq(v)
    rows=df.loc[mask]
    assert len(rows)==1,(kw,len(rows))
    return rows.iloc[0]
def point(ax,pos,r,style,vertical=False):
    color,marker,_=style
    b,l,h=[float(r[k]) for k in ['beta','school_low','school_high']]
    err=[[b-l],[h-b]]
    ax.errorbar(pos if vertical else b,b if vertical else pos,
                yerr=err if vertical else None,xerr=None if vertical else err,
                fmt=marker,ms=3.6,color=color,mec=color,mfc=color,mew=.7,
                elinewidth=.85,capsize=1.6,capthick=.7,zorder=3)
def forest(ax,labels,xlim,ticks,xlabel):
    ax.set(ylim=(len(labels)-.45,-.6),xlim=xlim,xticks=ticks,xlabel=xlabel)
    ax.set_yticks(range(len(labels)),labels)
    ax.tick_params(axis='y',length=0,pad=3)
    ax.spines['left'].set_visible(False)
    ax.axvline(0,color='#999999',lw=.6,zorder=0)
    ax.xaxis.set_major_formatter(FuncFormatter(fmt))
def handles(styles):
    return [Line2D([],[],ls='none',marker=m,ms=3.6,color=c,label=l) for c,m,l in styles.values()]
def draw_rasterized_legend_swatch(area,renderer):
    # OffsetBox computed this translation in the enclosing vector renderer's
    # display units. MixedModeRenderer changes DPI when rasterization begins;
    # carry the existing translation into that raster coordinate system.
    enclosing_scale=renderer.points_to_pixels(1.)
    @allow_rasterization
    def draw_in_raster_coordinates(artist,current_renderer):
        original=artist.offset_transform.get_matrix().copy()
        adjusted=original.copy()
        adjusted[:2,2]*=current_renderer.points_to_pixels(1.)/enclosing_scale
        artist.offset_transform.set_matrix(adjusted)
        try:DrawingArea.draw(artist,current_renderer)
        finally:artist.offset_transform.set_matrix(original)
    draw_in_raster_coordinates(area,renderer)
draw_rasterized_legend_swatch._supports_rasterization=True
def save(fig,name,panel_axes,panel_ids,row_groups=None,column_groups=None,exemptions=None):
    # Rasterize real data-bar hatches directly. A legend swatch needs its
    # DrawingArea rasterized as a unit so its DPI/offset transform is updated
    # before children draw. Rasterizing its child proxy alone displaces the
    # raster into the lower-left margin. Legend text remains vector/editable.
    for ax in fig.axes:
        for artist in ax.patches:
            if artist.get_hatch():artist.set_rasterized(True)
    for area in fig.findobj(DrawingArea):
        if any(isinstance(child,Patch) and child.get_hatch() for child in area.get_children()):
            area.draw=MethodType(draw_rasterized_legend_swatch,area)
            area.set_rasterized(True)
    fig.canvas.draw()
    alignment=require_matplotlib_panel_alignment(fig,axes=panel_axes,panel_ids=panel_ids,
        row_groups=row_groups or [],column_groups=column_groups or [],
        exemptions=exemptions or [],tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True,
        json_out=QA/(name+'.alignment.json'),overlay_svg=QA/(name+'.alignment.svg'))
    fig.savefig(FD/(name+'.png'),dpi=600)
    fig.savefig(FD/(name+'.pdf'),dpi=600)
    fig.savefig(FD/(name+'.svg'),dpi=600)
    fig.savefig(QA/(name+'_preview.png'),dpi=300)
    plt.close(fig)
    return alignment['verdict']

# Figure 2 (180 x 136 mm), equal plot areas for the upper diagnostic panels.
width_mm=180.0
fig=plt.figure(figsize=(width_mm/25.4,136/25.4))
label(fig,2,2,'a','Matched wave-3 records by grade')
ax_a=axes(fig,30,14,52,66)
for i,r in enumerate(ret.itertuples()):
    pct=float(r.mean)*100
    ax_a.barh(i,pct,height=.65,color='#8BCF8B',edgecolor=BLACK,lw=.45)
    ax_a.barh(i,100-pct,left=pct,height=.65,color='#CFCECE' if pct else '#FFFFFF',
              edgecolor=BLACK,lw=.45,hatch=None if pct else '////')
ax_a.set(ylim=(8.5,-.6),xlim=(0,100),xticks=[0,50,100],xlabel='Students (%)')
ax_a.set_yticks(range(9),[f'{int(r.grade)}   {int(r.sum):,}/{int(r.size):,}' for r in ret.itertuples()])
fig.text(4/180,1-11/136,'Grade  Matched / n',fontsize=6,ha='left',va='bottom')
ax_a.tick_params(axis='y',length=0)
ax_a.spines['left'].set_visible(False)
fig.legend(handles=[Patch(facecolor='#8BCF8B',label='Matched record'),Patch(facecolor='#CFCECE',label='No matched record'),
                    Patch(facecolor='white',edgecolor=BLACK,hatch='////',label='Grades 6, 9, 12')],
           loc='upper left',bbox_to_anchor=(.13,.365),ncol=1,fontsize=6)

label(fig,92,2,'b','Baseline matching profile')
ax_b=axes(fig,124,14,52,66)
order=smd.loc[smd.comparison.eq('All wave-2 students'),'variable'].tolist()
for i,v in enumerate(order):
    for comp_,dy,fill in [('All wave-2 students',-.15,True),('Continuing grades',.15,False)]:
        rows=smd.loc[smd.variable.eq(v)&smd.comparison.eq(comp_)]
        assert len(rows)==1,(v,comp_,smd.comparison.unique())
        val=float(rows.iloc[0].smd)
        ax_b.plot(val,i+dy,'o',ms=3.2,mec=BLUE,mfc=BLUE if fill else 'white',mew=.7)
ax_b.set_yticks(range(len(order)),order)
ax_b.set(ylim=(len(order)-.5,-.6),xlim=(-1.65,.15),xticks=[-1.5,-1,-.5,0],xlabel='Standardized mean difference')
assert smd.smd.min()>ax_b.get_xlim()[0] and smd.smd.max()<ax_b.get_xlim()[1]
ax_b.xaxis.set_major_formatter(FuncFormatter(fmt))
ax_b.tick_params(axis='y',length=0)
ax_b.spines['left'].set_visible(False)
for v in [-.1,.1]:ax_b.axvline(v,color='#AAAAAA',ls=':',lw=.7,zorder=0)
ax_b.axvline(0,color='#999999',lw=.6,zorder=0)
ax_b.legend(handles=[Line2D([],[],ls='none',marker='o',ms=3.2,color=BLUE,label='All baseline students'),
                      Line2D([],[],ls='none',marker='o',ms=3.2,mec=BLUE,mfc='white',label='Continuing grades')],
            loc='upper left',bbox_to_anchor=(-.02,-.13),fontsize=6)

label(fig,2,102,'c','General distress across cohort intervals')
ax_c=axes(fig,53,111,88,15)
forest(ax_c,['Waves 1–2\n2,022; 4 schools','Waves 2–3\n4,219; 10 schools'],(-.02,.27),[0,.1,.2],
       'Standardized coefficient (95% school CI)')
ax_c.text(1.055,1.04,'Wild p',transform=ax_c.transAxes,fontsize=6,va='bottom')
for i,t in enumerate(['1-2','2-3']):
    for ex,dy in [('neglect',-.15),('abuse',.15)]:
        r=hit(temp,transition=t,endpoint='distress',exposure=ex)
        point(ax_c,i+dy,r,EX[ex])
        ax_c.text(1.055,i+dy,pval(r.p_wild),transform=ax_c.get_yaxis_transform(),color=EX[ex][0],fontsize=6,va='center')
fig.legend(handles=handles(EX),loc='upper right',bbox_to_anchor=(.985,.245),ncol=1,fontsize=6)
align2=save(fig,'Figure2_followup_and_intervals',[ax_a,ax_b,ax_c],{ax_a:'a',ax_b:'b',ax_c:'c'},
    row_groups=[['a','b']],exemptions=[{'panels':['c'],'checks':['row','column','panel-width','horizontal-gutter','vertical-gutter'],
    'reason':'Temporal forest is a separate evidence role below the two baseline diagnostics'}])

if A.only_figure2:
    pd.DataFrame(REG).to_csv(SD/'figure_source_registry.csv',index=False)
    (QA/'figure_generation_summary.json').write_text(json.dumps(dict(
        model_fitting=False,source_rows=REG,new_figures=['Figure2_followup_and_intervals'],
        fix='Legend hatch DrawingArea rasterized with its own DPI/offset transform; text remains vector',
        alignment={'Figure2':align2}),indent=2)+'\n')
    print(json.dumps(dict(rendered=['Figure2_followup_and_intervals'],models_fitted=0,
        source_tables=len(REG),alignment={'Figure2':align2}),indent=2))
    sys.exit(0)

# Figure 3 (180 x 180 mm), source distributions and paired direct comparisons.
fig=plt.figure(figsize=(width_mm/25.4,180/25.4))
label(fig,2,2,'a','Joint follow-up content (r = 0.748)')
ax=axes(fig,14,18,47,47);top=axes(fig,14,9.5,47,7);side=axes(fig,62.5,18,9,47)
mat=np.zeros((10,13))
for r in joint.itertuples():mat[round(r.anhedonic_average*3),round(r.autonomic_average*4)]=r.percent
im=ax.imshow(np.ma.masked_where(mat==0,mat),origin='lower',extent=(-.125,3.125,-1/6,3+1/6),
             aspect='auto',cmap=counts,norm=LogNorm(vmin=.02,vmax=60),interpolation='none')
ax.set(xlabel='Autonomic item average (0–3)',ylabel='Anhedonic item average (0–3)',xticks=[0,1,2,3],yticks=[0,1,2,3])
ax.text(.05,-.14,f'{mat[0,0]:.1f}%',color='white',fontsize=6,ha='left',va='bottom')
for content,a,orient,col in [('autonomic',top,'v',TEAL),('anhedonic',side,'h',VIOLET)]:
    d=marg.loc[marg.content.eq(content)].sort_values('item_average')
    if orient=='v':
        a.bar(d.item_average,d.percent,width=.17,color=col);a.set(xlim=ax.get_xlim(),ylim=(0,75),yticks=[0,60],xticks=[])
        a.set_ylabel('%',rotation=0,rotation_mode='anchor',ha='right',va='center',labelpad=1)
        a.text(.12,61,f'{float(d.iloc[0].percent):.1f}% zero',color=col,fontsize=6,va='center')
    else:
        a.barh(d.item_average,d.percent,height=.23,color=col);a.set(ylim=ax.get_ylim(),xlim=(0,75),xticks=[0,60],yticks=[])
        a.set_xlabel('%',labelpad=1);a.text(5,2.25,f'{float(d.iloc[0].percent):.1f}%\nzero',color=col,fontsize=6)
cax=axes(fig,76,34,1.8,30);cb=fig.colorbar(im,cax=cax);cb.set_ticks([.1,1,10]);cb.set_ticklabels(['0.1','1','10'])
cb.ax.minorticks_off();cb.set_label('Students (%)',labelpad=1)

label(fig,92,2,'b','Item responses at follow-up')
ax_b3=axes(fig,122,15,50,45)
names={3:'Positive feeling',10:'Anticipation',16:'Enthusiasm',2:'Dry mouth',4:'Breathing difficulty',7:'Trembling',19:'Heart awareness'}
ypos=[0,1,2,3.7,4.7,5.7,6.7]
ramps={'anhedonic':['#F2F2F2','#D9B8D4','#B77AAE',VIOLET],'autonomic':['#F2F2F2','#B5DADE','#7AB8BF',TEAL]}
for y,r in zip(ypos,items.itertuples()):
    left=0
    for k in range(4):
        v=getattr(r,f'response_{k}_percent');ax_b3.barh(y,v,left=left,height=.7,color=ramps[r.content][k],edgecolor=BLACK,lw=.35);left+=v
    ax_b3.text(1.5,y,f'{r.response_0_percent:.1f}%',va='center',fontsize=6)
    ax_b3.text(101.5,y,f'{100-r.response_0_percent:.1f}%',va='center',fontsize=6,color=VIOLET if r.content=='anhedonic' else TEAL)
ax_b3.set_yticks(ypos,[f'{names[int(i)]} ({int(i)})' for i in items.item])
ax_b3.set(xlim=(0,100),ylim=(7.25,-1.25),xticks=[0,50,100],xlabel='Students (%)')
ax_b3.tick_params(axis='y',length=0);ax_b3.spines['left'].set_visible(False)
ax_b3.text(1.5,-.85,'0',fontsize=6,va='center');ax_b3.text(101.5,-.85,'≥1',fontsize=6,va='center')
for tick,r in zip(ax_b3.get_yticklabels(),items.itertuples()):tick.set_color(VIOLET if r.content=='anhedonic' else TEAL)
fig.text(94/180,1-72.5/180,'Response:',fontsize=6,va='center')
for g,(group,x0) in enumerate([('anhedonic',114),('autonomic',146)]):
    for k in range(4):
        fig.patches.append(plt.Rectangle(((x0+k*6)/180,1-74.5/180),5/180,4/180,transform=fig.transFigure,
                                         facecolor=ramps[group][k],edgecolor='none'))
        fig.text((x0+k*6+2.5)/180,1-72.5/180,str(k),ha='center',va='center',fontsize=5.8,
                 color='white' if k==3 else BLACK)

label(fig,2,79,'c','Standardized content contrast')
label(fig,92,79,'d','Common response-unit contrast')
ax_c3=axes(fig,31,90,40,29);ax_d3=axes(fig,123,90,40,29)
labs=['All grades\n4,219; 10 schools','Grades ≥7\n2,406; 7 schools']
for a,unit,xlim,ticks,xlab in [(ax_c3,'standardized',(-.13,.13),[-.1,0,.1],'Anhedonic − autonomic\n(baseline outcome SD)'),
                              (ax_d3,'raw_0to3',(-.065,.075),[-.04,0,.04],'Anhedonic − autonomic\n(item-average points)')]:
    forest(a,labs,xlim,ticks,xlab)
    a.text(1.06,1.02,'Wild p',transform=a.transAxes,fontsize=6,va='bottom')
    for i,d in enumerate([full,grade]):
        for ex,dy in [('neglect',-.16),('abuse',.16)]:
            r=hit(d,units=unit,exposure=ex);point(a,i+dy,r,EX[ex])
            a.text(1.06,i+dy,pval(r.p_wild),transform=a.get_yaxis_transform(),color=EX[ex][0],fontsize=6,va='center')
fig.legend(handles=handles(EX),loc='center',bbox_to_anchor=(.5,1-132/180),ncol=2,fontsize=6.3)

label(fig,2,140,'e','Full-sample common-prior component slopes')
ax_e=axes(fig,15,150,62,23)
for i,ex in enumerate(['neglect','abuse']):
    for content,dx in [('anhedonic_content',-.12),('autonomic_content',.12)]:
        r=hit(comp,exposure=ex,content=content);point(ax_e,i+dx,r,CO[content],vertical=True)
ax_e.set(xlim=(-.5,1.5),ylim=(-.01,.31),yticks=[0,.1,.2,.3],ylabel='β (95% school CI)')
ax_e.yaxis.set_major_formatter(FuncFormatter(fmt));ax_e.set_xticks([0,1],['Neglect','Non-sexual abuse'])
ax_e.axhline(0,color='#999999',lw=.6,zorder=0)
fig.legend(handles=handles(CO),loc='upper left',bbox_to_anchor=(.01,.245),fontsize=6,ncol=2)

label(fig,100,140,'f','Content measurement properties')
visits=[('1-2',1),('1-2',2),('2-3',2),('2-3',3)];xv=np.array([0,1,2.4,3.4]);f_axes=[]
for k,(col,title,ylim,ticks) in enumerate([('zero_percent','Zero score (%)',(0,85),[0,40,80]),
                    ('alpha_raw',"Cronbach’s α",(0,1),[0,.5,1]),('between_content_r','Content r',(0,1),[0,.5,1])]):
    a=axes(fig,112+k*23,150,17,23);f_axes.append(a)
    for content,dx in [('anhedonic_content',-.19),('autonomic_content',.19)]:
        vals=[float(hit(measure,transition=t,visit=v,content=content)[col]) for t,v in visits]
        if col=='between_content_r':
            if content=='anhedonic_content':a.bar(xv,vals,width=.62,color='#CFCECE',edgecolor=BLACK,lw=.4)
        else:a.bar(xv+dx,vals,width=.36,color=VIOLET if content=='anhedonic_content' else 'white',
                   edgecolor=BLACK if content=='anhedonic_content' else TEAL,lw=.6)
    a.set(xlim=(-.6,4),ylim=ylim,yticks=ticks);a.set_title(title,fontsize=6.2,pad=3)
    a.yaxis.set_major_formatter(FuncFormatter(fmt));a.set_xticks(xv,['W1','W2','W2','W3']);a.tick_params(axis='x',length=0)
    for x0,n in [(.5,'2,022'),(2.9,'4,219')]:a.text(x0,-.18,n,transform=a.get_xaxis_transform(),ha='center',va='top',fontsize=5.8,color=GREY)
align3=save(fig,'Figure3_common_prior_content',[ax_c3,ax_d3,*f_axes],
    {ax_c3:'c',ax_d3:'d',**{a:f'f{i+1}' for i,a in enumerate(f_axes)}},
    row_groups=[['c','d'],['f1','f2','f3']])

# Recreate the existing function-form plot without changing numerical geometry.
# Only the Q75 guide is shortened below the legend to avoid crossing its text.
plt.rcdefaults()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.5,'axes.labelsize':10,
                     'axes.titlesize':10,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none',
                     'axes.spines.top':False,'axes.spines.right':False})
fig,(sx,tx)=plt.subplots(2,1,figsize=(7.2,5.7),gridspec_kw={'height_ratios':[2.25,1]},
                         constrained_layout=True)
specs=[('all_linear','Linear exposures and baseline','#0072B2','-',1.7),
       ('exposure_linear_baseline_spline','Linear exposures; spline baseline','#009E73','--',1.2),
       ('exposure_splines_baseline_linear','Spline exposures; linear baseline','#CC79A7',':',1.7),
       ('exposure_and_baseline_splines','Spline exposures and baseline','#D55E00','-',1.9)]
assert len(curves)==92 and freq.n.sum()==8902
for key,title,color,linestyle,lw in specs:
    d=curves.loc[curves.specification.eq(key)].sort_values('neglect_score')
    assert len(d)==23 and d.neglect_score.min()==10 and d.neglect_score.max()==32
    if key in ['all_linear','exposure_and_baseline_splines']:
        sx.fill_between(d.neglect_score,d.school_low,d.school_high,color=color,alpha=.11,linewidth=0)
    sx.plot(d.neglect_score,d.estimate,color=color,linestyle=linestyle,lw=lw,label=title)
sx.axhline(0,color='#666666',lw=.7,zorder=0)
sx.axvline(17,ymin=0,ymax=.65,color='#888888',lw=.7,linestyle=(0,(2,3)),zorder=0)
sx.text(17.15,.55,'Q75 = 17',transform=sx.get_xaxis_transform(),va='center',ha='left',fontsize=8,color='#666666')
sx.set(xlim=(10,32),xticks=[10,15,20,25,30,32],xlabel='Baseline CTQ neglect score (curve range 10 to P95 = 32)',
       ylabel='Conditional general-distress difference\nfrom neglect score 10 (baseline SD)')
sx.legend(loc='upper left',frameon=False,borderaxespad=.5,handlelength=2.8)
sx.set_title('a  Model-dependent conditional differences',loc='left',fontweight='bold',pad=8)
tx.bar(freq.neglect_score,freq.percent,width=.88,color='#71889B',edgecolor='white',lw=.3)
tx.axvline(17,color='#777777',lw=.7,linestyle=(0,(2,3)))
tx.axvline(32,color='#777777',lw=.7,linestyle=(0,(5,3)))
tx.text(.98,.93,'Full eligible wave 2, N = 8,902\nScore 10: 37.4%',transform=tx.transAxes,va='top',ha='right',fontsize=8.5)
tx.set(xlim=(9.4,45.5),ylim=(0,42),xticks=[10,17,20,30,32,40,45],yticks=[0,10,20,30,40],
       xlabel='Baseline CTQ neglect score (full observed distribution)',ylabel='Students (%)')
tx.set_title('b  Reference-population score distribution',loc='left',fontweight='bold',pad=8)
fig.canvas.draw()
alignS3=require_matplotlib_panel_alignment(fig,axes=[sx,tx],panel_ids={sx:'a',tx:'b'},column_groups=[['a','b']],
       tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True,json_out=QA/'FigureS3_function_form.alignment.json',
       overlay_svg=QA/'FigureS3_function_form.alignment.svg')['verdict']
fig.savefig(FD/'FigureS3_function_form.png',dpi=600,facecolor='white',bbox_inches='tight')
fig.savefig(FD/'FigureS3_function_form.pdf',dpi=600,facecolor='white',bbox_inches='tight')
fig.savefig(FD/'FigureS3_function_form.svg',dpi=600,facecolor='white',bbox_inches='tight')
plt.close(fig)

copies=[]
for old,new in [('Figure1_broad_symptoms','Figure1_broad_symptoms'),('FigureS1_retention','FigureS1_paired_distress'),
                ('FigureS2_school_omission','FigureS2_school_omission')]:
    for ext in ['png','pdf','svg']:
        p=A.source_package/'paper/figures'/f'{old}.{ext}';q=FD/f'{new}.{ext}';shutil.copy2(p,q)
        assert sha(p)==sha(q)
        copies.append(dict(old_name=p.name,current_name=q.name,sha256=sha(q),byte_identical=True))
pd.DataFrame(REG).to_csv(SD/'figure_source_registry.csv',index=False)
(SD/'README_figure_sources.md').write_text('''# Current figure source data\n\nSix figures are active. Figures 1, S1 and S2 are byte-identical exports. Figures 2 and 3 were redrawn from existing aggregate results without fitting models. Renumbered S3 uses the same 92 curve rows, frequency distribution, four lines and two confidence bands; its Q75 guide was shortened below the legend to avoid crossing text.\n\nThe registry maps every active panel to its current source files, original relative source location, selection and SHA-256. The main Figure 3 contrasts and component slopes use common three-domain baseline adjustment; its component panel shows the full sample only. Earlier own-prior content results, adult comparisons, hash subsets and the removed mixed-specification heatmap are not active figure sources. Their historical results are retained separately.\n\nAll intervals retain the source definition. School-based model panels use CR1/t(G−1), with G=10, 7 or 4 as labeled. The function-form curve intervals are pointwise. Descriptive frequencies, alpha and correlations are summaries without inferential error bars. Full-baseline distributions are not paired-sample support counts.\n''')
with (SD/'README_figure_sources.md').open('a') as source_readme:
    source_readme.write('''\nThe Figure 1c unit and Figure S2 reference sources retain only sample=All (two of six original rows). Selected CSV source lines preserve their original numerical strings and empty cells; complete Discovery/Validation rows remain archived. Registry row counts, selection rules and current hashes reflect this selection.\n\nFigure 2 PDF/SVG exports rasterize the legend hatch DrawingArea with its offset converted into the raster DPI coordinate system, keeping the hatch at its intended legend location and all text editable. This fixes the displaced lower-left proxy mark without changing the 600-dpi PNG or any statistical geometry.\n''')
(QA/'figure_generation_summary.json').write_text(json.dumps(dict(model_fitting=False,
    source_rows=REG,byte_identical_copies=copies,new_figures=['Figure2_followup_and_intervals','Figure3_common_prior_content'],
    format_only_redraw='FigureS3_function_form: Q75 guide shortened; curve data unchanged',
    alignment={'Figure2':align2,'Figure3':align3,'FigureS3':alignS3},component_contrast_identity_checked=True),indent=2)+'\n')
print(json.dumps(dict(active_figure_count=6,exports=18,newly_drawn=2,format_only_redrawn=1,byte_identical_copied=3,
    models_fitted=0,source_tables=len(REG),alignment={'Figure2':align2,'Figure3':align3,'FigureS3':alignS3}),indent=2))
