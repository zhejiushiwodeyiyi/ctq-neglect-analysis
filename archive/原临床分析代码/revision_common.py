"""Numerical helpers for canonical 0-3 DASS revision analyses.

Functions are unchanged from the historical common.py; private raw-data loaders
are deliberately not exposed by this active package module.
"""
import numpy as np
from scipy import linalg
from statsmodels.stats.multitest import multipletests
SUB=['ctq_affective_ignore_score','ctq_body_ignore_score','ctq_affective_abuse_score','ctq_body_abuse_score']

def sex_code(series):
    """Normalize explicit English/Chinese labels without imputing unknown sex."""
    normalized=series.astype('string').str.strip().str.lower()
    return normalized.map({'female':1.,'male':0.,'女':1.,'男':0.})

def bh(p):
    p=np.asarray(p,float); q=np.full(p.shape,np.nan); ok=np.isfinite(p)
    if ok.any(): q[ok]=multipletests(p[ok],method='fdr_bh')[1]
    return q

def full_rank(X):
    """Remove redundant nuisance columns while preserving the exposure estimand."""
    others=[c for c in X.columns if c!='ctq_z']
    a=X[others].to_numpy(float)
    _,r,piv=linalg.qr(a,mode='economic',pivoting=True)
    tol=max(a.shape)*np.finfo(float).eps*abs(r).max()
    rank=int((abs(np.diag(r))>tol).sum())
    keep=[others[j] for j in sorted(piv[:rank])]
    residual=X.ctq_z.to_numpy()-a[:,piv[:rank]]@np.linalg.lstsq(a[:,piv[:rank]],X.ctq_z,rcond=None)[0]
    if np.linalg.norm(residual)<1e-8: raise ValueError('Exposure is aliased with nuisance covariates')
    return X[keep+['ctq_z']]
