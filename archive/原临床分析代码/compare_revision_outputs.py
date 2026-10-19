"""Compare re-fitted aggregate CSVs against frozen results, without modifying either."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

def compare_csv_set(actual, frozen, filenames=None, rtol=1e-10, atol=1e-10):
    actual=Path(actual);frozen=Path(frozen)
    paths=sorted(actual.glob('*.csv')) if filenames is None else [actual/f for f in filenames]
    if not paths:raise ValueError(f'No aggregate CSV files in {actual}')
    records=[]
    for p in paths:
        q=frozen/p.name
        if not q.is_file():raise FileNotFoundError(f'Frozen comparator missing: {q}')
        if p.read_bytes().strip()==b'' and q.read_bytes().strip()==b'':
            records.append(dict(file=p.name,rows=0,columns=0,max_numeric_absolute_difference=0.,passed=True));continue
        a=pd.read_csv(p);b=pd.read_csv(q)
        pd.testing.assert_frame_equal(a,b,check_dtype=False,rtol=rtol,atol=atol)
        numeric=[c for c in a.columns if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c])]
        maximum=0.
        for c in numeric:
            aa=a[c].to_numpy(float);bb=b[c].to_numpy(float)
            valid=np.isfinite(aa)&np.isfinite(bb)
            if valid.any():maximum=max(maximum,float(np.abs(aa[valid]-bb[valid]).max()))
        records.append(dict(file=p.name,rows=len(a),columns=len(a.columns),max_numeric_absolute_difference=maximum,passed=True))
    return dict(passed=True,rtol=rtol,atol=atol,n_files=len(records),max_numeric_absolute_difference=max(r['max_numeric_absolute_difference'] for r in records),files=records)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--actual',type=Path,required=True);p.add_argument('--frozen',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();report=compare_csv_set(a.actual,a.frozen)
    if a.report.exists():p.error('Report already exists; choose a new file.')
    a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
