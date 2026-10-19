"""Export descriptive score support only; this script fits no statistical models."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, type=Path)
parser.add_argument('--output', required=True, type=Path)
args = parser.parse_args()
data = pd.read_csv(args.input)
baseline = data.loc[data.wave.eq(2)].copy()
followup = data.loc[data.wave.eq(3)]
paired = baseline.merge(followup[['participant_code']], on='participant_code', validate='1:1')
assert len(baseline) == 8902 and len(paired) == 4219
paired['neglect_score'] = paired.ctq_affective_ignore_score + paired.ctq_body_ignore_score
assert paired.school_code.nunique() == 10
args.output.mkdir(parents=True, exist_ok=True)
support = []
for score in [10, 17, 25, 32]:
    for kind in ['exact', 'at_or_above']:
        mask = paired.neglect_score.eq(score) if kind == 'exact' else paired.neglect_score.ge(score)
        support.append({'reference_population':'paired_W2_W3','n_population':4219,
                        'baseline_score':score,'support_type':kind,'n_participants':int(mask.sum()),
                        'percent':100*float(mask.mean()),'schools_with_support':int(paired.loc[mask,'school_code'].nunique())})
pd.DataFrame(support).to_csv(args.output/'high_score_paired_support.csv',index=False)
distribution = paired.groupby(['school_code','neglect_score']).size().rename('n').reset_index()
distribution.to_csv(args.output/'paired_neglect_score_support_by_school.csv',index=False)
(args.output/'high_score_support_metadata.json').write_text(json.dumps({
    'input_sha256':sha256(args.input.read_bytes()).hexdigest(),
    'baseline_n':8902,'paired_n':4219,'schools':10,
    'purpose':'Descriptive support for already computed score contrasts; not a new regression or test.',
    'individual_identifiers_exported':False,'model_fits':0},ensure_ascii=False,indent=2))
print(pd.DataFrame(support).to_string(index=False))
