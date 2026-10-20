"""Independent aggregate-output checks for the locked CTQ specificity run.

No participant-level values are loaded. The private input is hashed in a stream.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EXPECTED = {
    'plan': '645da2b59fdab730fa604924cbf908b80066e5dfeda0440e2e3d733e0f208cda',
    'private_input': '5342bd167496989d0ebf083894d24547c50b210356715cdeceb49b76b5ad352b',
    'original_main_word': '6a7c61574f23fda96e1e4279560aceebe8e95a14c9ebfc44aeb4524214ca01cf',
    'original_supplement_word': '9a5312b27ee2a12e76d1b7765fda2674b8c97cd007ca4f30489f041b1f0be31f',
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def close(a: str | float, b: str | float, tol: float = 1e-9) -> bool:
    return abs(float(a) - float(b)) <= tol


def main() -> None:
    parser = argparse.ArgumentParser(description='Check aggregate outputs and frozen source identities')
    parser.add_argument('--run', type=Path, required=True,
                        help='Local directory with fractional/, ancillary/, and measurement/ outputs')
    parser.add_argument('--data', type=Path, required=True, help='Authorized participant table')
    parser.add_argument('--frozen', type=Path, required=True, help='Frozen Table S14 source table')
    parser.add_argument('--original-plan', type=Path)
    parser.add_argument('--original-main-word', type=Path)
    parser.add_argument('--original-supplement-word', type=Path)
    opts = parser.parse_args()
    run = opts.run.resolve()
    if run.is_relative_to(REPO):
        raise ValueError('Keep aggregate outputs outside the public code repository')
    source_checks = {'private_input': opts.data}
    for name, path in [('plan', opts.original_plan),
                       ('original_main_word', opts.original_main_word),
                       ('original_supplement_word', opts.original_supplement_word)]:
        if path is not None:
            source_checks[name] = path
    assert all(digest(path) == EXPECTED[name] for name, path in source_checks.items())
    grid = rows(run / 'fractional/four_model_grid.csv')
    frozen = rows(opts.frozen)
    assert len(grid) == 4
    assert {(x['model'], x['form']) for x in grid} == {
        (model, form) for model in ('OLS', 'fractional_logit')
        for form in ('linear', 'natural_cubic_spline_3df')}
    frozen_by_form = {x['function_form']: x for x in frozen if x['system'] == 'content'
                      and x['outcome_units'] == 'raw_0to3' and x['exposure'] == 'neglect'}
    for x in grid:
        assert int(x['n_participants']) == 4219 and int(x['n_rows']) == 8438
        assert int(x['schools']) == 10 and int(x['school_df']) == 9
        assert int(x['parameters']) == (48 if x['form'] == 'linear' else 68)
        assert close(float(x['delta']), float(x['a_difference']) - float(x['u_difference']))
        assert float(x['school_low']) < float(x['delta']) < float(x['school_high'])
        if x['model'] == 'fractional_logit':
            for col in ('predicted_at_10_A', 'predicted_at_17_A',
                        'predicted_at_10_U', 'predicted_at_17_U'):
                assert 0 < float(x[col]) < 3
        else:
            old = frozen_by_form[x['form']]
            for actual, reference in [('delta', 'estimate'), ('se_school', 'se_school'),
                                      ('school_low', 'school_low'), ('school_high', 'school_high')]:
                assert close(x[actual], old[reference])
    loso = rows(run / 'fractional/fractional_leave_one_school_out.csv')
    assert len(loso) == 20
    assert all(int(x['remaining_schools']) == 9 and x['converged'] == 'True' for x in loso)
    assert {x['form'] for x in loso} == {'linear', 'natural_cubic_spline_3df'}
    ancillary = rows(run / 'ancillary/ancillary_contrasts.csv')
    by_name = {x['check']: x for x in ancillary}
    assert len(ancillary) == len(by_name) == 4
    assert int(by_name['wave1_to_2_common_prior_linear_OLS']['n']) == 2022
    assert int(by_name['wave1_to_2_common_prior_linear_OLS']['schools']) == 4
    assert int(by_name['grade_ge7_existing_common_prior_linear_OLS_rescaled']['n']) == 2406
    assert int(by_name['wave2_to_3_IPW_p99_common_prior_linear_WLS']['n']) == 4219
    assert close(by_name['wave2_to_3_unweighted_identity_check']['delta_anhedonic_minus_autonomic'],
                 frozen_by_form['linear']['estimate'])
    freq = rows(run / 'measurement/item_category_frequencies.csv')
    assert len(freq) == 6 * 7 * 4
    counts: dict[tuple[str, str, str], int] = defaultdict(int)
    denominators: dict[tuple[str, str, str], int] = {}
    for x in freq:
        key = (x['wave'], x['stage'], x['item'])
        counts[key] += int(x['count'])
        assert x['category'] in {'0', '1', '2', '3'}
        denominators[key] = int(x['denominator'])
    assert len(counts) == 6 * 7
    assert all(counts[k] == denominators[k] for k in counts)
    factor = rows(run / 'measurement/ordinal_two_factor_diagnostics.csv')
    assert len(factor) == 6 and all(x['status'] == 'estimated' for x in factor)
    assert all(0 <= float(x['factor_r']) <= 1 for x in factor)
    assert len(rows(run / 'measurement/ordinal_two_factor_loadings.csv')) == 42
    assert len(rows(run / 'measurement/ordinal_thresholds.csv')) == 126
    assert len(rows(run / 'measurement/polychoric_pair_diagnostics.csv')) == 126
    assert len(rows(run / 'measurement/content_descriptives.csv')) == 12
    for relative in ('fractional/diagnostics.json', 'ancillary/ancillary_audit.json',
                     'measurement/audit.json'):
        assert json.loads((run / relative).read_text())['source_sha256'] == digest(opts.data)
    for output in run.rglob('*.csv'):
        assert output.name != opts.data.name
        with output.open(encoding='utf-8-sig') as stream:
            header = stream.readline().lower()
            assert 'participant_code' not in header and 'patient_id' not in header
    report = {
        'status': 'pass',
        'locked_source_sha256': {name: digest(path) for name, path in source_checks.items()},
        'four_grid_models': len(grid), 'loso_fits': len(loso),
        'ancillary_checks': len(ancillary), 'measurement_item_frequency_rows': len(freq),
        'data_privacy': 'No participant-level table copied to run; aggregate CSV headers screened.',
        'scope': 'Checks internal aggregate consistency and original-file integrity; does not independently validate clinical constructs or causal interpretation.',
    }
    output = run / 'validation'
    output.mkdir(exist_ok=True)
    (output / 'verify_run.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(report, ensure_ascii=False))


if __name__ == '__main__':
    main()
