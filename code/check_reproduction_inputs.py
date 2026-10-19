"""Read-only input/path/fingerprint preflight for this CTQ package.

Uses Python standard libraries only. Never imports analysis scripts, fits models,
enumerates bootstrap patterns, installs packages, or modifies package files.
"""
from pathlib import Path
import argparse
import ast
import csv
import hashlib
import importlib.metadata
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
# Frozen against the independently checked 16 October package; historical files
# and participant inputs are unchanged in the 17 October cleanup.
EXPECTED_SHA256 = {
    "archive/原临床分析代码/compare_revision_outputs.py": "bfe993eb3f5338ea9471820b985bf95d97e5b98b05164c388192b5403f5cc620",
    "archive/原临床分析代码/post_review_20261015/POST_REVIEW_DIAGNOSTIC_ADDENDUM_20261015.md": "4dfe10e58b24ea12956fa68828b491af77cb06b30944c512ba431065ce2b776a",
    "archive/原临床分析代码/post_review_20261015/POST_REVIEW_PLAN_20261015.md": "f1225b922caf9581d77ece1153d963511fc4800fe1dd3ee9b89c60a72ae322f4",
    "archive/原临床分析代码/post_review_20261015/adult_scaling_register_note.md": "6eaf69846e8117a6add6fe555dcb623677d72a150c176253b10dca61e065c44e",
    "archive/原临床分析代码/post_review_20261015/build_diagnostic_figure.py": "8f7dc4e35ed22a84b779d1301dadb2a14939f702f274db9ea05681a8c7d34cd5",
    "archive/原临床分析代码/post_review_20261015/dependencies/content_analysis_20261012.py": "eb0de7848322dfe66f89e652306f33cd155340515cf68988db1c36afd1bf5c07",
    "archive/原临床分析代码/post_review_20261015/diagnostic_addendum.py": "f6407bc3b22682454f49ab71649ebcfa5d2ce133598f9a3c8991305197f9b370",
    "archive/原临床分析代码/post_review_20261015/post_review_analysis.py": "2dcb31efe67e51b330ba0ab707f6e368519cca3c2cd4912b2196f4fa18cc2006",
    "archive/原临床分析代码/post_review_20261015/reference_results/content_raw_unit_slopes.csv": "44ad598ad06440022e045b7eda6b80128bfa059322e6224462eb41be018109fd",
    "archive/原临床分析代码/post_review_20261015/reference_results/content_revision_common_baseline_contrasts.csv": "d1549ca58cd3c087f4ea0436cb62ee9a37e0e96b95f6a89f1322db6c5612aa49",
    "archive/原临床分析代码/post_review_20261015/reference_results/external_clinical_sample.csv": "de089543fd44069b389d82d455e3f1fdc01f371c1dfd790b025865aabf812829",
    "archive/原临床分析代码/post_review_20261015/reference_results/external_clinical_scaling.json": "9157d32860eeb0e4fc0a4715f8127072929ad3bf2839120909e91392dbac2e1c",
    "archive/原临床分析代码/post_review_20261015/reference_results/manuscript_scaling.json": "65ef036efa19c9a519d882ea1c0c1821efec1414d5ed39ad61747634ef9e48ee",
    "archive/原临床分析代码/post_review_20261015/reference_results/source_clinical_content_contrast.csv": "ce8cf5567289ca083d5fc409c129f39d2db4369060f117cdd86903a8b7e4a748",
    "archive/原临床分析代码/post_review_20261015/scaling_register_reporting_corrected.csv": "845c0ce917af6372d448dfe86f9a020e0dc754812167c1c08b2465238a999ce0",
    "archive/原临床分析代码/revision_20261012/CONTENT_PLAN.md": "e50407135ee684b335c6521a074deefd9611e113c5a2a1ef1236f622d3839c35",
    "archive/原临床分析代码/revision_20261012/IPW_PLAN.md": "32537fd3ba2a1f298d1157c06433f63255e95ecc816364369635e224bb6d1cd9",
    "archive/原临床分析代码/revision_20261012/content_analysis.py": "b5edfd62bf5c72d01741d4ed4ba23cfe467156b2353cafc349bca8756571e4b0",
    "archive/原临床分析代码/revision_20261012/ipw_analysis.py": "4feb19046782a7a5694f5c93879768f1410dac775adef9de984407f7665a47e8",
    "archive/原临床分析代码/revision_20261012/reference_results/clinical_sensitivity.csv": "ee3141a10d5c50bbe99f92ccbed03cde70bc465c8e89efa1378c0ab2960e0502",
    "archive/原临床分析代码/revision_20261012/reference_results/clinical_structure_summary.json": "61ecb8df70274888c89e5aa8be778a7fc2135cb2c5096a85ff55c37d95123d51",
    "archive/原临床分析代码/revision_20261012/reference_results/content_raw_unit_contrast.csv": "2e9316d1b185688a069fe2eb1b4ca59a14fdfeb1d131cbe487109106f888a1fe",
    "archive/原临床分析代码/revision_20261012/reference_results/manuscript_scaling.json": "65ef036efa19c9a519d882ea1c0c1821efec1414d5ed39ad61747634ef9e48ee",
    "archive/原临床分析代码/revision_20261012/reference_results/source_clinical_content_contrast.csv": "ce8cf5567289ca083d5fc409c129f39d2db4369060f117cdd86903a8b7e4a748",
    "archive/原临床分析代码/revision_common.py": "99024c1d890a082916b02b87d08c4c3bc20304af22e60b26d6b8a5a8b6e2eb19",
    "archive/原临床分析代码/revision_runtime.py": "4b8fd7ad7c546ee11c0e2e4c11bbbbd8f3874f17e528418e81aa9e28682870ed",
    "data/questionnaire_analysis.csv": "5342bd167496989d0ebf083894d24547c50b210356715cdeceb49b76b5ad352b",
    "data/questionnaire_dictionary.json": "adfd6164b8db2483e01a214a921f47a3d1c3e333703233a0c57e48ce5d5cdd82",
    "archive/原分析汇总结果/revision_20261012/content_revision_audit.json": "87918f7b6f647e51a2648835d0a19ae5816023dd292c330650ef589e3d3d894b",
    "archive/原分析汇总结果/revision_20261012/ipw_summary.json": "3a60476b3c94c3ecc74490b24c6e3b9942798c90995934bcd47187d679f6fb49",
    "archive/原分析汇总结果/revision_20261015/post_review_analysis_audit.json": "9b35cd0eaa6267c64798d02b5767afff3db8521cf18344490732b5beb6df557e",
    "archive/原分析汇总结果/revision_20261015/spline_knots.json": "b61bb499485245ae7fb7aed955fd4c59537ececcccabb8733cee0b33d90d9d46",
    "archive/原分析汇总结果/revision_20261015/function_form_iqr_contrasts.csv": "eff60f76a986657c99164ff704b86227cbe77da544340c890e4ffda67ab640f7",
    "archive/原分析汇总结果/diagnostic_addendum_20261015/diagnostic_addendum_audit.json": "70b2008c032221f0032f5e1f4bedcad0714136ae76059c17c6195124a888bbb1",
    "source_data/analysis/core_26_model_catalog.csv": "45c940a01bb09815ba9559e00dc7d06931725d1a5cd889b26ba7c5ac5f0a2481",
    "source_data/analysis/scaling_register_full.csv": "1fe075fb4eed3098cf18ceee60dfa524562db765c8c02e3bbfca378fa14bff40"
}
ENTRY_POINTS = [
    'archive/原临床分析代码/revision_20261012/content_analysis.py',
    'archive/原临床分析代码/revision_20261012/ipw_analysis.py',
    'archive/原临床分析代码/post_review_20261015/post_review_analysis.py',
    'archive/原临床分析代码/post_review_20261015/diagnostic_addendum.py',
]
AUDITS = [
    'archive/原分析汇总结果/revision_20261012/content_revision_audit.json',
    'archive/原分析汇总结果/revision_20261012/ipw_summary.json',
    'archive/原分析汇总结果/revision_20261015/post_review_analysis_audit.json',
    'archive/原分析汇总结果/diagnostic_addendum_20261015/diagnostic_addendum_audit.json',
]
DISTRIBUTIONS = ['numpy', 'pandas', 'scipy', 'statsmodels', 'patsy',
                 'scikit-learn', 'threadpoolctl']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path,
                        help='Optional NEW JSON report outside this package; default writes only stdout.')
    args = parser.parse_args()
    checks = []
    for rel, expected in EXPECTED_SHA256.items():
        path = ROOT / rel
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        checks.append(dict(kind='fingerprint', path=rel, exists=path.is_file(),
                           expected_sha256=expected, actual_sha256=actual,
                           passed=actual == expected))
    for rel in ENTRY_POINTS:
        path = ROOT / rel
        try:
            tree = ast.parse(path.read_text())
            root_assignment = next(node.value for node in tree.body
                if isinstance(node, ast.Assign) and any(
                    isinstance(target, ast.Name) and target.id == 'ROOT'
                    for target in node.targets))
            root_expression = ast.unparse(root_assignment)
            expected_expression = 'Path(__file__).resolve().parents[3]'
            resolved = path.resolve().parents[3]
            ok = root_expression == expected_expression and resolved == ROOT
        except (OSError, SyntaxError, StopIteration) as exc:
            root_expression, resolved, ok = str(exc), None, False
        checks.append(dict(kind='entry_point_root', path=rel,
                           expression=root_expression,
                           resolved_root=str(resolved) if resolved else None,
                           passed=ok))
    initial = ROOT / 'archive/原分析汇总结果/revision_20261015'
    addendum = ROOT / 'archive/原临床分析代码/post_review_20261015/diagnostic_addendum.py'
    addendum_text = addendum.read_text() if addendum.is_file() else ''
    checks.append(dict(kind='addendum_explicit_initial_results',
                       path=str(initial.relative_to(ROOT)),
                       passed=all((initial/f).is_file() for f in [
                           'spline_knots.json', 'function_form_iqr_contrasts.csv'])
                           and "'--initial-results'" in addendum_text,
                       note='Use documented --initial-results; old default paper/evidence path is not present.'))
    scale = ROOT / 'source_data/analysis/scaling_register_full.csv'
    rows = list(csv.DictReader(scale.open())) if scale.is_file() else []
    checks.append(dict(kind='active_school_scaling', rows=len(rows),
                       populations=sorted({r['population'] for r in rows}),
                       passed=len(rows)==24 and {r['population'] for r in rows}=={'school'}))
    versions = {}
    for rel in AUDITS:
        path = ROOT / rel
        versions[rel] = json.loads(path.read_text()).get('versions', {}) if path.is_file() else None
    installed = {}
    for name in DISTRIBUTIONS:
        try:
            installed[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            installed[name] = None
    passed = all(check['passed'] for check in checks)
    report = dict(package_root=str(ROOT), scope='static inputs and path preflight only',
                  model_fitting_executed=False, bootstrap_executed=False,
                  analysis_modules_imported=False,
                  input_and_path_checks_passed=passed,
                  checks=checks, historical_recorded_versions=versions,
                  current_preflight_python=sys.version.split()[0],
                  current_installed_distribution_metadata=installed,
                  missing_current_distributions=[name for name, value in installed.items() if value is None],
                  current_environment_has_all_named_dependencies=all(installed.values()),
                  runtime_fit_validation_performed=False,
                  note='Input checks do not certify an executable numerical environment or numerical reproduction.')
    output = json.dumps(report, indent=2, ensure_ascii=False) + '\n'
    if args.report:
        target = args.report.expanduser().resolve()
        if target == ROOT or ROOT in target.parents:
            parser.error('--report must be outside this package.')
        if not target.parent.is_dir():
            parser.error('--report parent directory must already exist.')
        with target.open('x') as handle:
            handle.write(output)
    print(output, end='')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
