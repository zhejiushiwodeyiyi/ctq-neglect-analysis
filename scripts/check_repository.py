"""Check a code-only CTQ snapshot without importing analysis modules or data."""
from pathlib import Path
from hashlib import sha256
import ast
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    'archive/原临床分析代码/revision_common.py',
    'archive/原临床分析代码/revision_runtime.py',
    'archive/原临床分析代码/revision_20261012/content_analysis.py',
    'archive/原临床分析代码/revision_20261012/ipw_analysis.py',
    'archive/原临床分析代码/post_review_20261015/post_review_analysis.py',
    'archive/原临床分析代码/post_review_20261015/diagnostic_addendum.py',
    'archive/原临床分析代码/post_review_20261015/dependencies/content_analysis_20261012.py',
    'code/build_current_figures.py',
    'code/figure_qa_tools/audit_panel_alignment.py',
    'code/specificity_20261020/ANALYSIS_PLAN.md',
    'code/specificity_20261020/README.md',
    'code/specificity_20261020/measurement_audit.py',
    'code/specificity_20261020/bounded_models.py',
    'code/specificity_20261020/ancillary_checks.py',
    'code/specificity_20261020/verify_run.py',
    'code/specificity_20261020/build_figure_s4.py',
    'code/mnar_20261020/ANALYSIS_PLAN.md',
    'code/mnar_20261020/README.md',
    'code/mnar_20261020/run_mnar_sensitivity.py',
    'code/build_figureS5_followup_flow.py',
]
ALLOWED_SUFFIXES = {'.py', '.md', '.json', '.txt'}
SKIP_PARTS = {'.git', '.venv', 'venv', '__pycache__', 'matplotlib_cache', 'font_cache'}


def main():
    failures = []
    files = [p for p in ROOT.rglob('*') if p.is_file()
             and not any(part in SKIP_PARTS for part in p.relative_to(ROOT).parts)]
    for path in files:
        rel = path.relative_to(ROOT).as_posix()
        if path.is_symlink():
            failures.append(f'Symlink is outside the published file scope: {rel}')
        if path.name != '.gitignore' and path.suffix not in ALLOWED_SUFFIXES:
            failures.append(f'Unexpected file type in code-only snapshot: {rel}')
        if path.stat().st_size > 1024 * 1024:
            failures.append(f'Unexpectedly large file: {rel}')
        if path.suffix == '.py':
            try:
                ast.parse(path.read_text(encoding='utf-8'), filename=rel)
            except (SyntaxError, UnicodeError) as exc:
                failures.append(f'Syntax or encoding problem in {rel}: {exc}')
    for rel in REQUIRED:
        if not (ROOT / rel).is_file():
            failures.append(f'Missing code dependency: {rel}')
    for rel in REQUIRED[2:6]:
        path = ROOT / rel
        if path.is_file() and path.resolve().parents[3] != ROOT:
            failures.append(f'Entry-point directory depth is wrong: {rel}')
    manifest_path = ROOT / 'SOURCE_MANIFEST.json'
    if not manifest_path.is_file():
        failures.append('Missing SOURCE_MANIFEST.json')
        manifest = {'copied_files': []}
    else:
        manifest = json.loads(manifest_path.read_text())
    for record in manifest['copied_files']:
        path = ROOT / record['path']
        if not path.is_file() or sha256(path.read_bytes()).hexdigest() != record['sha256']:
            failures.append(f'Copied source fingerprint changed: {record["path"]}')
    for record in manifest.get('adapted_files_20261021', []):
        path = ROOT / record['path']
        if not path.is_file() or sha256(path.read_bytes()).hexdigest() != record['sha256']:
            failures.append(f'Adapted source fingerprint changed: {record["path"]}')
    for path in files:
        if path.suffix in ALLOWED_SUFFIXES and path.name != 'check_repository.py':
            text = path.read_text(encoding='utf-8')
            if '/Users/' in text or '-----BEGIN PRIVATE KEY-----' in text:
                failures.append(f'Machine-specific path or private key marker found: {path.relative_to(ROOT)}')
    report = {
        'passed': not failures,
        'checked_files': len(files),
        'python_scripts_parsed_without_import': sum(p.suffix == '.py' for p in files),
        'unchanged_copied_sources': len(manifest['copied_files']),
        'checked_adapted_sources_20261021': len(manifest.get('adapted_files_20261021', [])),
        'participant_data_included': False if not failures else 'inspect failures',
        'analysis_modules_imported': False,
        'models_fitted': 0,
        'bootstrap_executed': False,
        'numerical_reproduction_verified': False,
        'failures': failures,
    }
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if not failures else 1


if __name__ == '__main__':
    sys.exit(main())
