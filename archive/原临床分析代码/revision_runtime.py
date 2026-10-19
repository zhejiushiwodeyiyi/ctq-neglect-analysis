"""Portable output guards shared by current revision analysis entry points."""
from pathlib import Path
import argparse

def ensure_external_output(root, requested):
    root=Path(root).resolve(); output=Path(requested).expanduser().resolve()
    if output==root or root in output.parents:
        raise ValueError('Choose a NEW output directory outside this frozen submission package.')
    if output.exists():
        raise FileExistsError(f'Output already exists; choose a new directory: {output}')
    output.mkdir(parents=True)
    return output

def configure_output(root, description):
    parser=argparse.ArgumentParser(description=description)
    parser.add_argument('--output',type=Path,required=True,help='New output directory outside the submission package')
    args=parser.parse_args()
    try:return ensure_external_output(root,args.output)
    except (ValueError,FileExistsError) as exc:parser.error(str(exc))
