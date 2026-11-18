#!/usr/bin/env python3
"""Safely unpack downloaded OE3035 HydRA ZIPs into data/.

Examples:
    python src/import_data.py ~/Downloads/KCS_L1p000.zip ~/Downloads/KCS_L0p950.zip
    python src/import_data.py "C:\\Users\\me\\Downloads\\KCS_L1p000.zip" "C:\\Users\\me\\Downloads\\KCS_L0p950.zip"
"""
import sys
from pathlib import Path
from zipfile import ZipFile

ROOT=Path(__file__).resolve().parents[1]


def main():
    if len(sys.argv)!=3:
        raise SystemExit('Usage: python src/import_data.py KCS_L1p000.zip KCS_L0p950.zip')
    dest=ROOT/'data';dest.mkdir(parents=True,exist_ok=True)
    for p in map(Path,sys.argv[1:]):
        if not p.is_file():raise SystemExit(f'Missing downloaded file: {p}')
        if p.stem not in ('KCS_L1p000','KCS_L0p950'):
            raise SystemExit(f'Unexpected file: {p}. Expected both KCS hull ZIPs.')
        with ZipFile(p) as z:
            files=z.infolist()
            for info in files:
                name=info.filename.replace('\\','/')
                target=(dest/name).resolve()
                if not target.is_relative_to(dest.resolve()) or info.is_dir() and name.startswith('/'):
                    raise SystemExit(f'Unsafe ZIP path: {name}')
            z.extractall(dest)
        matches=list(dest.rglob(p.stem+'_hydra.json'))
        if len(matches)!=1:
            raise SystemExit(f'Unpacked {p} but could not locate precisely one {p.stem}_hydra.json. '
                             'Inspect the ZIP structure and place its matlab/ folder in data/.')
        print('Ready:',matches[0])
    print('Now run: python src/msi.py KCS_L0p950')

if __name__=='__main__':main()
