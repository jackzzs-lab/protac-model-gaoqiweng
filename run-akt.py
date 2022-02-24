#!/usr/bin/env python3
from pathlib import Path
import os, sys

input_path = Path.cwd()/'data'/'akt'
out_path = Path.cwd()/'data'/'akt_results'
cpus_each = 3
for i, ligand in enumerate(input_path.iterdir()):
    if ligand.is_dir():
        with open(ligand/'site_info.txt', 'r') as f:
            os.system(f"python main.py -irec {ligand/'receptor.pdb'} -ilig {ligand/'target.pdb'} -site={f.readlines()[0]} -ismi {ligand/'protac.smi'} -o {out_path} -cpu {cpus_each}")
            sys.exit(250)
