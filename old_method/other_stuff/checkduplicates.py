from ase.io import read
from asr.core import read_json
import numpy as np
import os
import sys
from asr.stack_bilayer import atomseq


folders = sys.argv[1:]

structures = [read(f"{folder}/structure.json") for folder in folders]


rmsd_tol = 0.3

for i1, (f1, s1) in enumerate(zip(folders, structures)):
    for i2, (f2, s2) in enumerate(zip(folders, structures)):
        if i1 == i2:
            continue
        
        if atomseq(s1, s2, full=True, rmsd_tol=rmsd_tol):
            print(f"{f1} == {f2}")
            
        
