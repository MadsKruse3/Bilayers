from pathlib import Path
from utils import listdirs
import os
from asr.core import read_json
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from gpaw import GPAW

def get_gap_type(folder):
    name_material = []
    atoms = read(f"{folder}/structure.json")
    dpath = f"{folder}/results-asr.gs.json"
    if os.path.exists(dpath):
        data = read_json(dpath)
        gap = data["gap"]
        if 'gap' in locals():
            name_material.append(("/".join((str(folder).split("/")[-2:])), gap ))
        else:
            gap = 'missing'
            name_material.append(("/".join((str(folder).split("/")[-2:])), gap ))
   
        return gap, name_material
    else:
        name_material.append(("/".join((str(folder).split("/")[-2:])), 'missing' ))
        return 'missing', name_material

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="Monolayer folders to analyse.")

    args = parser.parse_args()
    
    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]

    no_gaps = []
    gaps = []

    for folder in folders:
        if os.path.isfile(f'{folder}/results-asr.bilayer_magnetism.json'): 
#        if os.path.isfile(f'{folder}/results-asr.gs.json'): 
            gap, name_material = get_gap_type(folder)

            if gap == 0.0:
                no_gaps.append(name_material)
                print(name_material)
            else:
                gaps.append(name_material)
                #print(name_material)
                    
    print('non-gapped magnetic bilayers:',len(no_gaps))
    print('gapped magnetic bilayers:', len(gaps))
                
