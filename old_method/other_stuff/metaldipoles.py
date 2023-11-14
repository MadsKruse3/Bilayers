from pathlib import Path
from utils import listdirs, slide_equiv_p
from asr.core import read_json
import numpy as np
from dipolez import get_dipzs, get_dip
import os
from ase.io import read


"""
We want to check 
1. Whether there are pairs of bilayers which are both metals, can be switched between, and have large dipole switching
2. Whether there are bilayers that are metallic and have large dipole
"""


def is_metal(folder):
    fname = f"{folder}/results-asr.gs.json"
    if not os.path.exists(fname):
        raise ValueError("Results for GS does not exist for {folder}")

    data = read_json(fname)
    return np.allclose(data["gap"], 0.0) and np.allclose(data["gap_dir"], 0.0)


def metallic_dipole_switch(folder):
    names_dips = []
    atoms = read(f"{folder}/structure.json")
    area = np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))

    for sf in listdirs(folder):
        fname = f"{sf}/results-asr.gs.json"
        if not os.path.exists(fname):
            continue
        if not is_metal(sf):
            continue
        dip = get_dip(sf, area)
        if dip is None:
            continue
            
        names_dips.append((str(sf), dip))

    seen = []
    pairs = [] # (name1, name2, delta dip)
    for i1, (sf1, dip1) in enumerate(names_dips):
        seen.append(i1)
        for i2, (sf2, dip2) in enumerate(names_dips):
            if i2 in seen:
                continue
                
            if not slide_equiv_p(sf1, sf2):
                continue

            delta_dip = abs(dip1 - dip2)
            pairs.append((sf1, sf2, delta_dip))

    return pairs


def printstuff(folders):
    pairs = []
    N = len(folders)
    i = 0
    for folder in folders:
        i += 1
        print(f"Analysing folder {i}/{N} with science", end="\r")
        pairs += metallic_dipole_switch(folder)

    pairs = sorted(pairs, key=lambda t: t[2])

    def gname(f):
        return f.split("/")[-1]
    
    L1 = max(map(lambda t: len(gname(t[0])), pairs))
    L2 = max(map(lambda t: len(gname(t[1])), pairs))

    for n1, n2, d in pairs:
        if abs(d) < 0.002 * 1.602 * 1e-19 / 1e-10 * 0.2 / 3.204:
            continue
        f1 = gname(n1)
        f2 = gname(n2)
        ds = str(round(d, 15))
        print(f1.ljust(L1), " <-> ", f2.ljust(L2), " :: ", ds.ljust(10), " C / m")

    
            


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="Monolayer folders to analyse.")

    args = parser.parse_args()
    
    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]

    printstuff(folders)
