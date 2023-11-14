from pathlib import Path
from utils import listdirs
import os
from asr.core import read_json
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read

def get_dipzs(folder):
    dips = []
    large = []
    atoms = read(f"{folder}/structure.json")
    area = np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))
    for sf in listdirs(str(folder)):
        dpath = f"{sf}/results-asr.gs.json"
        if not os.path.exists(dpath):
            continue
        data = read_json(dpath)
        dip = data["dipz"] / area * 1.602 * 1e-19 / 1e-10
        if abs(dip) > 0.002 * 1.602 * 1e-19 / 1e-10 * 2 / 3.204:
            large.append(("/".join((str(sf).split("/")[-2:])), dip))
        dips.append(dip)

    return dips, large

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="Monolayer folders to analyse.")

    args = parser.parse_args()
    
    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]
        
    dips = []
    larges = []
    for folder in folders:
        dip, large = get_dipzs(folder)
        dips.extend(dip)
        larges.extend(large)

    L = max(map(lambda t: len(t[0]), larges))
    for n, d in larges:
        print(f"{n}".ljust(L), f":{d}")
    print(f"Number of materials : {len(dips)}")
    print(f"Number of large dips: {len(larges)}")
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4), dpi=200)
    axes[0].hist(dips, bins=40)
    axes[0].set_xlabel(r"Dipole density [C m$^{-1}$]")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Full range")
    axes[1].hist(list(filter(lambda d: abs(d) > 2 / 3.204 * 0.002 * 1.602 * 1e-19 / 1e-10, dips)), bins=40)
    axes[1].set_xlabel(r"Dipole density [C m$^{-1}$]")
    axes[1].set_ylabel("Count")
    axes[1].set_title(r"|d| > 2e-12  [C m$^{-1}$]")
    plt.savefig("dipolehistogram.png")
    # plt.show()
