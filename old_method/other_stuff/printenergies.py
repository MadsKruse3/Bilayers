from pathlib import Path
from asr.core import read_json
from ase.io import read
from argparse import ArgumentParser
import os
from collectlayerdata import bindingenergy, construct_label

parser = ArgumentParser()
parser.add_argument("-f", "--folders", nargs="*")

args = parser.parse_args()

if args.folders is not None:
    args.folders = list(filter(lambda x: os.path.isdir(x), args.folders))
else:
    args.folders = ["."]



labels = []
energies = []
for folder in args.folders:
    p = Path(folder)
    monolayer = read(f"{folder}/structure.json")
    ml_energy = monolayer.get_potential_energy()
    cell = monolayer.get_cell()
    for x in p.glob("*/"):
        if not x.is_dir():
            continue
        x = x.absolute()
        energy = bindingenergy(str(x), ml_energy, cell)
        label = construct_label(str(x), 1)
        labels.append(label)
        energies.append(energy)


maxl = max(map(len, labels))
labs_es = sorted(list(zip(labels, energies)), key=lambda t: -t[1])
for label, energy in labs_es:
    estr = str(round(energy, 6))
    estr = estr.ljust(8, "0")
    print(f"{label}".ljust(maxl) + " binding energy " + f"= {estr} eV / Å^2")
        
    
