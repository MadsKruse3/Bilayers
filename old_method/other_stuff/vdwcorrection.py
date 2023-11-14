from ase.calculators.calculator import get_calculator_class
from ase.io import read
import numpy as np 
import os


def calc_vdw(folder):
    Calculator = get_calculator_class("dftd3")
    calc = Calculator()

    atoms = read(f"{folder}/structure.json")
    atoms.set_calculator(calc)

    e = atoms.get_potential_energy()
    return e


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        folders = sys.argv[1:]
    else:
        folders = ["."]

    N = len(folders)
    for i, folder in enumerate(folders):
        if N > 1:
            print(f"Calculating vdW correction for material {i}/{N}", flush=True, end="\r")
        if os.path.exists(f"{folder}/vdw_e.npy"):
            continue
        vdw_e = calc_vdw(folder)
        np.save(f"{folder}/vdw_e.npy", vdw_e)
    print(f"Calculating vdW correction for material {i}/{N}", flush=True)

