"""
This script accepts a list of folders that contains bilayers
as argument and collects the data from these folders.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from asr.core import read_json
from pathlib import PosixPath

def bindingenergy(folder, monolayer_energy, cell):
    p = f"{folder}/results-asr.relax_bilayer.json"
    if not os.path.exists(p):
        return None
    data = read_json(p)
    energy = data["energy"]
    cellarea = np.linalg.norm(np.cross(cell[0], cell[1]))
    return (2 * monolayer_energy - energy) / cellarea


def bindinglength(folder, atoms):
    p = f"{folder}/results-asr.relax_bilayer.json"
    if not os.path.exists(p):
        return None
    data = read_json(p)
    w = np.max(atoms.positions[:, 2]) - np.min(atoms.positions[:, 2])
    return data["optimal_height"] - w


def get_topfolder(folder, n):
    """Return topfolder from a path.

    /home/niflheim/asbra/tree/AB2/MoS2/MoS2-b3b/MoS2-2-1_0_0_1-0_0/
    --->
    n = 1 -> /home/niflheim/asbra/tree/AB2/MoS2/MoS2-b3b
    n = 2 ->  /home/niflheim/asbra/tree/AB2/MoS2
    and so on.
    n = 0 returns folder.
    n negative returns path starting from the root.
    n = -1 -> /home
    n = -2 -> /home/niflheim
    and so on.
    """
    if n == 0:
        return folder

    return "/" + "/".join(list(filter(lambda x: x != "", folder.split("/")))[:-n])


def construct_label(folder, label_level):
    return list(filter(lambda x: x!= "", folder.split("/")))[-(label_level)]


def get_energy_length(folder, vdw_correction=True):
    if type(folder) == PosixPath:
        folder = str(folder.absolute())
    topfolder = get_topfolder(folder, 1)
    monolayer = read(f"{topfolder}/structure.json")
    monolayer_energy = monolayer.get_potential_energy()
    if vdw_correction:
        vdwp = f"{topfolder}/vdw_e.npy"
        if not os.path.exists(vdwp):
            return None, None
        vdw_e = np.load(vdwp).item()
        monolayer_energy += vdw_e


    cell = monolayer.get_cell()
    energy = bindingenergy(folder, monolayer_energy, cell)
    length = bindinglength(folder, monolayer)
    return energy, length


def collectdata(folders, label_level, vdw_correction=True):
    """Collect bilayer data from folders.
    
    label_level indicates how the label is determined
    from the folder structure.

    e.g. 

    /home/niflheim/asbra/tree/AB2/MoS2/MoS2-b3b/MoS2-2-1_0_0_1-0_0/
    --->
    label_level = 1 -> MoS2-2-1_0_0_1-0_0
                = 2 -> MoS2-b3b
                = 3 -> MoS2
    and so on.
    """
    # data has label -> list of binding length/binding energy tuples
    data = {}

    for folder in folders:
        energy, length = get_energy_length(folder, vdw_correction)
        label = construct_label(folder, label_level)

        if energy is None or length is None:
            continue

        if label not in data:
            data[label] = {"lengths": [],
                           "energies": []}
            
        data[label]["lengths"].append(length)
        data[label]["energies"].append(energy)
        
    return data
        

def flatten(data):
    """Takes a data dict and flattens it.

    returns energies, lengths
    """
    energies = []
    lengths = []
    for label, dct in data.items():
        energies.extend(dct["energies"])
        lengths.extend(dct["lengths"])

    return energies, lengths
        

if __name__ == "__main__":
    import os
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*")
    parser.add_argument("-n", "--label_level", type=int, default=0)

    args = parser.parse_args()

    args.folders = list(filter(lambda x: os.path.isdir(x), args.folders))

    
    data = collectdata(args.folders, args.label_level)
    print(data)
    
