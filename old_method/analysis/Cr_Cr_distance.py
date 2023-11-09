from ase.neighborlist import NeighborList
import numpy as np
from asr.core import read_json
from asr.utils.bilayerutils import translation
from ase.io import read    
from asr.core import read_json
    
def mindist(folder):

    distances = []

    atom = read(f'{folder}/structure.json') 
    pos = atom.get_positions()
    distances.append(np.linalg.norm(pos[0] - pos[8]))
    distances.append(np.linalg.norm(pos[1] - pos[8]))
    distances.append(np.linalg.norm(pos[0] - pos[9]))
    distances.append(np.linalg.norm(pos[1] - pos[9]))
    print(atom)
    
    return distances

if __name__ == "__main__":
    from argparse import ArgumentParser
    from pathlib import Path
    
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*")
    
    args = parser.parse_args()
    
    if len(args.folders) > 0:
        folders = [Path(f) for f in args.folders if Path(f).is_dir()]
    else:
        folders = [Path(".")]


    for folder in folders:
        dist = mindist(folder)

        
        print(dist)
