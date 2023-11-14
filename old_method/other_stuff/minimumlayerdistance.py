from ase.neighborlist import NeighborList
import numpy as np
from asr.core import read_json
from asr.utils.bilayerutils import translation

def construct_bilayer(path, h=None):
    from ase.io import read
    from asr.core import read_json
    
    top_layer = read(f'{path}/toplayer.json')
    minz = np.min(top_layer.positions[:, 2])
    tol = 0.1
    lowest_ind = np.abs(top_layer.positions[:, 2] - minz) < tol
    tags = np.zeros(len(top_layer))
    tags[lowest_ind] = 1
    top_layer.set_tags(tags)
    base = read(f'{path}/../structure.json')
    maxz = np.max(base.positions[:, 2])
    highest_ind = np.abs(base.positions[:, 2] - maxz) < tol
    tags = np.zeros(len(base))
    tags[highest_ind] = -1
    base.set_tags(tags)

    t = np.array(read_json(f'{path}/translation.json')['translation_vector']).astype(float)
    if h is None:
        h = read_json(f'{path}/results-asr.relax_bilayer.json')['optimal_height']

    return translation(t[0], t[1], h, top_layer, base), h


def mindist(folder):
    bilayer, h = construct_bilayer(folder)

    radii = np.zeros(len(bilayer))
    tags = bilayer.get_tags()
    inds = np.logical_not(np.isclose(tags, 0.0))
    radii[inds] = 1.2 * h

    nl = NeighborList(radii)
    nl.update(bilayer)

    d = 100
    for index in [i for i, b in enumerate(inds) if b]:
        indices, offsets = nl.get_neighbors(index)
        
        pos1 = bilayer.positions[index]
        tag = tags[index]
        assert not np.allclose(tag, 0.0)
        for i, o in zip(indices, offsets):
            if tags[i] == tag or np.allclose(tags[i], 0.0):
                continue

            pos2 = bilayer.positions[i] + o @ bilayer.get_cell()
            dist = np.linalg.norm(pos1 - pos2)
            if dist < d:
                d = dist
    
    return d




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
        mindist(folder)

