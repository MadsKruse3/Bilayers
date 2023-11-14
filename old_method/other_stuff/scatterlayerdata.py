import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from asr.core import read_json
from collectlayerdata import collectdata
import matplotlib.pyplot as plt

def scatter(data, include_label=True):
    plt.figure(figsize=(8,3 * 8 / 5), dpi=200)
    for label, dct in data.items():
        energies = dct["energies"]
        lengths = dct["lengths"]
        if include_label:
            plt.scatter(energies, lengths, label=label, s=5)
        else:
            plt.scatter(energies, lengths, s=5, c="tab:blue")
        
    plt.xlabel("Energies [eV / Å^2]", fontsize=12)
    plt.ylabel("Interlayer distance [Å]", fontsize=12)
    if include_label:
        plt.legend()
        


if __name__ == "__main__":
    import os
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*")
    parser.add_argument("-n", "--label_level", type=int, default=0)
    parser.add_argument("-s", "--save", type=str, default="")
    args = parser.parse_args()
    args.folders = list(filter(lambda x: os.path.isdir(x), args.folders))

    data = collectdata(args.folders, args.label_level)
    
    scatter(data, args.label_level!=-1)
    if args.save != "":
        if args.save.endswith(".png"):
            name = args.save
        else:
            name = args.save + ".png"
        plt.savefig(name)
    else:
        plt.show()
