"""
Take as input a number of monolayer folders.

Analyse every subfolder (bilayer) and determine if the
band gap or band character has changed.
"""
import sys
import os
from asr.core import read_json
from pathlib import Path
from enum import Enum, auto
from collections import namedtuple
import numpy as np
from minimumlayerdistance import mindist
from utils import listdirs, all_relax_done
from stackingworkflow import stability_criterion
from utils import has_doi

class GapType(Enum):
    DIRECT   = 1 # auto()
    INDIRECT = 2 # auto()

class GapTypeChange(Enum):
    NOCHANGE        = 1 # auto()
    DIRECT2INDIRECT = 2 # auto()
    INDIRECT2DIRECT = 3 # auto()

GapData = namedtuple("GapData", ["bilayername"     ,
                                 "gap_change"      ,
                                 "gapdir_change"   ,
                                 "gaptype_change"  ,
                                 "minimum_distance",
                                 "monofolder"])


def makesummary(data):
    gapchanges = [d.gap_change for d in data]
    gapdirchanges = [d.gapdir_change for d in data]
    distances = [d.minimum_distance for d in data]

    Nnochange = len([d for d in data if d.gaptype_change == GapTypeChange.NOCHANGE])
    N2indirect = len([d for d in data if d.gaptype_change == GapTypeChange.DIRECT2INDIRECT])
    N2direct = len([d for d in data if d.gaptype_change == GapTypeChange.INDIRECT2DIRECT])
    
    dg = pf(np.corrcoef(gapchanges, distances)[0, 1])
    dgd = pf(np.corrcoef(gapdirchanges, distances)[0, 1])

    b = [d for d in data if d.gaptype_change != GapTypeChange.NOCHANGE]
    b1 = [d.gap_change for d in b]
    d1 = [d.minimum_distance for d in b]
    dg1 = pf(np.corrcoef(b1, d1)[0, 1])
    

    return {"Total #bilayers        ": Nnochange + N2indirect + N2direct,
            "No Change              ": Nnochange,
            "Direct   -> Indirect   ": N2indirect,
            "Indirect -> Direct     ": N2direct,
    }



def bandgapchanges(folder, verbose=False, fast=False):
    """This function analysis band gap changes for a given monolayer.

    First checks if the monolayer gs data is available in the cmr tree.
    If it isn't, an empty list is returned and the missing data is reported.

    Calculates how much the band gaps (minimum gap and direct gap) have changes.
    Calculates whether the minimum gap has changed character (direct/indirect).

    If any data is missing from some of the bilayers, this is reported.

    """
    monolayer_path = Path("/home/niflheim2/cmr/C2DB-ASR/tree/" + "/".join(str(folder).split("/")[-3:]) + "/results-asr.gs.json")
    if not monolayer_path.exists():
        if verbose:
            print(f"No data for {monolayer_path}")
        return []

    monodata    = read_json(monolayer_path)
    monogap     = monodata["gap"]
    monogap_dir = monodata["gap_dir"]

    def getcharacter(data):
        if data["gap"] >= data["gap_dir"]:
            return GapType.DIRECT
        else:
            return GapType.INDIRECT

    def changedescriptor(char1, char2):
        if char1 == char2:
            return GapTypeChange.NOCHANGE  # "No Change"
        elif char1 == GapType.DIRECT:
            return GapTypeChange.DIRECT2INDIRECT  # "Direct -> Indirect"
        else:
            return GapTypeChange.INDIRECT2DIRECT  # "Indirect -> Direct"
    
    monocharacter = getcharacter(monodata)
    
    gapdata = []
    alldone = True
    selected = [str(s) for s in stability_criterion(folder)]
    for subfolder in selected:
        matname = "/".join(subfolder.split("/")[-2:])
        gsdatapath = f"{subfolder}/results-asr.gs.json"
        if not os.path.exists(gsdatapath):
            alldone = False
            continue
        data = read_json(gsdatapath)

        gap       = data["gap"]
        gap_dir   = data["gap_dir"]
        character = getcharacter(data)
        dist      = -1 if fast else mindist(subfolder) 

        gapdata.append(GapData(matname, gap - monogap, gap_dir - monogap_dir, changedescriptor(monocharacter, character), dist, folder))

    if not alldone and verbose:
        print(f"Not all bilayer gs calculations were done for {'/'.join(str(monolayer_path).split('/')[6:-1])}")

    return gapdata


def pf(f):
    s = str(round(f, 4))
    if f >= 0:
        s = "+" + s
    return s


def printtables(data):
    # print(f"Number of materials {c}")

    gapchangeds = [t for t in data if t.gaptype_change != GapTypeChange.NOCHANGE]
    
    # gapchangeds = data
    desclen = max(len(t.bilayername) for t in gapchangeds) + 10
    # for d in gapchangeds:
    #     print("#" * (desclen + 5))
    #     namelen = len(d.bilayername)
    #     delta = desclen - namelen
    #     pre = delta // 2
    #     post = delta - pre
    #     print(":::".ljust(pre), d.bilayername.ljust(namelen + post), ":::")
    #     print("Type change       :", d.gaptype_change)
    #     print("Gap change        :", d.gap_change)
    #     print("Direct gap change :", d.gapdir_change)
    #     print("#" * (desclen + 5))
    #     print("")

    desclen = max(len(t.bilayername) for t in gapchangeds) + 3

    def change2str(t):
        if t == GapTypeChange.DIRECT2INDIRECT:
            return "Direct   -> Indirect"
        elif t == GapTypeChange.INDIRECT2DIRECT:
            return "Indirect -> Direct"
        elif t == GapTypeChange.NOCHANGE:
            return "     No Change    "
        else:
            raise ValueError()
            

    w1 = desclen
    w2 = 25
    w3 = 18
    w4 = 25
    def rowprint(name, t, f1, f2, f3):
        print(name.ljust(w1), change2str(t).ljust(w2), pf(f1).ljust(w3), pf(f2).ljust(w4), pf(f3))
        
    print("")
    print("")
    header = "".join(["Descriptor".ljust(w1), "Type change".ljust(w2), "Gap change [eV]".ljust(w3), "Direct gap change [eV]".ljust(w4), "Minimum distance [Å]"])
    print(header)
    print("-"*len(header))
    for d in gapchangeds:
        rowprint(d.bilayername.ljust(desclen), d.gaptype_change,
                 d.gap_change, d.gapdir_change, d.minimum_distance)


def doidata(data):
    dat = [d for d in data if has_doi(d[-1])]

    return dat
        
        
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="Monolayer folders to analyse.")
    parser.add_argument("-t", "--table", action="store_true", help="Show table of results.")
    parser.add_argument("-p", "--plot", action="store_true", help="Plot gap changes in histogram.")
    parser.add_argument("-f", "--fast", action="store_true", help="Do fast analysis w/o dist-gap change correlation.")
    parser.add_argument("-s", "--save", type=str, default="", help="Name of file to save histogram in. If no string is given, don't save.")

    args = parser.parse_args()
    
    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]
        
    data = []
    c = 0

    N = len(folders)
    i = 0
    for folder in folders:
        i += 1
        print(f"Processing folder {i}/{N}", end="\r")
        if not all_relax_done(folder):
            continue
        s = bandgapchanges(folder, fast=args.fast)
        if len(s) > 0:
            c += 1
        data.extend(s)

    gapchanges = [d.gap_change for d in data]
    gapdirchanges = [d.gapdir_change for d in data]
    distances = [d.minimum_distance for d in data]


    def makesummary(data):
        Nnochange = len([d for d in data if d.gaptype_change == GapTypeChange.NOCHANGE])
        N2indirect = len([d for d in data if d.gaptype_change == GapTypeChange.DIRECT2INDIRECT])
        N2direct = len([d for d in data if d.gaptype_change == GapTypeChange.INDIRECT2DIRECT])

        dg = pf(np.corrcoef(gapchanges, distances)[0, 1])
        dgd = pf(np.corrcoef(gapdirchanges, distances)[0, 1])

        b = [d for d in data if d.gaptype_change != GapTypeChange.NOCHANGE]
        b1 = [d.gap_change for d in b]
        d1 = [d.minimum_distance for d in b]
        dg1 = pf(np.corrcoef(b1, d1)[0, 1])


        return {"Total                  ": Nnochange + N2indirect + N2direct,
                "No Change              ": Nnochange,
                "Direct   -> Indirect   ": N2indirect,
                "Indirect -> Direct     ": N2direct,
                "C(Dist, Gap change)    ": dg,
                "C(Dist, D. Gap change) ": dgd,
                "C2(Dist, Gap change)   ": dg1
                }

    summ = makesummary(data)

    for k, v in summ.items():
        print(f"{k} :: {v}")

    
    print("")
    dat = [d for d in data if has_doi(d.monofolder)]
    summ = makesummary(dat)
    for k, v in summ.items():
        print(f"{k} :: {v}")

    for d in dat:
        if d.gaptype_change != GapTypeChange.NOCHANGE:
            print(d.monofolder)

    if args.table:
        printtables(data)

    def histofunc(data, ax):
        counts, bins = np.histogram(data)
        ax.hist(bins[1:], bins, weights=counts / np.sum(counts) * 100)

    if args.plot or args.save != "":
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(nrows=1, ncols=2)
        histofunc(gapchanges, axes[0])
        axes[0].set_title("Gap Change")
        axes[0].set_ylabel("Count (%)")
        axes[0].set_xlabel("Energy [eV]")
        histofunc(gapdirchanges, axes[1])
        axes[1].set_title("Direct Gap Change")
        axes[1].set_xlabel("Energy [eV]")
    if args.save != "":
        plt.savefig(args.save)
    if args.plot:
        plt.show()
