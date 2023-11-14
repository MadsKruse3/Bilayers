import os
from ase.io import read
from asr.core import read_json
import numpy as np
from pathlib import Path
from asr.utils.slidingequivalence import test_slide_equiv, slide_vector_for_bilayers


def parse_folders():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="Monolayer folders to analyse.")
    args = parser.parse_args()
    
    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]

    return folders

def listdirs(path):
    if type(path) != str:
        path = str(path)
    return [f"{path}/" + x for x in os.listdir(path) if os.path.isdir(f"{path}/{x}")]

def all_relax_done(folder):
    for sf in listdirs(folder):
        if not os.path.exists(f"{sf}/results-asr.relax_bilayer.json"):
            return False

    return True


def has_id(monolayer_folder):
    name = "/".join(str(monolayer_folder).split("/")[-3:])
    fullname = f"/home/niflheim2/cmr/C2DB-ASR/tree/{name}/info.json"
    if not os.path.exists(fullname):
        return False
    info = read_json(fullname)
    return "icsd_id" in info or "cod_id" in info or "doi" in info


def has_doi(monolayer_folder):
    name = "/".join(str(monolayer_folder).split("/")[-3:])
    fullname = f"/home/niflheim2/cmr/C2DB-ASR/tree/{name}/info.json"
    if not os.path.exists(fullname):
        return False
    info = read_json(fullname)
    return "doi" in info



def slide_equiv_p(bilayer1, bilayer2):
    _k = bilayer1.replace("--", "-/")
    _k = _k.replace("_-", "_/")
    _ksplit = _k.split("-")[2:]
    _kIz = len(_ksplit) == 3
    mat = _ksplit[0]
    vec = _ksplit[-1]
    
    _k2 = bilayer2.replace("--", "-/")
    _k2 = _k2.replace("_-", "_/")
    _k2split = _k2.split("-")[2:]
    _k2Iz = len(_k2split) == 3
    mat2 = _k2split[0]
    vec2 = _k2split[-1]

    slide_equiv = slide_vector_for_bilayers(bilayer1, bilayer2) is not None    
    return (mat == mat2 and _k2Iz == _kIz) or slide_equiv

def get_slide_equivs(monolayer_folder):
    sfs = listdirs(monolayer_folder)
    
    pairs = []
    for i1, sf1 in enumerate(sfs):
        for sf2 in sfs[i1+1:]:
            try:
                if slide_equiv_p(sf1, sf2):
                    pairs.append((sf1, sf2))
            except Exception:
                pass
    return pairs


def gap_type(results):
    gap = results["gap"]
    gap_dir = results["gap_dir"]
    
    if gap_dir <= gap:
        return "direct"
    else:
        return "indirect"

def gap(results):
    return min(results["gap"], results["gap_dir"])


def gap_change(results1, results2):
    
    gt1 = gap_type(results1)
    gt2 = gap_type(results2)

    gap1 = gap(results1)
    gap2 = gap(results2)

    if gt1 == gt2:
        return "no change", gap2 - gap1
    else:
        return f"{gt1} -> {gt2}", gap2- gap1


def get_dip(folder, area):
    dpath = f"{folder}/results-asr.gs.json"
    if not os.path.exists(dpath):
        raise ValueError(f"asr.gs does not exist for {folder}")

    data = read_json(dpath)
    dip = data["dipz"] / area * 1.602 * 1e-19 / 1e-10
    return dip


def get_dipole_switches(folder):
    data = []
    atoms = read(f"{folder}/structure.json")
    area = np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))

    for sf in listdirs(folder):
        fname = f"{sf}/results-asr.gs.json"
        if not os.path.exists(fname):
            continue
        
        dip = get_dip(sf, area)
        if dip is None:
            continue

        data.append((str(sf), dip))


    seen = []
    pairs = []
    for i1, (sf1, dip1) in enumerate(data):
        seen.append(i1)
        for i2, (sf2, dip2) in enumerate(data):
            if i2 in seen:
                continue
            try:
                b = slide_equiv_p(sf1, sf2)
            except Exception as e:
                print(f"WARNING: {sf1} OR {sf2} could not be parsed")
                continue
            if not slide_equiv_p(sf1, sf2):
                continue

            delta_dip = abs(dip1 - dip2)
            pairs.append((str(sf1), str(sf2), delta_dip))

    return pairs

        

def is_metal(folder):
    fname = f"{folder}/results-asr.gs.json"
    if not os.path.exists(fname):
        raise ValueError("Results for GS does not exist for {folder}")

    data = read_json(fname)
    return np.allclose(data["gap"], 0.0) and np.allclose(data["gap_dir"], 0.0)


def get_metal_dipole_switches(folder):
    data = []
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

        data.append((str(sf), dip))


    seen = []
    pairs = []
    for i1, (sf1, dip1) in enumerate(data):
        seen.append(i1)
        for i2, (sf2, dip2) in enumerate(data):
            if i2 in seen:
                continue
            try:
                b = slide_equiv_p(sf1, sf2)
            except Exception as e:
                print(f"WARNING: {sf1} OR {sf2} could not be parsed")
                continue
            if not slide_equiv_p(sf1, sf2):
                continue

            delta_dip = abs(dip1 - dip2)
            pairs.append((str(sf1), str(sf2), delta_dip))

    return pairs

        
