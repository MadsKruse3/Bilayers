#from analyses import analysis_objs
from analyses_mag import analysis_objs
from ase.io import read
from utils import slide_equiv_p
from typing import List
from pathlib import Path
from collectlayerdata import get_energy_length
from asr.utils.slidingequivalence import test_slide_equiv

def filter_monolayers(folders, fname: str = "results-asr.gs.json"):

    def file_exists(monolayer_folder):
        path = "/".join([x for x in monolayer_folder.split("/") if x != ""][-3:])
        c2db_path = "/home/niflheim2/cmr/C2DB-ASR/tree/" + path + "/" + fname
        c2db_path2 = "/home/niflheim2/cmr/C2DB-ASR/icsd_cod_materials/tree/" + path + "/" + fname
        return Path(c2db_path).is_file() or Path(c2db_path2).is_file()

    def atmost_eight_atoms(monolayer_folder):
        atoms = read(f"{monolayer_folder}/structure.json")
        return len(atoms) <= 8

    thefilter = lambda f: file_exists(f) and atmost_eight_atoms(f)

    return list(filter(thefilter, folders))

def accepted_bilayers(monolayer_folder: str, only_most_stable: bool,
                      fname: str = "results-asr.gs.json", return_all=False) -> List[str]:
    """Calculate list of bilayers to analyse.

    If #atoms in monolayer > 8 return empty list.

    Else apply stability criterion.

    Also filter on required data.
    """
    ml_atoms_path = f"{monolayer_folder}/structure.json"
    assert Path(ml_atoms_path).is_file()
    ml_atoms = read(ml_atoms_path)

    vdw_path = f"{monolayer_folder}/vdw_e.npy"
    assert Path(vdw_path).is_file(), f"vdw correction does not exist for {monolayer_folder}. This is equired!"

    if len(ml_atoms) > 8 and not return_all:
        return []

    bilayer_folders = [x for x in Path(monolayer_folder).iterdir() if x.is_dir()]
    selected_bilayers = stability_criterion(bilayer_folders, only_most_stable, return_all=return_all)

    return selected_bilayers

def stability_criterion(bilayer_folders: List[Path], only_most_stable: bool,
                        fname: str = "results-asr.gs.json", return_all: bool = False) -> List[Path]:
    """Return those folders that satisfy stability criterion."""
    nmats = 5
    deltaE = 0.002
    cutoff = 0.15

    energies = []
    for folder in bilayer_folders:
        # Read stability data
        # Check for existence of file
        relax_path = f"{folder}/results-asr.relax_bilayer.json"
        
        # Skip this bilayer if relaxation is not done
        if not Path(relax_path).is_file():
            continue

        if fname is not None:
            data_path = f"{folder}/{fname}"
            if not Path(data_path).is_file():
                continue

        energy, length = get_energy_length(folder)
        assert energy is not None

        energies.append((folder, energy))

    if len(energies) == 0:
        return []

    if return_all and only_most_stable:
        return [max([(str(a), b) for a, b in energies], key=lambda t:t[1])[0]]


    maxE = max([e for f, e in energies])
    if maxE > cutoff:
        return []

    selected = filter(lambda t: abs(t[1] - maxE) < deltaE, energies)
    selected = sorted(list(selected), key=lambda t: t[1])[:nmats]
    if only_most_stable:
        selected = [selected[-1]]

    return [str(f.resolve()) for f, e in selected]

def get_switchables(bilayer_folders: List[str]) -> List[str]:
    switches = []
    for i, bl1 in enumerate(bilayer_folders):
        for bl2 in bilayer_folders[i+1:]:
            if slide_equiv_p(bl1, bl2):
                switches.append((bl1, bl2))

    return switches

def is_invertable(folder):
    v = test_slide_equiv(bilayer)
    if v is None:
        return False
    elif np.allclose(v, 0.0):
        return False
    else:
        return True

def main(monolayer_folders, analyses):
    for monolayer_folder in monolayer_folders:
        bilayers = accepted_bilayers(monolayer_folder, only_most_stable=False)
        switchables = get_switchables(bilayers)
        invertables = [bilayer for bilayer in bilayers if is_invertable(bilayer)]

        for analysis in analyses:
            analysis.run_monolayer(monolayer_folder)

            analysis.run_bilayers(bilayers)

            analysis.run_switchables(switchables)

            analysis.run_invertables(invertables)

            analysis.finalize_monolayer()

    for analysis in analyses:
        analysis.save()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*")
    parser.add_argument("-p", "--plot", action="store_true")
    parser.add_argument("-s", "--save", action="store_true")
    parser.add_argument("-n", "--plotindex", type=int, default=-1)
    args = parser.parse_args()
    
    if not args.plot:
        assert len(args.folders) > 0
        monolayer_folders = filter_monolayers(args.folders)
        main(monolayer_folders, analysis_objs)
    else:
        if args.plotindex != -1:
            analysis_objs[args.plotindex].plot(args.save)
        else:
            for obj in analysis_objs:
                obj.plot(args.save)
