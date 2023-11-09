from abc import ABC, abstractmethod
from asr.core import read_json, ASRResult
from typing import List
from pathlib import Path
import numpy as np
from ase.io import read
from collectlayerdata import get_energy_length
from utils import slide_equiv_p
from asr.utils.slidingequivalence import test_slide_equiv 
import os
import matplotlib.pyplot as plt
from myqueue.task import task as mqtask
from myqueue.task import State

def task(*args, **kwargs):
    """Get MyQueue task instance."""
    name = kwargs.get("name") or args[0]
    if "creates" not in kwargs:
        kwargs["creates"] = [f"results-{name}.json"]
    return mqtask(*args, **kwargs)

def accepted_bilayers(monolayer_folder: str, only_most_stable: bool,
                      fname: str = "results-asr.gs.json", return_all=False) -> List[str]:
    """Calculate list of bilayers to analyse.

    If #atoms in monolayer > 8 return empty list.

    Else apply stability criterion.

    Also filter on required data.
    """
    ml_atoms_path = f"{monolayer_folder}/structure.json"
    #assert Path(ml_atoms_path).is_file()
    ml_atoms = read(ml_atoms_path)

    vdw_path = f"{monolayer_folder}/vdw_e.npy"
    assert Path(vdw_path).is_file(), f"vdw correction does not exist for {monolayer_folder}. This is equired!"

    if len(ml_atoms) > 8 and not return_all:
        return []

    bilayer_folders = [x for x in Path(monolayer_folder).iterdir() if x.is_dir()]
    most_stable_bilayers = stability_criterion(bilayer_folders, only_most_stable, return_all=return_all)
    selected_bilayers = dynamical_stability(most_stable_bilayers)

    #return most_stable_bilayers
    return selected_bilayers

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

def dynamical_stability(most_stable_bilayers: str) -> List[str]:
    dynamically_stable_bilayers = []
    for folder in most_stable_bilayers:
        if Path(f'{folder}/Stability.txt').is_file():
            #stability = np.load(f'{folder}/Stability.txt', allow_pickle=True)
            stability = np.loadtxt(f'{folder}/Stability.txt', dtype=str)
            if stability == 'Stable':
                dynamically_stable_bilayers.append(folder)
            if stability == 'Unstable-shifted':
                dynamically_stable_bilayers.append(folder)

    return dynamically_stable_bilayers

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
                        
def get_failed_calculations(monolayer_folder, selected_bilayers):
    
    tasks_interlayer = []
    tasks_anisotropy = []
    tasks_anisotropy_U = []

    for bilayer_folder in selected_bilayers:
        interlayer_task = None
        anisotropy_task = None
        anisotropy_U_task = None

        interlayer_task = task("asr.interlayer_magnetic_exchange",
                               resources="24:8h", 
                               folder=bilayer_folder)
        anisotropy_task = task("asr.magnetic_anisotropy",
                                   resources="24:8h", 
                                   folder=bilayer_folder)
        anisotropy_U_task= task("asr.magnetic_anisotropy_U",
                                    resources="24:8h", 
                                    folder=bilayer_folder)

        if interlayer_task.read_state_file() == State.FAILED:
            tasks_interlayer.append(bilayer_folder)
            print(bilayer_folder)
        if anisotropy_task.read_state_file() == State.FAILED:
            tasks_anisotropy.append(bilayer_folder)
            print(bilayer_folder)
        if anisotropy_U_task.read_state_file() == State.FAILED:
            tasks_anisotropy_U.append(bilayer_folder)
            print(bilayer_folder)

    #anisotropy_monolayer_task = task("asr.magnetic_anisotropy",
    #                       resources="24:8h", 
    #                       folder=monolayer_folder)

    anisotropy_U_monolayer_task= task("asr.magnetic_anisotropy_U",
                            resources="24:8h", 
                            folder=monolayer_folder)


    return tasks_interlayer, tasks_anisotropy, tasks_anisotropy_U


def analyze_interlayer_exchange(bilayer_interlayer_exchanges):
    signs = []
    for x in bilayer_interlayer_exchanges:
        signs.append(np.sign(x["exchange_energy"]))
        
    if all(sign == signs[0] for sign in signs):
        if len(signs) > 1:
            if signs[0] == -1:
                magnetic_variety = 'all FM'
            if signs[0] == 1:
                magnetic_variety = 'all AFM'
        else:
            magnetic_variety = 'None'
    if not all(sign == signs[0] for sign in signs):
        magnetic_variety = 'Mixed'

    return dict(magnetic_variety=magnetic_variety)


def analyze_anisotropies(monolayer_properties, bilayer_anisotropies):
    
    def compare_bilayers(bilayer_interlayer_exchanges, bilayer_anisotropies):
        
        def compare_anisotropy(bilayer_anisotropies):
            
            directions = []
            directions_U3 = []
            for bilayer in bilayer_anisotropies:
                dE_zx = bilayer["bilayer_anisotropy_zx"]
                dE_zy = bilayer["bilayer_anisotropy_zy"]
                if max(dE_zx, dE_zy) > 0:
                    directions.append(-1)
                else:
                    directions.append(1)

                dE_zx_U3 = bilayer["bilayer_anisotropy_zx_U3"]
                dE_zy_U3 = bilayer["bilayer_anisotropy_zy_U3"]
                if max(dE_zx_U3, dE_zy_U3) > 0:
                    directions_U3.append(-1)
                else:
                    directions_U3.append(1)

            if all(direction == directions[0] for direction in directions):
                if len(directions) > 1:
                    if directions[0] == -1:
                        anisotropic_variety = 'all in plane'
                    if directions[0] == 1:
                        anisotropic_variety = 'all out of plane'
                else:
                    anisotropic_variety = 'None'
            if not all(direction == directions[0] for direction in directions):
                anisotropic_variety = 'Mixed'
                
            if all(direction_U3 == directions_U3[0] for direction_U3 in directions_U3):
                if len(directions_U3) > 1:
                    if directions_U3[0] == -1:
                        anisotropic_variety_U3 = 'all in plane'
                    if directions_U3[0] == 1:
                        anisotropic_variety_U3 = 'all out of plane'
                else:
                    anisotropic_variety_U3 = 'None'
            if not all(direction_U3 == directions_U3[0] for direction_U3 in directions_U3):
                anisotropic_variety_U3 = 'Mixed'

            return dict(anisotropic_variety=anisotropic_variety, 
                        anisotropic_variety_U3=anisotropic_variety_U3)
        
        #magnetic_variety = compare_bilayers(bilayer_properties)
        anisotropic_variety = compare_anisotropy(bilayer_anisotropies)

        return dict(anisotropic_variety=anisotropic_variety)
        #return dict(magnetic_variety=magnetic_variety, anisotropic_variety=anisotropic_variety)

    def compare_monolayer_to_bilayer(monolayer_properties, bilayer_anisotropies):
        
        def compare_anisotropies(monolayer_properties, bilayer_anisotropies):

            monolayer_to_bilayer = []
            monolayer_to_bilayer_U3 = []

            mono_dE_zx = monolayer_properties["monolayer_anisotropy_zx"]
            mono_dE_zy = monolayer_properties["monolayer_anisotropy_zy"]
            mono_dE_zx_U3 = monolayer_properties["monolayer_anisotropy_zx_U3"]
            mono_dE_zy_U3 = monolayer_properties["monolayer_anisotropy_zy_U3"]
            
            bil_dE_zx = bilayer_anisotropies["bilayer_anisotropy_zx"]
            bil_dE_zy = bilayer_anisotropies["bilayer_anisotropy_zy"]
            bil_dE_zx_U3 = bilayer_anisotropies["bilayer_anisotropy_zx_U3"]
            bil_dE_zy_U3 = bilayer_anisotropies["bilayer_anisotropy_zy_U3"]
            
            if max(mono_dE_zx, mono_dE_zy) > 0 and max(bil_dE_zx, bil_dE_zy) > 0:
                monolayer_to_bilayer = 'in plane to in plane'
            if max(mono_dE_zx, mono_dE_zy) > 0 and max(bil_dE_zx, bil_dE_zy) < 0:
                monolayer_to_bilayer = 'in plane to out of plane'
            if max(mono_dE_zx, mono_dE_zy) < 0 and max(bil_dE_zx, bil_dE_zy) > 0:
                monolayer_to_bilayer = 'out of plane to in plane'
            if max(mono_dE_zx, mono_dE_zy) < 0 and max(bil_dE_zx, bil_dE_zy) < 0:
                monolayer_to_bilayer = 'out of plane to out of plane'
                    
            if max(mono_dE_zx_U3, mono_dE_zy_U3) > 0 and max(bil_dE_zx_U3, bil_dE_zy_U3) > 0:
                monolayer_to_bilayer_U3 = 'in plane to in plane'
            if max(mono_dE_zx_U3, mono_dE_zy_U3) > 0 and max(bil_dE_zx_U3, bil_dE_zy_U3) < 0:
                monolayer_to_bilayer_U3 = 'in plane to out of plane'
            if max(mono_dE_zx_U3, mono_dE_zy_U3) < 0 and max(bil_dE_zx_U3, bil_dE_zy_U3) > 0:
                monolayer_to_bilayer_U3 = 'out of plane to in plane'
            if max(mono_dE_zx_U3, mono_dE_zy_U3) < 0 and max(bil_dE_zx_U3, bil_dE_zy_U3) < 0:
                monolayer_to_bilayer_U3 = 'out of plane to out of plane'

            return monolayer_to_bilayer, monolayer_to_bilayer_U3

        compare_monolayer_to_bilayer = []
        compare_monolayer_to_bilayer_U3 = []
        for bilayer_anisotropy in bilayer_anisotropies:
            anisotropy_mono_to_bilayer, anisotropy_mono_to_bilayer_U3 = compare_anisotropies(monolayer_properties, bilayer_anisotropy)
            compare_monolayer_to_bilayer.append(anisotropy_mono_to_bilayer)
            compare_monolayer_to_bilayer_U3.append(anisotropy_mono_to_bilayer_U3)

        return dict(compare_monolayer_to_bilayer=compare_monolayer_to_bilayer, compare_monolayer_to_bilayer_U3=compare_monolayer_to_bilayer_U3)

    comparison_bilayer_results = compare_bilayers(bilayer_interlayer_exchanges, bilayer_anisotropies)
    comparison_monolayer_to_bilayer_results = compare_monolayer_to_bilayer(monolayer_properties, bilayer_anisotropies)
    
    return dict(comparison_bilayer_results=comparison_bilayer_results, comparison_monolayer_to_bilayer_results=comparison_monolayer_to_bilayer_results)

def get_properties(monolayer_folder, bilayer_folders):    

    def get_monolayer_properties(monolayer_folder):
        dpath = f"{monolayer_folder}/results-asr.magnetic_anisotropy_U.json"
        if os.path.exists(dpath):        
            atoms = read(f'{monolayer_folder}/structure.json')
    
            path = "/".join([x for x in monolayer_folder.split("/") if x != ""][-3:])

            try:
                monolayer_anisotropy_data = read_json("/home/niflheim2/cmr/C2DB-ASR/tree/" + path + "/" + "results-asr.magnetic_anisotropy.json")
            except:
                monolayer_anisotropy_data = read_json("/home/niflheim2/cmr/C2DB-ASR/icsd_cod_materials/tree/" + path + "/" + "results-asr.magnetic_anisotropy.json")

            monolayer_anisotropy_U3_data = read_json(f"{monolayer_folder}/results-asr.magnetic_anisotropy_U.json")
            monolayer_anisotropy_zx = monolayer_anisotropy_data["dE_zx"]
            monolayer_anisotropy_zy = monolayer_anisotropy_data["dE_zy"]
            monolayer_anisotropy_zx_U3 = monolayer_anisotropy_U3_data["dE_zx"]
            monolayer_anisotropy_zy_U3 = monolayer_anisotropy_U3_data["dE_zy"]

            return dict(monolayer_anisotropy_zx=monolayer_anisotropy_zx, monolayer_anisotropy_zy=monolayer_anisotropy_zy,
                        monolayer_anisotropy_zx_U3=monolayer_anisotropy_zx_U3, monolayer_anisotropy_zy_U3=monolayer_anisotropy_zy_U3)
        else:
            return []

    def get_bilayer_properties(bilayer):
        
        def get_interlayer_exchange(bilayer_folder):
            atoms = read(f"{bilayer_folder}/structure.json")
            dpath = f"{bilayer_folder}/results-asr.interlayer_magnetic_exchange.json"
            magnetic_atoms = []

            TM3d_atoms = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
               
            if os.path.exists(dpath):
                for atom in atoms:
                    if atom.symbol in TM3d_atoms:
                        magnetic_atoms.append(1)
                    else:
                        magnetic_atoms.append(0)
                
                mag_atoms = [x for x, z in enumerate(magnetic_atoms) if z == 1]
                data = read_json(dpath)

                magmoms_FM = data["magmoms_fm"]
                magmom_FM = data["M_FM"]
                magmoms_AFM = data["magmoms_afm"]
                magmom_AFM = data["M_AFM"]
                exchange_energy = data["eDIFF"]
 
                data = read_json(f"{bilayer_folder}/results-asr.bilayersummary.json")
                exfoliation_energy = data["exfoliation_energy"]
            
                gap = read_json(f"{bilayer_folder}/results-asr.gs.json")["gap"]
                
                deviation_matrix_fm = []
                for x in mag_atoms:
                    deviation_matrix_fm.append([ (abs(magmoms_FM[x]) - abs(magmoms_FM[y]))/abs(magmoms_FM[x]) for y in mag_atoms])
                deviation_matrix_fm = np.array(deviation_matrix_fm)
        
                deviation_matrix_afm = []
                for x in mag_atoms:
                    deviation_matrix_afm.append([ (abs(magmoms_AFM[x]) - abs(magmoms_AFM[y]))/abs(magmoms_AFM[x]) for y in mag_atoms])
                deviation_matrix_afm = np.array(deviation_matrix_afm)
            
                if not np.allclose((magmom_AFM/len(mag_atoms)),0, atol=0.01):
                    deviation_matrix_fm = 'None'
                    deviation_matrix_afm = 'None'
            else:
                deviation_matrix_fm = 'None'
                deviation_matrix_afm = 'None'

            check_values = []
            if not deviation_matrix_fm == 'None':
                if not deviation_matrix_afm == 'None':
                    for x in deviation_matrix_fm, deviation_matrix_afm:
                        for y in x:
                            for z in y:
                                if abs(z) > 0.05:
                                    check_values.append(z)
                                    print(bilayer)
                            
            if deviation_matrix_fm == 'None' or deviation_matrix_afm == 'None':
                check_values.append('None')
            if len(check_values) == 0:
                return dict(exchange_energy=exchange_energy, exfoliation_energy=exfoliation_energy,gap=gap,
                            magmoms_FM=magmoms_FM, magmom_FM=magmom_FM, magmoms_AFM=magmoms_AFM, magmom_AFM=magmom_AFM)

            if not len(check_values) == 0:
                return []

        def get_anisotropies(bilayer_folder):
            dpath = f"{bilayer_folder}/results-asr.magnetic_anisotropy_U.json"            
            if os.path.exists(dpath):
                bilayer_anisotropy_data = read_json(f"{bilayer_folder}/results-asr.magnetic_anisotropy.json")
                bilayer_anisotropy_zx = bilayer_anisotropy_data["dE_zx"]
                bilayer_anisotropy_zy = bilayer_anisotropy_data["dE_zy"]

                bilayer_anisotropy_U3_data = read_json(f"{bilayer_folder}/results-asr.magnetic_anisotropy_U.json")
                bilayer_anisotropy_zx_U3 = bilayer_anisotropy_U3_data["dE_zx"]
                bilayer_anisotropy_zy_U3 = bilayer_anisotropy_U3_data["dE_zy"]
            
                return dict(bilayer_anisotropy_zx=bilayer_anisotropy_zx, bilayer_anisotropy_zy=bilayer_anisotropy_zy, 
                            bilayer_anisotropy_zx_U3=bilayer_anisotropy_zx_U3, bilayer_anisotropy_zy_U3=bilayer_anisotropy_zy_U3)
            else:
                return []
        
        interlayer_exchange = get_interlayer_exchange(bilayer)
        anisotropies = get_anisotropies(bilayer)
        
        return interlayer_exchange, anisotropies

    bilayer_exchange_constants = []
    bilayer_anisotropies = []
    monolayer_properties = get_monolayer_properties(monolayer_folder)

    for bilayer in bilayer_folders:
        bilayer_interlayer_exchange, bilayer_anisotropy = get_bilayer_properties(bilayer)

        if not len(bilayer_interlayer_exchange) == 0:
            bilayer_exchange_constants.append(bilayer_interlayer_exchange)
        if not len(bilayer_anisotropy) == 0:
            bilayer_anisotropies.append(bilayer_anisotropy)

    return monolayer_properties, bilayer_exchange_constants, bilayer_anisotropies

def make_plots(monolayer_dE_zx_distribution, monolayer_dE_zy_distribution, 
               bilayer_dE_zx_distribution, bilayer_dE_zy_distribution,
               metal_bilayer_exchange_vs_exfoliation_energy, insulator_bilayer_exchange_vs_exfoliation_energy):
    
    plt.figure()
    plt.scatter([xy[0] for xy in monolayer_dE_zx_distribution[:]], [xy[1] for xy in monolayer_dE_zx_distribution[:]], c='b', alpha=0.5)
    plt.axline((0, 0), slope=1, color="black", linestyle=(0, (5, 5)))
    plt.axhline(y=0, color='black')
    plt.axvline(x=0, color='black')
    plt.xlabel("monolayer magnetic anisotropy zx")
    plt.ylabel("monolayer magnetic anisotropy zx with U3")

    plt.tight_layout()
    plt.savefig('anisotropy_monolayer_zx_distribution.pdf')

    plt.figure()
    plt.scatter([xy[0] for xy in monolayer_dE_zy_distribution[:]], [xy[1] for xy in monolayer_dE_zy_distribution[:]], c='b', alpha=0.5)
    plt.axline((0, 0), slope=1, color="black", linestyle=(0, (5, 5)))
    plt.axhline(y=0, color='black')
    plt.axvline(x=0, color='black')
    plt.xlabel("monolayer magnetic anisotropy zy")
    plt.ylabel("monolayer magnetic anisotropy zy with U3")

    plt.tight_layout()
    plt.savefig('anisotropy_monolayer_zy_distribution.pdf')

    plt.figure()
    plt.scatter([xy[0] for xy in bilayer_dE_zx_distribution[:]], [xy[1] for xy in bilayer_dE_zx_distribution[:]], c='b', alpha=0.5)
    plt.axline((0, 0), slope=1, color="black", linestyle=(0, (5, 5)))
    plt.axhline(y=0, color='black')
    plt.axvline(x=0, color='black')
    plt.xlabel("bilayer magnetic anisotropy zx")
    plt.ylabel("bilayer magnetic anisotropy zx with U3")

    plt.tight_layout()
    plt.savefig('anisotropy_bilayer_zx_distribution.pdf')

    plt.figure()
    plt.scatter([xy[0] for xy in bilayer_dE_zy_distribution[:]], [xy[1] for xy in bilayer_dE_zy_distribution[:]], c='b', alpha=0.5)
    plt.axline((0, 0), slope=1, color="black", linestyle=(0, (5, 5)))
    plt.axhline(y=0, color='black')
    plt.axvline(x=0, color='black')
    plt.xlabel("bilayer magnetic anisotropy zy")
    plt.ylabel("bilayer magnetic anisotropy zy with U3")
    
    plt.tight_layout()
    plt.savefig('anisotropy_bilayer_zy_distribution.pdf')


    plt.figure()
    plt.scatter([xy[0] for xy in metal_bilayer_exchange_vs_exfoliation_energy[:]], [xy[1] for xy in metal_bilayer_exchange_vs_exfoliation_energy[:]], c='b', alpha=0.5)
    plt.scatter([xy[0] for xy in insulator_bilayer_exchange_vs_exfoliation_energy[:]], [xy[1] for xy in insulator_bilayer_exchange_vs_exfoliation_energy[:]], c='r', alpha=0.5)
    plt.axhline(y=0, color='black')
    plt.xlabel("exfoliation energy")
    plt.ylabel("interlayer exchange energy")

    plt.tight_layout()
    plt.savefig('interlayer_exchange_vs_exfoliation_energy.pdf')


def absolute_value(arr, val):
    a  = arr[ np.abs(arr - val/100.*arr.sum()).argmin() ]
    return a

def make_piecharts(bilayer_to_bilayer_magnetic_states, monolayer_to_bilayer_anisotropy,
                   monolayer_to_bilayer_anisotropy_U3, bilayer_to_bilayer_anisotropy, bilayer_to_bilayer_anisotropy_U3):

    all_FM = []
    all_AFM = []
    mixed = []
    
    for x in bilayer_to_bilayer_magnetic_states:
        if x == 'all FM':
            all_FM.append(1)
        if x == 'all AFM':
            all_AFM.append(1)
        if x == 'Mixed':
            mixed.append(1)

    stables = [len(all_FM), len(all_AFM), len(mixed)]
    stables = np.array(stables)
    labels = ["All FM", "All AFM", "Mixed"]
    
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=200, subplot_kw=dict(aspect="equal"))

    explode = [0.05] * len(stables)

    wedges, texts, autotexts = ax.pie(stables, autopct=lambda val: absolute_value(stables, val), shadow=False)
    for autotext in autotexts:
        autotext.set_color('white')
  
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.2*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, **kw)

    fig.savefig('piechart_bilayer_to_bilayer_magnetic_states.pdf')

    all_in_plane = []
    all_out_of_plane = []
    mixed = []
    
    for x in bilayer_to_bilayer_anisotropy:
        if x == 'all in plane':
            all_in_plane.append(1)
        if x == 'all out of plane':
            all_out_of_plane.append(1)
        if x == 'Mixed':
            mixed.append(1)

    #stables = [len(all_in_plane), len(all_out_of_plane), len(mixed)]
    #stables = np.array(stables)
    #labels = ["All in plane", "All out of plane", "Mixed"]
    stables = [len(all_in_plane), len(mixed), len(all_out_of_plane)]
    stables = np.array(stables)
    labels = ["All in plane", "Mixed", "All out of plane"]

    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=200, subplot_kw=dict(aspect="equal"))

    explode = [0.05] * len(stables)

    wedges, texts, autotexts = ax.pie(stables, autopct=lambda val: absolute_value(stables, val), shadow=False)
    for autotext in autotexts:
        autotext.set_color('white')
  
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.2*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, **kw)

    fig.savefig('piechart_anisotropy_bilayer_to_bilayer.pdf')

    all_in_plane = []
    all_out_of_plane = []
    mixed = []    

    for x in bilayer_to_bilayer_anisotropy_U3:
        if x == 'all in plane':
            all_in_plane.append(1)
        if x == 'all out of plane':
            all_out_of_plane.append(1)
        if x == 'Mixed':
            mixed.append(1)

    stables = [len(all_in_plane), len(all_out_of_plane), len(mixed)]
    stables = np.array(stables)
    labels = ["All in plane", "All out of plane", "Mixed"]
    
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=200, subplot_kw=dict(aspect="equal"))

    explode = [0.05] * len(stables)

    wedges, texts, autotexts = ax.pie(stables, autopct=lambda val: absolute_value(stables, val), shadow=False)
    for autotext in autotexts:
        autotext.set_color('white')
  
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.2*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, **kw)

    fig.savefig('piechart_anisotropy_bilayer_to_bilayer_U3.pdf')

    in_plane_to_in_plane = []
    in_plane_to_out_of_plane = []
    out_of_plane_to_in_plane = []
    out_of_plane_to_out_of_plane = []

    for x in monolayer_to_bilayer_anisotropy:
        for y in x:
            if y == 'in plane to in plane':
                in_plane_to_in_plane.append(1)
            if y == 'in plane to out of plane':
                in_plane_to_out_of_plane.append(1)
            if y == 'out of plane to in plane':
                out_of_plane_to_in_plane.append(1)
            if y == 'out of plane to out of plane':
                out_of_plane_to_out_of_plane.append(1)

    #stables = [len(in_plane_to_in_plane), len(in_plane_to_out_of_plane), len(out_of_plane_to_in_plane), len(out_of_plane_to_out_of_plane)]
    #stables = np.array(stables)
    #labels = ["in plane \n to in plane", "in plane to \n out of plane", "out of plane \n to in plane", "out of plane \n to out of plane"]
    stables = [len(in_plane_to_in_plane), len(in_plane_to_out_of_plane),  len(out_of_plane_to_out_of_plane) , len(out_of_plane_to_in_plane)]
    stables = np.array(stables)
    labels = ["in plane \n to in plane", "in plane to \n out of plane", "out of plane \n to out of plane", "out of plane \n to in plane"]

    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=200, subplot_kw=dict(aspect="equal"))

    explode = [0.05] * len(stables)

    wedges, texts, autotexts = ax.pie(stables, autopct=lambda val: absolute_value(stables, val), shadow=False)
    for autotext in autotexts:
        autotext.set_color('white')
  
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")
    
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.2*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, **kw)

    fig.savefig('piechart_anisotropy_monolayer_to_bilayer.pdf')

    in_plane_to_in_plane = []
    in_plane_to_out_of_plane = []
    out_of_plane_to_in_plane = []
    out_of_plane_to_out_of_plane = []
    
    for x in monolayer_to_bilayer_anisotropy_U3:
        for y in x:
            if y == 'in plane to in plane':
                in_plane_to_in_plane.append(1)
            if y == 'in plane to out of plane':
                in_plane_to_out_of_plane.append(1)
            if y == 'out of plane to in plane':
                out_of_plane_to_in_plane.append(1)
            if y == 'out of plane to out of plane':
                out_of_plane_to_out_of_plane.append(1)

    stables = [len(in_plane_to_in_plane), len(in_plane_to_out_of_plane), len(out_of_plane_to_in_plane), len(out_of_plane_to_out_of_plane)]
    stables = np.array(stables)
    labels = ["in plane to \n in plane", "in plane to \n out of plane", "out of plane \n to in plane", "out of plane \n to out of plane"]
    
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=200, subplot_kw=dict(aspect="equal"))

    explode = [0.05] * len(stables)

    wedges, texts, autotexts = ax.pie(stables, autopct=lambda val: absolute_value(stables, val), shadow=False)
    for autotext in autotexts:
        autotext.set_color('white')
  
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.2*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, **kw)

    fig.savefig('piechart_anisotropy_monolayer_to_bilayer_U3.pdf')

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*")
    args = parser.parse_args()
   
    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]

    monolayer_folders = filter_monolayers(args.folders)

    magnetic_bilayers = []
    magnetic_monolayers = []
    select_bilayers = []
    select_ani_bilayers = []

    monolayer_dE_zy_distribution = []
    monolayer_dE_zx_distribution = []
    bilayer_dE_zy_distribution = []
    bilayer_dE_zx_distribution = []

    monolayer_to_bilayer_anisotropy = []
    monolayer_to_bilayer_anisotropy_U3 = []
    bilayer_to_bilayer_anisotropy = []
    bilayer_to_bilayer_anisotropy_U3 = []

    bilayer_to_bilayer_magnetic_states = []
    #bilayer_exchange_vs_exfoliation_energy = []
    metal_bilayer_exchange_vs_exfoliation_energy = []
    insulator_bilayer_exchange_vs_exfoliation_energy = []

    tasks1 = []
    tasks2 = []
    tasks3 = []
    
    for monolayer_folder in monolayer_folders:
        selected_bilayers = accepted_bilayers(monolayer_folder, only_most_stable=False)    
        magnetic_bilayers.append(selected_bilayers)

        if not len(selected_bilayers) == 0:
            tasks_interlayer, tasks_anisotropy, tasks_anisotropy_U = get_failed_calculations(monolayer_folder, selected_bilayers)

            if not len(tasks_interlayer) == 0:
                tasks1.append(len(tasks_interlayer))
            if not len(tasks_anisotropy) == 0:
                tasks2.append(len(tasks_anisotropy))
            if not len(tasks_anisotropy_U) == 0:
                tasks3.append(len(tasks_anisotropy_U))

            monolayer_properties, bilayer_interlayer_exchanges, bilayer_anisotropies = get_properties(monolayer_folder, selected_bilayers)
            select_bilayers.append(len(bilayer_interlayer_exchanges))
            select_ani_bilayers.append(len(bilayer_anisotropies))
            if not len(bilayer_interlayer_exchanges) == 0:
                comparison_exchange = analyze_interlayer_exchange(bilayer_interlayer_exchanges)
                bilayer_to_bilayer_magnetic_states.append(comparison_exchange["magnetic_variety"])
                for bilayer_interlayer_exchange in bilayer_interlayer_exchanges:
                    if bilayer_interlayer_exchange["gap"] == 0:
                        metal_bilayer_exchange_vs_exfoliation_energy.append([bilayer_interlayer_exchange["exfoliation_energy"], bilayer_interlayer_exchange["exchange_energy"]])
                    if not bilayer_interlayer_exchange["gap"] == 0:
                        insulator_bilayer_exchange_vs_exfoliation_energy.append([bilayer_interlayer_exchange["exfoliation_energy"], bilayer_interlayer_exchange["exchange_energy"]])

            if not len(bilayer_anisotropies) == 0:
                comparison_anisotropies = analyze_anisotropies(monolayer_properties, bilayer_anisotropies)

                monolayer_dE_zy_distribution.append([monolayer_properties["monolayer_anisotropy_zy"], monolayer_properties["monolayer_anisotropy_zy_U3"]])
                monolayer_dE_zx_distribution.append([monolayer_properties["monolayer_anisotropy_zx"], monolayer_properties["monolayer_anisotropy_zx_U3"]])
                monolayer_to_bilayer_anisotropy.append(comparison_anisotropies["comparison_monolayer_to_bilayer_results"]["compare_monolayer_to_bilayer"])
                monolayer_to_bilayer_anisotropy_U3.append(comparison_anisotropies["comparison_monolayer_to_bilayer_results"]["compare_monolayer_to_bilayer_U3"])

                bilayer_to_bilayer_anisotropy.append(comparison_anisotropies["comparison_bilayer_results"]["anisotropic_variety"]["anisotropic_variety"])
                bilayer_to_bilayer_anisotropy_U3.append(comparison_anisotropies["comparison_bilayer_results"]["anisotropic_variety"]["anisotropic_variety_U3"])

                for bilayer_anisotropy in bilayer_anisotropies:
                    bilayer_dE_zy_distribution.append([bilayer_anisotropy["bilayer_anisotropy_zy"], bilayer_anisotropy["bilayer_anisotropy_zy_U3"]])
                    bilayer_dE_zx_distribution.append([bilayer_anisotropy["bilayer_anisotropy_zx"], bilayer_anisotropy["bilayer_anisotropy_zx_U3"]])
                    
    make_plots(monolayer_dE_zx_distribution, monolayer_dE_zy_distribution, 
               bilayer_dE_zx_distribution, bilayer_dE_zy_distribution,
               metal_bilayer_exchange_vs_exfoliation_energy, insulator_bilayer_exchange_vs_exfoliation_energy)

    make_piecharts(bilayer_to_bilayer_magnetic_states, monolayer_to_bilayer_anisotropy,
                   monolayer_to_bilayer_anisotropy_U3, bilayer_to_bilayer_anisotropy, bilayer_to_bilayer_anisotropy_U3)

    print(sum(select_bilayers))
    print(sum(select_ani_bilayers))
    print('interlayer tasks', sum(tasks1))
    print('anisotropy tasks', sum(tasks2))
    print('anisotropy U tasks', sum(tasks3))
