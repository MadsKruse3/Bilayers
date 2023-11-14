from myqueue.task import task as mqtask
from pathlib import Path
from asr.core import read_json
from ase.io import read
import numpy as np
import os

### UTILS ###


def listdirs(path):
    if type(path) != str:
        path = str(path.absolute())
    return [f"{path}/" + x for x in os.listdir(path) if os.path.isdir(f"{path}/{x}")]


def task(*args, **kwargs):
    """Get MyQueue task instance."""
    name = kwargs.get("name") or args[0]
    if "creates" not in kwargs:
        kwargs["creates"] = [f"results-{name}.json"]
    return mqtask(*args, **kwargs)


def getresources():
    """Randomly generate a resource.

    Use this to distribute jobs evenly across multiple
    partitions.
    """
    r = np.random.rand()
    if r < 0.6:
        resources = "40:5h"
    elif r < 0.8:
        resources = "16:1d"
    else:
        resources = "24:10h"

    return resources


### /UTILS ###


"""
Workflow:
0. Calculate vdw correction to monolayer energy
1. Create bilayers
2. Relax bilayers
3. Select relevant bilayers
4. Calculate GS
5. Calculate BS
6a. Calculate emasses for insulators
7a. Calculate magnetic state via custom bilayer mag recipe.
"""


def vdw(folder):
    """Calculate vdw correction.

    Requires that vdwtask.py is in PYTHONPATH.
    """
    if os.path.exists(f"{folder}/vdw_e.npy"):
        return []

    return [task("vdwtask", resources="1:10m", folder=folder)]


def create_bilayers(folder):
    tasks = [task("asr.stack_bilayer", resources="1:30m", folder=folder,
                  restart=2)]
    return tasks


def relax_bilayers(folder):
    restart = 3

    tasks = []
    for subfolder in listdirs(folder):
        res = getresources()
        tasks += [task("asr.relax_bilayer", resources=res,
                       restart=restart, folder=subfolder)]
    return tasks


def stability_criterion(folder):
    """Apply stability criterion."""
    from collectlayerdata import get_energy_length
    energies = []
    for subf in listdirs(folder):
        resultsfname = f"{subf}/results-asr.relax_bilayer.json"
        if not os.path.exists(resultsfname):
            continue

        energy, length = get_energy_length(subf)

        if energy is None:
            continue

        energies.append((subf, energy))

    if len(energies) == 0:
        return []

    nmats = 5
    deltaE = 0.002
    cutoff = 0.15
    maxE = max([e for s, e in energies])
    # Select materials within window of maxE
    selected = list(filter(lambda t: abs(t[1] - maxE) <= deltaE, energies))
    # Select materials that pass exfoliable criterion
    selected = list(filter(lambda t: t[1] <= cutoff, selected))
    # Select a max number of materials
    selected = sorted(selected, key=lambda t: t[1])[:nmats]
    return [t[0] for t in selected]


def basic_properties(bilayer_folder):
    res = getresources()
    tasks = []

    tasks += [task("asr.gs", resources=res, folder=bilayer_folder,
                   restart=2)]
    tasks += [task("asr.magstate", resources="1:10m", folder=bilayer_folder,
                   deps=["asr.gs"])]
    return tasks


def magnetism(bilayer_folder):
    """Submit asr.bilayer_magnetism based on criteria.

    Only submit if material contains one or more of the follow
    atoms:
    TM3d_atoms = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
    Only submit if magstate is magnetic
    """
    TM3d_atoms = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
    # Relax should be done here so we can read structure
    atoms = read(f"{bilayer_folder}/structure.json")
    if not any(x in atoms.symbols for x in TM3d_atoms):
        return []
    # magstate should also be done
    magdata = read_json(f"{bilayer_folder}/results-asr.magstate.json")
    if not magdata["is_magnetic"]:
        return []

    # Use 3 eV for all atoms
    res = getresources()
    return [task("asr.bilayer_magnetism -u 3", resources=res, folder=bilayer_folder,
                 restart=2)]


def emasses(bilayer_folder):
    """Run emasses if material has a gap"""
    gap = read_json(f"{bilayer_folder}/results-asr.gs.json").get("gap")
    if gap > 0.01:
        res = getresources()
        return [task("asr.emasses", resources=res,
                     folder=bilayer_folder, restart=2)]
    else:
        return []


def bandstructure(bilayer_folder):
    res = getresources()
    return [task("asr.bandstructure", resources=res,
                 folder=bilayer_folder, restart=2)]


def all_done(list_of_tasks):
    """Determine if all tasks in list_of_tasks are done."""
    return all([task.is_done() for task in list_of_tasks])


def create_tasks():
    folder = Path(".").absolute()
    atoms = read(f"{folder}/structure.json")
    if len(atoms) > 8:
        return []

    tasks = []

    tasks += vdw(folder) + create_bilayers(folder)
    if not all_done(tasks):
        return tasks

    tasks += relax_bilayers(folder)
    if not all_done(tasks):
        return tasks

    selected = stability_criterion(folder)
    for bilayer_folder in selected:
        tasks += basic_properties(bilayer_folder)
    if not all_done(tasks):
        return tasks

    for bilayer_folder in selected:
        tasks += emasses(bilayer_folder)
        tasks += magnetism(bilayer_folder)
    if not all_done(tasks):
        return tasks

    for bilayer_folder in selected:
        tasks += bandstructure(bilayer_folder)

    return tasks


if __name__ == "__main__":
    tasks = create_tasks()

    for ptask in tasks:
        print(ptask, ptask.is_done())
