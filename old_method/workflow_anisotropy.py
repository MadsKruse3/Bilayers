from myqueue.task import task as mqtask
from pathlib import Path
from asr.core import read_json
from ase.io import read
import numpy as np
import os, re

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


def getresources(folder):
    atom = read(f'{folder}/structure.json')
    stochiometry = len(atom)

    if stochiometry < 5:
        resources = "16:2d"
        
    if 4 < stochiometry < 8:
        resources = "24:2d"
        
    if 7 < stochiometry < 12:
        resources = "40:2d"

    if stochiometry > 11:
        resources = "40:2d"
        
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
        res = getresources(folder)
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

def stability_criterion_monolayer(folder):
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

    if not len(selected) == 0:
        t = selected[0]
        t = t[0]
        return t
    else:
        return []

def basic_properties(bilayer_folder):
    res = getresources(bilayer_folder)
    tasks = []

    tasks += [task("asr.gs", resources=res, folder=bilayer_folder,
                   restart=2)]
    tasks += [task("asr.magstate", resources="1:10m", folder=bilayer_folder,
                   deps=["asr.gs"])]
    return tasks


def magnetic_anisotropy_U(folder):
    """Submit asr.magnetic_anisotropy_U based on criteria

    Only submit if asr.interlayer_magnetic_exchange_exists.
    """
    tasks = []


    tasks += [task("asr.magnetic_anisotropy_U", resources="1:20m", folder=folder,
                   restart=2, creates="results-asr.magnetic_anisotropy_U.json")]
    return tasks

def interlayer_magnetic_exchange_new(folder):
    """Submit asr.interlayer_magnetic_exchange_new in selected folders.
    
    """
    tasks = []


    tasks += [task("asr.interlayer_magnetic_exchange_new -u 3", resources="40:2d", folder=folder,
                       restart=2, creates="results-asr.interlayer_magnetic_exchange.json")]
    return tasks



def magnetic_anisotropy(folder):
    """Submit asr.magnetic_anisotropy based on criteria

    Only submit if asr.interlayer_magnetic_exchange_exists.
    """
    tasks = []

    tasks += [task("asr.magnetic_anisotropy", resources="1:10m", folder=folder,
                       restart=2, creates="results-asr.magnetic_anisotropy.json")]
    return tasks

def magnetic_anisotropy_U_monolayer(monolayer_folder):
    """Submit asr.magnetic_anisotropy_U based on criteria

    Only submit if asr.interlayer_magnetic_exchange_exists.
    """
    tasks = []

    tasks += [task("asr.magnetic_anisotropy_U", resources="1:10m", folder=monolayer_folder,
                       restart=2, creates="results-asr.magnetic_anisotropy_U.json")]
    
    return tasks


#def magnetic_anisotropy_U(folder):
#    """Submit asr.magnetic_anisotropy_U based on criteria
#
#    Only submit if asr.interlayer_magnetic_exchange_exists.
#    """
#    tasks = []

#    if Path(f"{folder}/results-asr.interlayer_magnetic_exchange.json").is_file():
#        tasks += [task("asr.magnetic_anisotropy_U", resources="1:20m", folder=folder,
#                      restart=2, creates="results-asr.magnetic_anisotropy_U.json")]
#        return tasks
#    else:
#        return tasks

#def interlayer_magnetic_exchange_new(folder):
#    """Submit asr.interlayer_magnetic_exchange_new in selected folders.
    
#    """
#    tasks = []

#    tasks += [task("asr.interlayer_magnetic_exchange_new -u 3", resources="40:2d", folder=folder,
#                       restart=2, creates="results-asr.interlayer_magnetic_exchange.json")]
#    return tasks



#def magnetic_anisotropy(folder):
#    """Submit asr.magnetic_anisotropy based on criteria

#    Only submit if asr.interlayer_magnetic_exchange_exists.
#    """
#    tasks = []

#    if Path(f"{folder}/results-asr.interlayer_magnetic_exchange.json").is_file():
#        tasks += [task("asr.magnetic_anisotropy", resources="1:10m", folder=folder,
#                       restart=2, creates="results-asr.magnetic_anisotropy.json")]
#        return tasks
#    else:
#        return tasks

#def magnetic_anisotropy_U_monolayer(monolayer_folder, bilayer_folder):
#    """Submit asr.magnetic_anisotropy_U based on criteria
#
#    Only submit if asr.interlayer_magnetic_exchange_exists.
#    """
#    tasks = []
#    #magnetic_bilayer_present = []

#    bilayer_folder = f'{bilayer_folder}'
#    if Path(f"{monolayer_folder}/results-asr.monolayer_magnetism.json").is_file():
#        
#        tasks += [task("asr.magnetic_anisotropy_U", resources="1:10m", folder=monolayer_folder,
#                       restart=2, creates="results-asr.magnetic_anisotropy_U.json")]
#        return tasks
#    else:
#        return tasks

def magnetic_monolayer(monolayer_folder):
    """Submit asr.bilayer_anisotropy based on criteria.

    Only submit if material contains one or more of the follow
    atoms:
    TM3d_atoms = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
    Only submit if asr.interlayer_magnetic_exchange exists.
    """
    
    tasks = []

    #if Path(f"{bilayer_folder}/results-asr.interlayer_magnetic_exchange.json").is_file():
    #res = getresources(monolayer_folder)
    tasks += [task("asr.monolayer_magnetism -u 3", resources='24:2d', folder=monolayer_folder,
                       restart=2, creates="results-asr.monolayer_magnetism.json")]
    
    
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
    if Path(f"{bilayer_folder}/results-asr.magstate.json").is_file():
        magdata = read_json(f"{bilayer_folder}/results-asr.magstate.json")
        if not magdata["is_magnetic"]:
            return []

    if Path(f"{bilayer_folder}/results-asr.interlayer_magnetic_exchange_U0.json").is_file():
        return []
    
    return [task("asr.interlayer_magnetic_exchange_U0", resources='48:2d', folder=bilayer_folder,
                 restart=2, creates="results-asr.interlayer_magnetic_exchange_U0.json")]
        

def emasses(bilayer_folder):
    """Run emasses if material has a gap"""
    gap = read_json(f"{bilayer_folder}/results-asr.gs.json").get("gap")
    if gap > 0.01:
        res = getresources(bilayer_folder)
        return [task("asr.emasses", resources=res,
                     folder=bilayer_folder, restart=2)]
    else:
        return []


def bandstructure(bilayer_folder):
    res = getresources(bilayer_folder)
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

    #tasks += vdw(folder) + create_bilayers(folder)
    #if not all_done(tasks):
    #    return tasks

    #tasks += relax_bilayers(folder)
    #if not all_done(tasks):
    #    return tasks

    #selected_monolayer = stability_criterion_monolayer(folder)
    #monolayer_folder = re.split('[/]',f'{selected_monolayer}')
    #monolayer_folder.remove(monolayer_folder[-1])
    #monolayer_folder = "/".join(monolayer_folder)

    selected = stability_criterion(folder)
    
    for bilayer_folder in selected:
        tasks += magnetism(bilayer_folder)
   
        #tasks += magnetic_anisotropy(bilayer_folder)
        #tasks += magnetic_anisotropy_U_monolayer(monolayer_folder, bilayer_folder)
        #tasks += magnetic_anisotropy_monolayer(monolayer_folder, bilayer_folder)
        #if len(tasks) == 0:
        #tasks += magnetic_monolayer(monolayer_folder, bilayer_folder)
 
            
    #tasks = []
    #wrong_states = np.loadtxt('/home/niflheim2/cmr/WIP/stacking/wrong_initial_states.npy', dtype=np.str)
    #for bilayer_folder in wrong_states:
    #    tasks += interlayer_magnetic_exchange_new(bilayer_folder)


    #tasks = []
    #monolayers = np.loadtxt('/home/niflheim2/cmr/WIP/stacking/monolayer_folders.npy', dtype=np.str)
    #for monolayer_folder in monolayers:
    #    tasks += magnetic_monolayer(monolayer_folder)
    
    

    #if not all_done(tasks):
    #    return tasks
    
    #tasks = []
    #monolayers_anisotropy_U = np.loadtxt('/home/niflheim2/cmr/WIP/stacking/monolayer_anisotropy_U.npy', dtype=np.str)
    #for monolayer_folder in monolayers_anisotropy_U:
    #    tasks += magnetic_anisotropy_U_monolayer(monolayer_folder)

    #monolayers_anisotropy = np.loadtxt('/home/niflheim2/cmr/WIP/stacking/monolayer_anisotropy.npy', dtype=np.str)
    #for monolayer_folder in monolayers_anisotropy:
    #    tasks += magnetic_anisotropy(monolayer_folder)

    #if not all_done(tasks):
    #    return tasks

    #tasks = []
    #bilayers = np.loadtxt('/home/niflheim2/cmr/WIP/stacking/bilayer_anisotropy.npy', dtype=np.str)
    #for monolayer_folder in bilayers:
    #    tasks += magnetic_anisotropy(monolayer_folder)

    #tasks = []
    #bilayers_U = np.loadtxt('/home/niflheim2/cmr/WIP/stacking/bilayer_anisotropy_U.npy', dtype=np.str)
    #for monolayer_folder in bilayers_U:
    #    tasks += magnetic_anisotropy_U(monolayer_folder)

    #if not all_done(tasks):
    #    return tasks
    
    
    return tasks

if __name__ == "__main__":
    tasks = create_tasks()

    for ptask in tasks:
        print(ptask, ptask.is_done())
