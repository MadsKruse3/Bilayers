"""Workflow for C2DB.

Important utility functions in this script:

  create_tasks
  get_mag_states
  get_most_stable_magstate
  task
  is_thermodynamically_stable
  is_dynamically_stablex
  all_done

The workflow is factored into four components:

  magstate_workflow
  therm_stability_workflow
  dynamical_stability_workflow
  property_workflow

"""
from myqueue.task import task as mqtask
from pathlib import Path
from asr.core import read_json, chdir, write_json
from asr.setup.strains import main as setupstrains
from asr.setup.strains import get_strained_folder_name, get_relevant_strains
from ase.io import read
import numpy as np
import os

VERBOSE = os.environ.get('MQVERBOSE', False)

def listonlydirs(path):
    return [f"{path}/" + x for x in os.listdir(path) if os.path.isdir(f"{path}/{x}")]


def task(*args, **kwargs):
    """Get MyQueue task instance."""
    name = kwargs.get("name") or args[0]
    if "creates" not in kwargs:
        kwargs["creates"] = [f"results-{name}.json"]
    return mqtask(*args, **kwargs)


def get_cwd():
    """Get current working directory."""
    return Path('.').absolute()

def stacking_setup(folder):
    tasks = []
    if not "results-asr.stack_bilayer.json" in list(map(str, folder.iterdir())):
        tasks += [task("asr.stack_bilayer", resources="1:30m",
                       folder=folder)]
    return tasks


def vdw_task(folder):
    tasks = []
    if not "vdw_e.npy" in list(map(str, folder.iterdir())):
        tasks += [task("vdwtask", resources="1:10m",
                       folder=folder)]
    return tasks



def is_magnetic(folder):
    folder = str(folder.absolute())
    name_components = folder.split("/")
    name_w_id = ""
    for comp in name_components[::-1]:
        name_w_id = f"/{comp}" + name_w_id
        if comp == "tree":
            break

    c2db_name = f"/home/niflheim2/cmr/C2DB-ASR/{name_w_id}/results-asr.magstate.json"
    import os
    if os.path.exists(c2db_name):
        data = read_json(c2db_name)
        return data["is_magnetic"]
    else:
        magnetics = ["Fe", "Ni", "Co",]
        return any(x in folder for x in magnetics)


def add_settings(params):
    d = {"d3": True,
         "xc": "PBE",
         "PWE": 800,
         "kpts": {"density": 6.0, "gamma": True},
         "mixer": {"type": "default",
                   "beta": None,
                   "nmaxold": None,
                   "weight": None}}
    params["asr.relax_bilayer"]["settings"] = d
    

def stacking_workflow(folder):
    restart = 2

    tasks = []
    for subfolder in listonlydirs(folder):
        # resources = "48:4h" if np.random.rand() > 0.5 else "40:5h"
        r = np.random.rand()
        if r < 0.45:
            resources = "40:xeon40_clx:5h"
        elif r < 0.8:
            resources = "48:5h"
        else:
            resources = "40:5h"

        # err_exists = any("asr.relax_bilayer." in str(x) and ".err" in str(x) for x in Path(subfolder).iterdir())
        if os.path.exists(str(folder) + f"/{subfolder}/asr.relax_bilayer.FAILED"):
            paramname = f"{str(folder)}/{subfolder}/params.json"
            if os.path.exists(paramname):
                params = read_json(paramname)
            else:
                params = {}
                params["asr.relax_bilayer"] = {}
            params["asr.relax_bilayer"]["restart"] = True
            if "settings" not in params["asr.relax_bilayer"]:
                add_settings(params)
            if params["asr.relax_bilayer"]["settings"]["mixer"] == "mixerdif":
                params["asr.relax_bilayer"]["settings"]["mixer"] = {"type": "mixerdif",
                                                                    "beta":0.02,
                                                                    "nmaxold": 5,
                                                                    "weight": 100}
            elif is_magnetic(folder):
                params["asr.relax_bilayer"]["settings"]["mixer"] = "mixerdif"
            else:
                params["asr.relax_bilayer"]["settings"]["mixer"] = {"type": "mixer",
                                                                    "beta":0.02,
                                                                    "nmaxold": 5,
                                                                    "weight": 100}

            write_json(f"{str(folder)}/{subfolder}/params.json", params)

        tasks += [task("asr.relax_bilayer", resources=resources,
                       restart=restart,
                       folder=str(folder) + "/" + subfolder)]
    return tasks


def basic_workflow(folder):
    """Generate tasks related finding most stable magnetic structures."""
    tasks = []
    # resources = "48:5h" if np.random.rand() > 0.2 else "40:5h" 
    r = np.random.rand()
    if r < 0.7:
        resources = "40:xeon40_clx:5h"
    elif r < 0.85:
        resources = "48:5h"
    else:
        resources = "40:5h"

    tasks += [task("asr.structureinfo", resources="1:10m", folder=folder)]
    tasks += [task("asr.gs@calculate", resources=resources, folder=folder)]
    tasks += [task("asr.gs", folder=folder, resources=resources,
                   deps=["asr.structureinfo", "asr.gs@calculate"]),
              task("asr.magstate", deps=["asr.gs"], folder=folder)]

    if not Path('results-asr.magstate.json').is_file():
        return tasks
    else:
        results = read_json("results-asr.magstate.json")

        ismag = results["is_magnetic"]
        if ismag:
            tasks += [task("asr.magnetic_anisotropy", resources=resources,
                   folder=folder, deps="asr.gs@calculate")]

    return tasks


def all_done(list_of_tasks):
    """Determine if all tasks in list_of_tasks are done."""
    return all([task.is_done() for task in list_of_tasks])


def magnetic_workflow(folder):
    tasks = []
    resources = "48:5h" if np.random.rand() > 0.2 else "40:5h" 
    tasks += []
    return tasks
    

def property_workflow(folder):
    """Generate tasks for various material properties."""
    verbose_print('Executing property_workflow()')
    # resources = "48:5h" if np.random.rand() > 0.2 else "40:5h" 
    r = np.random.rand()
    if r < 0.45:
        resources = "40:xeon40_clx:5h"
    elif r < 0.8:
        resources = "48:5h"
    else:
        resources = "40:5h"

    gres = lambda : "48:5h" if np.random.rand() > 0.2 else "40:5h" 
    atoms = read(folder / 'structure.json')
    tasks = [task("asr.bandstructure",
                  folder=folder,
                  resources=resources,
                  deps=["asr.gs"]),
             task("asr.projected_bandstructure",
                  folder=folder,
                  resources="1:10m",
                  deps=["asr.bandstructure"]),
             # task("asr.polarizability",
             #      folder=folder,
             #      diskspace=100,
             #      resources="120:20h",
             #      deps=["asr.gs"]),
             task('asr.pdos', folder=folder,
                  resources=resources,
                  deps=['asr.gs'])]

    gsresults = read_json(folder / "results-asr.gs.json")
    structureinfo = read_json(folder / "results-asr.structureinfo.json")
    gap = gsresults.get("gap")

    if gap > 0.01:
        verbose_print(get_cwd(), 'Has band gap.')
        tasks += [task("asr.emasses",
                       folder=folder,
                       resources=resources,
                       deps=["asr.gs"]),
                  task("asr.formalpolarization",
                       folder=folder,
                       resources=gres())]
    else:
        tasks += [
            task(
                "asr.plasmafrequency",
                folder=folder,
                resources="40:20h",
                diskspace=100,
                creates=["results-asr.plasmafrequency.json"],
                deps=["asr.gs"],
            )]

    return tasks


def return_tasks(tasks):
    """Wrap function for returning tasks."""
    if VERBOSE:
        print(get_cwd(), tasks)
    return tasks


def verbose_print(*args):
    """Only print if VERBOSE."""
    if VERBOSE:
        print(*args)


def stability_criterion(folder):
    """Apply stability criterion."""
    from collectlayerdata import get_energy_length
    energies = []
    subs = list(folder.iterdir())
    for subf in subs:
        if not subf.is_dir():
            continue
        resultsfname = f"{subf}/results-asr.relax_bilayer.json"
        if not os.path.exists(resultsfname):
            continue
        data = read_json(resultsfname)
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


def create_tasks():
    """Create MyQueue Task list for the Bilayer workflow.

    Note that this workflow relies on the folder layout so be
    careful.
    """
    verbose_print(get_cwd())
    folder = Path('.')
    tasks = []
    
    # Comment this part out if pymatgen is causing trouble.
    # Run stacking manually instead.
    tasks += stacking_setup(folder)
    tasks += vdw_task(folder)
    if not all_done(tasks):
        return tasks

    tasks += stacking_workflow(folder)
    if not all_done(tasks):
        return tasks

    selected = stability_criterion(folder)

    for subf in selected:
        subtasks = []
        if not subf.is_dir():
            continue
        subfolder = subf.absolute()
        subtasks += basic_workflow(subfolder)
        if not all_done(subtasks):
            tasks += subtasks
            continue

        subtasks += property_workflow(subfolder)
        tasks += subtasks

    return tasks


if __name__ == '__main__':
    tasks = create_tasks()

    for ptask in tasks:
        print(ptask, ptask.is_done())
