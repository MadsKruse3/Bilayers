from pathlib import Path
from asr.core import read_json
import numpy as np
import sys


p = Path("/home/niflheim2/cmr/C2DB-ASR/tree")
folders = sys.argv[1:]
fname = "stablematerials.npy"

stables = []
i = 0

for z in p.glob("*"):
    if not any(f == z.absolute().name and z.is_dir() for f in folders):
        continue

    for x in z.glob("*/*"):
        i += 1
        print(f"Analysing material {i}", flush=True, end="\r")
        stable1 = False
        stable2 = False
        stable3 = False
        afm     = False
        for y in x.glob("results-asr.*.json"):
            if "convex_hull.json" in str(y):
                data = read_json(str(y))
                stable1 = data["thermodynamic_stability_level"] == 3
                if not stable1:
                    break
            elif "phonons.json" in str(y):
                data = read_json(str(y))
                stable2 = data["dynamic_stability_phonons"] == "high"
                if not stable2:
                    break
            elif "stiffness.json" in str(y):
                data = read_json(str(y))
                stable3 = data["dynamic_stability_stiffness"] == "high"
                if not stable3:
                    break
            elif "magstate.json" in str(y):
                data = read_json(str(y))
                afm = data["magstate"].lower() == "afm"
                if afm:
                    break
                    
        if not (stable1 and stable2 and stable3) or afm:
            continue

        stables.append(str(x))


print(f"Analysing material {i}", flush=True)
print(f"Number of stable materials: {len(stables)}")
np.save(fname, stables)
