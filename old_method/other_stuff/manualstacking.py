import sys
from pathlib import Path
from subprocess import call
import os

p = Path(".")

if len(sys.argv) > 1:
    folders = sys.argv[1:]
    N = len(folders)
    for i, x in enumerate(folders):
        print(f"Stacking folder {i}/{N}")
        if "results-asr.stack_bilayer.json" in os.listdir(str(x)):
            print(f"Skipping {x}", flush=True)
            print("")
            continue

        call(["asr", "run", "stack_bilayer", str(x)])
        print("")
else:
    if "results-asr.stack_bilayer.json" in os.listdir("."):
        print(f"Skipping", flush=True)
        print("")
    else:
        call(["asr", "run", "stack_bilayer"])
        print("")
