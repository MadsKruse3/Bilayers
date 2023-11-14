"""
This script calls ase gui in a helpful way on a number
of structure-files.

Useful for looking at and comparing bilayers.
"""

import sys
from subprocess import Popen

folders = sys.argv[1:]


for folder in folders:
    name = folder.replace("/", "_")
    # First copy the file. ase gui uses the filename as the window title,
    # so we do this to avoid every window being called 'structure.json'.
    cmd = ["cp", "-v", f"{folder}/structure.json", f"{folder}/{name}_structure.json"]
    cmd2 = ["ase", "gui", f"{folder}/{name}_structure.json"]
    Popen(cmd)
    Popen(cmd2)



    
