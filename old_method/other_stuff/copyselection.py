import numpy as np
import os
import sys
from pathlib import Path
from shutil import copyfile
"""
This is a utility script that copies selected materials.

These materials are assumed to be stored in a list in the
format that selectstablematerials.py produces.

The file containing the materials is given as an
argument to the script.

The materials are copied in such a way as to preserve
the folder structure in the original tree.

It is assumed to be run on materials selected from the
C2DB-ASR tree (with the structure it had in September 2020)
so if it is run on a tree with different structure it may fail.
"""

if len(sys.argv) > 1:
    filename = sys.argv[1] 
else:
    filename = "stablematerials.npy"

selection = np.load(filename)
    

def easymake(f):
    if not os.path.exists(f):
        os.mkdir(f)

easymake("tree")

for folder in selection:
    s, f, m = folder.split("/")[-3:]
    s = "tree/" + s
    f = s + "/" + f
    m = f + "/" + m
    
    easymake(s)
    easymake(f)
    easymake(m)

    copyfile(folder + "/structure.json", m + "/structure.json")

