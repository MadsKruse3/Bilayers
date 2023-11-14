# Bilayer Workflow

Tools and guides for performing the bilayer workflow.

To work with the bilayer workflow you must be on the ASR branch "stacking". You must also have myqueue installed.

## Getting started

0. Add the location of this repository to your PYTHONPATH.

1. Use the script selectstablematerials.py to get a list of highly stable materials from C2DB-ASR. For example "python3 \<PathTobilayerworkflow Repo>/selectstablematerials.py A AB ABC" will collect a list of stable materials from C2DB with stoichiometries A, AB, and ABC.

2. Copy the selected materials into your own folder structure using "python3 \<PathTobilayerworkflow Repo>/copyselection.py". copyselection.py creates a folder structure that mimics the C2DB folder structure, e.g. it will create a folder "tree" with subfolders "A", "AB", "ABC", etc.

3. Run the workflow with "mq workflow \<PathTobilayerworkflow Repo>/stackingworkflow.py tree\/\*\/\*\/\*\/".

    3.1 If pymatgen is causing trouble, create a venv using the gpaw-venv.sh script.

    3.2 If it still doesn't work, you can manually run the stacking part by using manualstacking.py



## Utilities
collectlayerdata.py contains the function collectdata which is a handy utility to collect binding energies and binding lengths of the bilayers and group them in a dict according to stoichiometry, monolayer chemical formula, or bilayer descriptor.

viewatoms.py takes a list of folders and uses ase gui to plot the structures with helpful titles. Useful for comparing and viewing many bilayers.

## Info
- Before collecting data or computing results, be sure to run vdwcorretion.py on all the monolayers to get the vdw correction. The correction is saved in a separate file, so you may need to add the correction yourself if you are running a custom script.

## Running
Currently being run is
- Everything from C2DB
- From Thomas mats: A2*, A3*, A, AB, AB2, ABC, ABCD, AB3, AB4, AB5 ABC2 AB2C2
