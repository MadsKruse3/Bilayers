from pathlib import Path
from utils import listdirs
import os, re
from asr.core import read_json
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from gpaw import GPAW

def get_IDs(monolayer_folder):
    fname = 'info.json'
    fname_structure = 'structure.json'
    path = "/".join([x for x in monolayer_folder.split("/") if x != ""][-3:])
    c2db_path = "/home/niflheim2/cmr/C2DB-ASR/tree/" + path + "/" + fname
    icsd_cod_path = "/home/niflheim2/cmr/C2DB-ASR/icsd_cod_materials/tree/" + path + "/" + fname
    excluded_AFM_path = "/home/niflheim2/cmr/C2DB-ASR/excluded_AFM_tree/" + path + "/" + fname
    c2db_path_structure = "/home/niflheim2/cmr/C2DB-ASR/tree/" + path + "/" + fname_structure
    icsd_cod_path_structure = "/home/niflheim2/cmr/C2DB-ASR/icsd_cod_materials/tree/" + path + "/" + fname_structure
    excluded_AFM_path_structure = "/home/niflheim2/cmr/C2DB-ASR/excluded_AFM_tree/" + path + "/" + fname_structure
    if Path(c2db_path).is_file():
        #data = read_json(f'{c2db_path}')
        n = len(read(f'{c2db_path_structure}'))
        if n > 8:
            print(n)

    if Path(icsd_cod_path).is_file():
        #data = read_json(f'{icsd_cod_path}')
        n = len(read(f'{icsd_cod_path_structure}'))
        if n > 8:
            print(n)

    if Path(excluded_AFM_path).is_file():
        #data = read_json(f'{excluded_AFM_path}')
        n = len(read(f'{excluded_AFM_path_structure}'))
        if n > 8:
            print(n)

    if not Path(c2db_path).is_file() and not Path(icsd_cod_path).is_file() and not Path(excluded_AFM_path).is_file():
        n = 'missing'
        print('missing n:', monolayer_folder)

    return n         
    


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="Monolayer folders to analyse.")

    args = parser.parse_args()
    
    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]

    for folder in folders:
        if os.path.isfile(f'{folder}/results-asr.interlayer_magnetic_exchange.json'):
           
            monolayer_folder = re.split('[/]',f'{folder}')
            monolayer_folder.remove(monolayer_folder[-1])
            monolayer_folder = "/".join(monolayer_folder)
     
            n = get_IDs(monolayer_folder)
            
