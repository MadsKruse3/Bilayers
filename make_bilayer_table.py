from pathlib import Path
import numpy as np
import os, re, sys, glob, json

from asr.core import read_json
from ase.io import read

from tabulate import tabulate
from texttable import Texttable
import latextable

from bs4 import BeautifulSoup

### call UTILS ######
code_path = '/home/niflheim2/cmr/WIP/stacking/tree-sahar/Workflow_revised2'
sys.path.insert(1,code_path)
from collectlayerdata_revised import get_energy_length
from bilayer_utils import *


def get_descriptor(folder):
    translation = read_json(f'{folder}/translation.json')['translation_vector']
    transform = read_json(f'{folder}/transformdata.json')

    rotation = transform['rotation']

    t_c = transform['translation'][:2] + translation[:2]

    p = "'" if not np.allclose(rotation, np.eye(3)) else ""
    B = 'B' if not np.allclose(t_c, 0.0) else 'A'

    descriptor = 'A' + B + p
    return descriptor

def get_full_descriptor(folder, atoms=None):
    from asr.bilayerdescriptor import get_matrix_descriptor
    if atoms is None:
        if not Path(f"{folder}/structure.json").is_file():
            p = Path(folder).resolve().parents[0]
            atoms = read(f"{p}/structure.json")
        else:
            p = Path(f'{folder}/structure.json').absolute()
            atoms = read(str(p))

    folder = [x for x in folder.split("/") if x != ""][-1]
    folder = folder.replace("--", "-M").replace("_-", "_M")
    desc = "-".join(folder.split("-")[2:])

    # Extract matrix
    def tofloat(x):
        if x.startswith("M"):
            return - float(x[1:])
        else:
            return float(x)

    parts = [p for p in desc.split("-") if p != ""]
    (a, b, c, d) = parts[0].split("_")
    matrix = np.array([[tofloat(a), tofloat(b)], [tofloat(c), tofloat(d)]])
    (tx, ty) = parts[-1].split("_")
    tx = tofloat(tx)
    ty = tofloat(ty)
    iz = "Iz" in desc

    descriptor = get_matrix_descriptor(atoms, matrix)

    if iz:
        descriptor = f'({descriptor}, Iz)'
    else:
        descriptor = f'({descriptor})'

    descriptor += f'_({tx:0.2f}, {ty:0.2f})'

    descriptor = descriptor.replace("_", "  ")
    return descriptor

def get_proper_name(name):
    new_word = ''
    for x in name:
        if x == '_':
            new_word += '$'
            new_word += '\\'
            new_word += '_'
            new_word += '$'
        else:
            new_word += f'{x}'

    return new_word

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="folders to analyse.")

    args = parser.parse_args()
    
    if len(args.folders) > 0:
        mlfolders = [Path(x).absolute() for x in args.folders]
    else:
        mlfolders = [Path(".").absolute()]

    mlfolders = parse_folders()[-1]
    mlfolders = [mlfolders]

    name = re.split('[/]',f'{mlfolders}')[-1]
    name = re.split('[)]', f'{name}')[0]
    name = name[:-1]

    os.chdir('/home/niflheim2/cmr/WIP/stacking/tree-mads/tables/generated_tables')
    with open(f'{name}_table_generated.txt','w') as f1:
        for ml in mlfolders:
            rows = []
            row = []
          
            row.append(name)
            row.append('BL descriptor')
            row.append('Type')
            row.append('Type ref.')
            row.append('Dyn. stab.')
            row.append('zscan')
            row.append('gs')
            row.append('E$_{b}$ ref.')
            row.append('Code')
            row.append('xc')
            row.append('vdW')
            row.append('lat. c.')
            row.append('lat. c. ref.')
            row.append('IL')
            row.append('X-X dist')
            row.append('IL ref.')
            row.append('X-X ref.')
            
            rows.append(row)
         
            bilayers = []       
            for bl in listdirs(ml):
                if '-2-' in bl:
                    if os.path.isfile(f"{bl}/results-asr.relax_bilayer.json"):
                        bilayers.append(bl)
                        stab = read_json(f"{bl}/Stability.json")
                        print(stab)

            zscan_selected = stability_criterion(bilayers, 0.010, source='zscan')
            gs_selected = stability_criterion(zscan_selected, 0.003, source='gs')
            #for bl in zscan_selected:
            for bl in gs_selected:
                row = []
            
                blname = re.split('[/]',f'{bl}')[-1]
                full_descriptor = get_full_descriptor(bl)
                full_descriptor = BeautifulSoup(full_descriptor, features="html5lib").get_text()
                stacking_type = get_descriptor(bl)

                try:
                    with open(f'{bl}/Stability.txt') as f:
                        dyn_stab = f.readlines()[0]
                except:
                    dyn_stab = 'Not Available'

                blname = get_proper_name(blname)
                structure = read(f'{bl}/structure.json')
                
                cell = structure.get_cell()
                a = cell[0]
                b = cell[1]
                a = np.linalg.norm(a)
                b = np.linalg.norm(b)

                #print('Prefactor:', np.linalg.norm(np.cross(cell[0], cell[1]) )/(np.linalg.norm(cell[0])*np.linalg.norm(cell[1]) ) )
                #print('Area:', np.linalg.det(2*cell[:2,:2])

                eb_zscan, length_closest = get_energy_length(bl, 'zscan')
                try: 
                    eb_gs, _ = get_energy_length(bl, 'gs')
                except:
                    eb_gs = 'Not Available'
               
                if os.path.isfile(f"{bl}/results-asr.zscan.json"):
                    length_sameatom = read_json(f"{bl}/results-asr.zscan.json")["optimal_height"]
                else:
                    length_sameatom = read_json(f"{bl}/results-asr.relax_bilayer.json")["optimal_height"]    #Magnetic: zscn.json  


                row.append(blname) #bilayer folder
                row.append(full_descriptor) #full descriptor
                row.append(stacking_type) # Our type
                row.append('\\textcolor{blue}{ref}') #Their type
                row.append(dyn_stab) #Dynamical stability 
                row.append(eb_zscan*1000) #binding energy with zscan recipe
                row.append(eb_gs*1000) #binding energy with gs recipe
                row.append('\\textcolor{blue}{ref}') # Their binding energy 
                row.append('\\textcolor{blue}{code}') #DFT code used 
                row.append('\\textcolor{blue}{xc}') #XC functional used
                row.append('\\textcolor{blue}{vdw}') #vdW correction used
                row.append((str("{:#.3g}".format(a)).rstrip('.'),str("{:#.3g}".format(b)).rstrip('.'))) # Our lattice constant
                row.append(('-','-')) ## Their lattice constants
                row.append(length_sameatom) # Interlayer distance 
                row.append(length_closest) # Distance between closest atoms
                row.append('\\textcolor{blue}{ref}') # Interlayer distance ref.
                row.append('\\textcolor{blue}{ref}') # Distance between closest atoms ref.

                rows.append(row)

        sys.stdout = f1
        table = Texttable()
        table.set_cols_align(["c"] * 17)
        table.set_deco(Texttable.HEADER | Texttable.VLINES)
        format_array = ["t", "t", "t", "t", "t", "f", "f", "t", "t", "t", "t", "f", "t", "f", "f", "t", "t"] 
        table.set_cols_dtype(format_array)
        table.set_precision(2)
        table.add_rows(rows)
        print(latextable.draw_latex(table))
        print(latextable.draw_latex(table, caption=f"""The table contains comparisons between our results for {name} 
        with those found elsewhere in literature. The columns contain: bilayer stackings, descriptors, 
        stacking type, stacking type of the reference paper, dynamical stability, binding energy using the zscan recipe,
        binding energy using the gs recipe, binding energy found in the reference paper, DFT code used in the reference paper, 
        XC-functional used in the reference paper, vdW-correction (if any) used in the reference paper, 
        lattice constants, lattice constants found in the reference paper,
        interlayer distances measured between the same atoms in the middle of each layer,
        the shortest interlayer distance measured between the highest atom in the bottom layer 
        and the lowest hanging atom in the top layer and the same two measures (if computed) found in the reference paper."""))
    
    f1.close()
